"""
ParaRead is a module for processing sequencing reads (bam or sam files) in
parallel. You use it by writing a child class that extends the included
ParaReadProcessor class. It should implement the __call__() function at a minimum,
defining what to do with each chromosome.

Then you simply call object.registerFiles(),
then object.runProcess(), object.combine().
"""

import pysam
import multiprocessing
import sys
from collections import OrderedDict
import shutil # to remove temporary folder
import subprocess # to run shell command combining results
import os
import tempfile # to create temporary directories
import atexit

"""
This module-level variable acts like a global. The ParaReadProcessor child
classes will be function classes that will be called by the
multiprocessing pool map function, which must pickle everything in the call.
To work around this, I stick any unpicklable or large items into this paraReadFiles
dict, so they aren't part of the class, and thus will not cause map to fail.
"""
paraReadFiles = {}


class ParaReadProcessor(object):
	"""
	An abstract class. If writing a child that inherits from this class, you
	should include:

	self.resultAcronym -- a short identifier summarizing the type of result
	produced by your class. It will be appended to temp files, default outnames,
	etc. __call__() -- a function that takes a chrom, and writes results to the
	tempFolder/chrom.tsv file registerFiles() -- a function that adds elements to
	the variable ParaRead.paraReadFiles for anything that is required by your
	function.
	"""

	def __init__(self, bamFileName, nProc, outFile, tempParent, limit=None, allow_unaligned=True):
		"""
		@param limit Limit to only certain chromosomes (given as a list)
		"""
		if not hasattr(self, 'resultAcronym'):
			raise NotImplementedError("The subclass should assign resultAcronym before calling the ParaReadProcessor superclass constructor.")
		self.bamFileName = bamFileName
		self.nProc=int(nProc)
		self.setOutFile(bamFileName, outFile)
		self.makeTempFolder(bamFileName, tempParent)
		self.limit=limit
		self.allow_unaligned = allow_unaligned

	def setOutFile(self, bamFileName, outFile):
		"""
		Sets the outFile, creating something sensible from the bamFileName
		if nothing is provided.
		"""
		if outFile == "" or outFile is None:
			outFile = os.path.splitext(os.path.basename(bamFileName))[0]
			outFile += "_" + self.resultAcronym + ".txt"
		self.outFile = outFile

	def makeTempFolder(self, bamFileName, tempParent):
		"""
		Given a 'parent' directory (CWD by default), this generates a
		randomized folder name and creates it for temporary files. The name of
		the folder includes the bamFileName plus random characters.
		"""
		if tempParent =="" or tempParent is None:
			tempParent=os.path.dirname(self.outFile)

		prefix = "tmp_" + self.resultAcronym + "_"
		prefix += os.path.splitext(os.path.basename(bamFileName))[0] + "_"
		self.tempFolder = tempfile.mkdtemp(prefix=prefix, dir=tempParent)

		# when we make the temp folder, we should also register to delete it
		# upon termination.
		atexit.register(self.exit_handler)

	def exit_handler(self):
		"""
		Clean up properly
		"""
		shutil.rmtree(self.tempFolder)

	def __call__(self, chrom):
		"""
		The SubClass the inherits ParaReadProcessor must implement __call__.
		Your __call__ function should take a chrom string, and process all the
		reads on that chromosome, producing a file called tempFolder/chrom.txt
		and then return either chrom, if the process succeeded, or 'None', if
		it failed. It may also call a subfunction.
		"""
		raise NotImplementedError( "Child class must implement __call__" )

	def registerFiles(self, bamFormat=None):
		"""
		registerFiles() should add to the paraReadFiles dict any variables that
		will be required by your __call__ function, but are too big or unpicklable,
		and thus cannot be included in the class itself.
		It can be overwritten by the child class, but if so, you
		should use super to call this parent version as well, to populate the
		variable with bamFile.
		"""

		print("Using input file:" + self.bamFileName)
		# pysam can read either sam or bam input, but you must specify which
		# I can try to guess based on extension, though.
		if (bamFormat is None):
			if self.bamFileName.endswith(".bam"):
				bamFormat = True
			elif self.bamFilename.endswith(".sam"):
				bamFormat = False
			else:
				raise NotImplementedError("ParaRead requires sam or bam input")

		if (bamFormat is True):
			mode="rb"
		else:
			mode="r"

		# Here, check_sq is necessary so that ParaRead can process unaligned files,
		# which is occasionally desirable.
		paraReadFiles['bamFile'] = pysam.AlignmentFile(self.bamFileName , mode=mode, check_sq = False)

	def unbufferedWrite(self, txt):
		""" Writes unbuffered output by flushing after each stdout.write call """
		sys.stdout.write(txt)
		sys.stdout.flush()

	def getChromListFromBamHeader(self, bamFile):
		"""Get a list of chromosomes (and lengths) in this bamfile from header"""
		# chrlist = [item['SN'] for item in bamFile.header['SQ']]
		# Dict doesn't preserve order:
		# chrlist = {item['SN']: item['LN'] for item in bamFile.header['SQ']}
		# it's important to preserve order, since the chr ID in individual reads
		# indexes into this chromosome order.
		chrlist=OrderedDict()
		for item in bamFile.header['SQ']:
			if self.limit is not None:
				if item['SN'] in self.limit:
					chrlist[item['SN']] = item['LN']
			else:
				chrlist[item['SN']] = item['LN']

		print("Chrom limit: " + str(self.limit))

		#if (len(chrlist) == 0):
		#	raise NotImplementedError("Is this file aligned?")
		#chrlist.keys() # chromosome list
		return(dict(chrlist))

	def chrlist(self):
		#global paraReadFiles
		return self.getChromListFromBamHeader(paraReadFiles['bamFile'])

	def runProcess(self):
		"""
		Function actually runs the parallel processing, calling your function
		on each chromosome independently
		"""
		# Because this class is a function class (implemenets __call__), I can
		# call "self()" as a function, which is what runs the match function
		# on a single chromosome. By mapping self() across multiple chroms,
		# I get the parallel version.

		print("Temporary files stored in: " + self.tempFolder)
		print("Cores used: " + str(self.nProc))
		self.unbufferedWrite("Processing...")
		chrlist = self.chrlist()
		self.chrs = chrlist.keys() #e.g. ['chr1', 'chr2', 'chr3', 'chr4']

		# Unaligned files lack any chrom header lines; we just pass "all"
		# and at the moment, can't process these in parallel (but this could
		# change if someone wants to implement it).
		if (len(self.chrs) == 0):
			if self.allow_unaligned:
				self.chrs=["all"]
			else:
				raise Exception("No chroms in header;  This file is empty or unaligned. Aligned reads are required.")
		# Process chromosomes in parallel

		if self.nProc == 1:
			r = map(self, self.chrs)
		else:
			p = multiprocessing.Pool(self.nProc)
			# The typical call to map fails to acknowledge KeyboardInterrupts,
			# this fix, http://stackoverflow.com/a/1408476/946721, seems to help
			#r = p.map(self, self.chrs)
			r = p.map_async(self, self.chrs).get(9999999)

		self.goodChrs = [x for x in r if x is not None]
		self.badChrs = [self.chrs[i] for i, x in enumerate(r) if x == None]
		print ("\nFailed chroms: " + str(self.badChrs))
		print ("\nSuccessful chroms: " + str(self.goodChrs))

	def combine(self):
		"""
		After running the process in parallel, this 'reduce' step will simply
		merge all the temporary files into one, and rename it to the output
		file name.
		"""
		# Build a command to concat all the output files.
		if (len(self.goodChrs) < 1):
			print("No successful chromosomes, so no combining.")
		else:
			cmd = "cat " + \
				" ".join(
				[os.path.join(self.tempFolder, chr+".txt") for chr in self.goodChrs]) + \
				" > " + self.outFile
		subprocess.call(cmd, shell=True)

		print("Outfile: " + self.outFile)

	def reverse_complement(self, sequence):
		# Credits: crazyhottommy
		"""
		Returns reverse complement of a DNA string.
		"""
		seq_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
		return "".join([seq_dict[base] for base in reversed(sequence)])

	def complement(self, sequence):
		"""
		Returns the complement of a DNA sequence without reversing the string.
		"""
		seq_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
		return "".join([seq_dict[base] for base in sequence])
