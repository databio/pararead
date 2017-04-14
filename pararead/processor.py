"""
Process sequencing reads in SAM or BAM format in parallel.

Do so by deriving from ParaReadProcessor, critically implementing __call__().
That defining what to do with each chromosome.
Then you simply call object.register_files(),
then object.run(), object.combine().
"""

import abc
import atexit
import itertools
import multiprocessing
import os
import shutil
import tempfile

import pysam

from .exceptions import FileTypeException
from .utils import *


"""
This module-level variable acts like a global.
The ParaReadProcessor child classes will be
function classes that will be called by the
multiprocessing pool map function, which must
pickle everything in the call. To work around this,
I stick any unpicklable or large items into this
PARA_READ_FILES dict, so they aren't part of the class,
and thus will not cause map to fail.
"""
PARA_READ_FILES = {}
READS_FILE_KEY = "readsfile"



class ParaReadProcessor(object):
    """
    An abstract class for parallel processing of sequencing reads.
    Implement __call__ to define what will be done with each batch
    (e.g., chromosome) of reads.
    """

    __metaclass__ = abc.ABCMeta

    SUPPORTED_FILETYPES = ["BAM", "SAM"]     # Restrict input types.


    def __init__(self,
                 path_reads_file, cores, chunksize,
                 outfile=None, action=None,
                 temp_folder_parent_path=None,
                 limit=None, allow_unaligned=True,
                 require_new_outfile=False, output_type="txt"):
        """
        Regardless of subclass behavior, there is a
        set of fields that an instance should have.

        Parameters
        ----------
        path_reads_file : str
            Data location (aligned BAM/SAM file).
        cores : int | str
            Number of processors to use.
        chunksize : int
            Number of reads per processing chunk.
        outfile : str, optional
            Path to location for output file. Either this
            or `action` is required.
        action : str, optional
            What the child class is doing, used to
            derive outfile name if unspecified. If
            `outfile` is unspecified, this is required.
        temp_folder_parent_path : str, optional
            Temporary folder. If unspecified, this will
            match the folder containing the output file.
        limit : list of str, optional
            Which chromosomes to process, process all by default.
        allow_unaligned : bool, default True
            Whether to allow unaligned reads.
        require_new_outfile : bool
            Whether to raise an exception if output file already exists.
        output_type : str, optional
            Type of output file(s) generated. This is used by both 
            intermediate files that are created and by the combine() 
            step that creates final output.

        Raises
        ------
        ValueError
            If given neither `outfile` path nor `action` action name,
            or if output file already exists and a new one is required.

        """

        # Initial path setup and filetype handling.
        name_reads_file = os.path.basename(path_reads_file)
        readsfile_basename, _ = os.path.splitext(name_reads_file)
        self.path_reads_file = path_reads_file

        # Use given output path or derive it from processing behavior.
        if outfile:
            self.outfile = outfile
        elif action:
            self.outfile = make_outfile_name(readsfile_basename, action)
        else:
            raise ValueError("Either path to output file or "
                             "name of processing action is required.")

        if os.path.exists(self.outfile):
            if require_new_outfile:
                raise ValueError(
                        "Outfile already exists: '{}'".format(self.outfile))
            else:
                print("WARNING: Output file already exists and "
                      "will be overwritten: '{}'".format(self.outfile))

        # Create temp folder that's deleted upon exit.
        if not temp_folder_parent_path:
            temp_folder_parent_path = os.path.dirname(self.outfile)

        # Handling of temporary folder.
        prefix = "tmp_{}_".format(readsfile_basename)
        if action:
            prefix += "{}_".format(action)
        self.temp_folder = tempfile.mkdtemp(
                prefix=prefix, dir=temp_folder_parent_path)

        # Add a couple lines so that tests can execute quietly.
        # Tests handle cleanup separately, so if this existence
        # check is omitted, exceptions can be squelched, but their
        # messages can still appear due to rmtree()'s use of sys.exc_info().
        def clean():
            if os.path.exists(self.temp_folder):
                shutil.rmtree(self.temp_folder)
        atexit.register(clean)

        # Behavior/execution parameters.
        self.cores = int(cores)
        self.chunksize = chunksize
        self.limit = limit
        self.allow_unaligned = allow_unaligned
        self.output_type = output_type


    @abc.abstractmethod
    def __call__(self, reads_chunk, chunk_id):
        """
        Perform 'chunk'-wise processing implemented in the subclass.
        
        A concrete implementation operates on a chunk of sequencing reads. 
        By default, a 'chunk' simply consists of contiguous reads from the 
        input file, but a subclass that overrides the partitioning mechanism 
        can change that. In that case, the provided 'chunk_id' will likely 
        be meaningful, but otherwise it can be ignored. If processing fails 
        for the given chunk, the concrete implementation should communicate 
        this by returning None.

        Parameters
        ----------
        reads_chunk : Iterable of pysam.AlignedSegment
            Chunk of sequencing reads to process.
        chunk_id : int or str
            Reads chunk identifier; this is passed in case it's of use 
            in the client implementation, but for many use cases it 
            can be ignored.
        
        Returns
        -------
        None or object
            None if the chunk's processing failed, otherwise a non-null result.
        
        """
        pass


    def _tempf(self, chrom):
        """
        Derive name for temporary file from chromosome name.

        Parameters
        ----------
        chrom : str
            Name of chromosome (single processing partition).

        Returns
        -------
        str
            Name for tempfile corresponding to given unit name.

        """
        return os.path.join(self.temp_folder,
                            "{}.{}".format(chrom, self.output_type))


    def register_files(self):
        """
        Add to module map any large/unpicklable variables required by __call__.
        
        This can be overridden by the child class, but if so you should call
        this version as well, to populate the mapping with needed file(s).
        
        Raises
        ------
        FileTypeException
            If the path to the reads file given doesn't appear 
            to match one of the supported file types.
        """

        print("Using input file: '{}'", self.path_reads_file)
        _, extension = os.path.splitext(self.path_reads_file)
        filetype = extension[1:].upper()

        if filetype not in self.SUPPORTED_FILETYPES:
            raise FileTypeException(self.path_reads_file,
                                    self.SUPPORTED_FILETYPES)

        # TODO: pysam docs say 'u' for uncompressed BAM.
        mode = 'rb' if filetype == 'BAM' else 'r'

        # Here, check_sq is necessary so that ParaRead can process
        # unaligned files, which is occasionally desirable.
        readsfile = pysam.AlignmentFile(
                self.path_reads_file, mode=mode, check_sq=False)
        PARA_READ_FILES[READS_FILE_KEY] = readsfile

        def ensure_closed():
            if readsfile.is_open:
                readsfile.close()
        atexit.register(ensure_closed)


    def partition(self, readsfile):
        """
        Partition sequencing reads into equally-sized 'chunks' for 
        parallel processing. This treats all reads equally, in that 
        they are grouped by contiguous occurrence within the given 
        input file. This facilitates processing of unaligned reads, 
        but it means that reads from the same chromosome will not be 
        processed together. This can be overridden if that's desired.
        
        Parameters
        ----------
        readsfile

        Returns
        -------

        """
        return itertools.groupby(enumerate(readsfile),
                                 key=lambda (i, r): i // self.chunksize)


    def run(self):
        """
        Do the processing defined partitioned across each unit (chromosome).

        Returns
        -------
        collections.Iterable of str
            Names of chromosomes for which result is non-null.

        Raises
        ------
        MissingHeaderException
            If attempting to run with an unaligned reads file 
            in the context of an aligned file requirement.

        """
        # Because this class is a function class (implements __call__), I can
        # call "self()" as a function, which is what runs the match function
        # on a single chromosome. By mapping self() across multiple chroms,
        # I get the parallel version.

        print("Temporary files will be stored in: '{}'".
              format(self.temp_folder))
        print("Processing with {} cores...".format(self.cores))

        # e.g. ['chr1', 'chr2', 'chr3', 'chr4']
        try:
            chunks = self.partition(PARA_READ_FILES[READS_FILE_KEY])
        except KeyError:
            print("No '{}' has been established; call "
                  "'register_files' before 'run'".format(READS_FILE_KEY))
            raise

        """
        # Unaligned files lack any chrom header lines.
        # If so, we'll get back a null to disambiguate empty filter case if
        # we're permitting unaligned input. If requiring aligned input, the
        # call will have already generated an exception to that effect.
        if chroms is None:
            # If permitting unaligned input, set chroms to an indicator.
            # At the moment, such a case won't parallelize,
            # but this could change if someone wants to implement it.
            raise MissingHeaderException(
                    filepath=PARA_READ_FILES[READS_FILE_KEY].filename)
        if len(chroms) == 0:
            # This is the case in which filtration left us with no chroms.
            print("No chromosomes retrieved from '{}' header when "
                  "filtering to: {}".format(self.path_reads_file,
                                            self.limit))
            return []
        """

        # Some implementors may have a strand mode attribute.
        # If so, log it here to avoid duplicate messaging, as it
        # will remain constant across processed chunks (chromosomes).
        try:
            print("STRAND MODE: {}".format(self.use_strand))
        except AttributeError:
            pass

        chunk_ids, chunks = zip(*chunks)

        # Process chromosomes in parallel
        if self.cores == 1:
            results = map(self, chunks)
        else:
            p = multiprocessing.Pool(self.cores)
            # The typical call to map fails to acknowledge KeyboardInterrupts.
            # This fix helps: http://stackoverflow.com/a/1408476/946721
            results = p.map_async(self, chunks).get(9999999)

        result_by_chunk = zip(chunk_ids, results)
        bad_chunks, good_chunks = \
                partition_chunks_by_null_result(result_by_chunk)

        """
        result_by_chromosome = zip(chroms, results)
        bad_chroms, good_chroms = \
                partition_chunks_by_null_result(result_by_chromosome)
        """

        print("Discarding {} chunk(s) of reads: {}".
              format(len(bad_chunks), bad_chunks))
        print("Keeping {} chunk(s) of reads: {}".
              format(len(good_chunks), good_chunks))

        return good_chunks


    def combine(self, good_chromosomes):
        """
        After running the process in parallel, this 'reduce' step will simply
        merge all the temporary files into one, and rename it to the output
        file name.
        """
        if not good_chromosomes:
            print("No successful chromosomes, so no combining.")
            return
        else:
            print("Merging {} files into output file: '{}'".
                  format(len(good_chromosomes), self.outfile))
            with open(self.outfile, 'w') as outfile:
                for chrom in good_chromosomes:
                    with open(self._tempf(chrom), 'r') as tmpf:
                        for line in tmpf:
                            outfile.write(line)
