"""
Process sequencing reads in SAM or BAM format in parallel.

Do so by deriving from ParaReadProcessor, critically implementing __call__().
That defining what to do with each chromosome.
Then you simply call object.register_files(),
then object.run(), object.combine().
"""

import abc
import atexit
import logging
import multiprocessing
import os
import shutil
import tempfile

import pysam

from .exceptions import *
from .utils import *


_LOGGER = logging.getLogger(__name__)


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

    SUPPORTED_FILETYPES = ["BAM", "SAM"]

    def __init__(self,
                 path_reads_file, cores,
                 outfile=None, action=None,
                 temp_folder_parent_path=None,
                 limit=None, allow_unaligned=True,
                 require_new_outfile=False):
        """
        Regardless of subclass behavior, there is a
        set of fields that an instance should have.

        Parameters
        ----------
        path_reads_file : str
            Data location (aligned BAM/SAM file).
        cores : int | str
            Number of processors to use.
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

        Raises
        ------
        FileTypeException
            If the path to the reads file given doesn't appear 
            to match one of the supported file types.
        ValueError
            If given neither `outfile` path nor `action` action name,
            or if output file already exists and a new one is required.

        """

        # Initial path setup and filetype handling.
        name_reads_file = os.path.basename(path_reads_file)
        readsfile_basename, filetype = os.path.splitext(name_reads_file)
        self.filetype = filetype.upper().strip(".")
        if self.filetype not in self.SUPPORTED_FILETYPES:
            raise FileTypeException(got=self.filetype,
                                    known=self.SUPPORTED_FILETYPES)
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
                _LOGGER.warn("Outfile already exists and will be overwritten: "
                             "'%s'", self.outfile)

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
        self.limit = limit
        self.allow_unaligned = allow_unaligned


    @abc.abstractmethod
    def __call__(self, chrom):
        """
        Perform the processing desired by the base class, per-chromosome.
        This should take a chrom string, and process all the reads on
        that chromosome, producing a file called temp_folder/<chrom>.<ext>
        and then return either chrom, if the process succeeded,
        or 'None', if it failed.

        Parameters
        ----------
        chrom : str
            Chromosome to process
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
        return os.path.join(self.temp_folder, "{}.txt".format(chrom))


    def register_files(self):
        """
        register_files() should add to the PARA_READ_FILES dict any variables
        that are required by a child class's __call__ function, but are too
        big or unpicklable and thus cannot be included in the class itself.
        It can be overridden by the child class, but if so you should call
        this parent version as well, to populate the variable with readsfile.
        """
        path_reads_file = self.path_reads_file
        _LOGGER.info("Using input file: '%s'", path_reads_file)
        mode = 'rb' if self.filetype == 'BAM' else 'r'
        # Here, check_sq is necessary so that ParaRead can process
        # unaligned files, which is occasionally desirable.
        PARA_READ_FILES[READS_FILE_KEY] = pysam.AlignmentFile(
                path_reads_file, mode=mode, check_sq=False)


    def run(self):
        """
        Do the processing defined partitioned across each unit (chromosome).

        Returns
        -------
        collections.Iterable of str
            Names of chromosomes for which result is non-null.

        Raises
        ------

        """
        # Because this class is a function class (implements __call__), I can
        # call "self()" as a function, which is what runs the match function
        # on a single chromosome. By mapping self() across multiple chroms,
        # I get the parallel version.

        _LOGGER.info("Temporary files will be stored in: '{}'".
                     format(self.temp_folder))
        _LOGGER.info("Cores: {}".format(str(self.cores)))
        _LOGGER.info("Processing...")

        # e.g. ['chr1', 'chr2', 'chr3', 'chr4']
        chroms = chromosomes_from_bam_header(
            chroms=self.limit, readsfile=PARA_READ_FILES[READS_FILE_KEY])

        # Unaligned files lack any chrom header lines; we just pass "all"
        # and at the moment, can't process these in parallel (but this could
        # change if someone wants to implement it).
        if len(chroms) == 0:
            if self.allow_unaligned:
                chroms = ["all"]
            else:
                raise MissingHeaderException(PARA_READ_FILES[READS_FILE_KEY])

        # Some implementors may have a strand mode attribute.
        # If so, log it here to avoid duplicate messaging, as it
        # will remain constant across processed chunks (chromosomes).
        try:
            _LOGGER.info("STRAND MODE: {}".format(self.use_strand))
        except AttributeError:
            pass

        # Process chromosomes in parallel
        if self.cores == 1:
            results = map(self, chroms)
        else:
            p = multiprocessing.Pool(self.cores)
            # The typical call to map fails to acknowledge KeyboardInterrupts.
            # This fix helps: http://stackoverflow.com/a/1408476/946721
            results = p.map_async(self, chroms).get(9999999)

        result_by_chromosome = zip(chroms, results)
        bad_chroms, good_chroms = self.filter_chroms(result_by_chromosome)

        _LOGGER.info("Discarding {} chromosomes: {}".
                     format(len(bad_chroms), bad_chroms))
        _LOGGER.info("Keeping {} chromosomes: {}".
                     format(len(good_chroms), good_chroms))

        return good_chroms


    @staticmethod
    def filter_chroms(result_by_chromosome):
        """
        Bin chromosome name by whether processing result was null.

        Parameters
        ----------
        result_by_chromosome : Sequence of (str, object)
            Mapping from chromosome name to processing result

        Returns
        -------
        (list of str, list of str)
            Sequence of names of chromosomes for which result was null,
            and an analogous sequence for those with a null result.

        """
        bad_chroms, good_chroms = [], []
        for c, r in result_by_chromosome:
            chroms = bad_chroms if r is None else good_chroms
            chroms.append(c)
        return bad_chroms, good_chroms


    def combine(self, good_chromosomes):
        """
        After running the process in parallel, this 'reduce' step will simply
        merge all the temporary files into one, and rename it to the output
        file name.
        """
        if not good_chromosomes:
            _LOGGER.info("No successful chromosomes, so no combining.")
            return
        else:
            _LOGGER.debug("Merging %d files into output file: '%s'",
                          len(good_chromosomes), self.outfile)
            with open(self.outfile, 'w') as outfile:
                for chrom in good_chromosomes:
                    _LOGGER.debug("Fetching results: %s", chrom)
                    with open(self._tempf(chrom), 'r') as tmpf:
                        _LOGGER.debug("Adding contents from: '%s'", tmpf.name)
                        for line in tmpf:
                            outfile.write(line)