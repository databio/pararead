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
import logging
import multiprocessing
import os
import shutil
import tempfile

from ubiquerg import is_command_callable
from .exceptions import \
    CommandOrderException, IllegalChunkException, \
    MissingOutputFileException, UnknownChromosomeException
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
CHUNKS_PER_CORE = 5
CORES_PARAM_NAME = "cores"


_LOGGER = logging.getLogger(__name__)


class ParaReadProcessor(object):
    """
    Base class for parallel processing of sequencing reads.
    
    Implement __call__ to define work for each reads chunk, (e.g., chromosome).
    Unaligned reads are permitted, but the work then cannot rely on any sort 
    of biologically meaningful chunking of the reads unless a partition() 
    function is implemented. If unaligned reads are used and no partition() is 
    implemented, reads will be arbitrarily split into chunks.
    
    """

    __metaclass__ = abc.ABCMeta

    def __init__(
            self, path_reads_file, cores, outfile=None, action=None,
            temp_folder_parent_path=None, limit=None, allow_unaligned=False,
            require_new_outfile=False, by_chromosome=True,
            intermediate_output_type="txt", output_type="txt",
            retain_temp=False):
        """
        :param str path_reads_file: data location (aligned BAM/SAM file).
        :param int | str cores: number of processors to use.
        :param str outfile: path to location for output file; either this
            or action is required.
        :param str action: name for what the child class is doing, used to
            derive outfile name if unspecified; if outfile is unspecified,
            this is required.
        :param str temp_folder_parent_path: temporary folder; if unspecified,
            this will match the folder containing the output file.
        :param list[str] limit: which chromosomes to process, process all by
            default.
        :param bool allow_unaligned: whether to allow unaligned reads.
        :param bool require_new_outfile: whether to raise an exception if
            output file already exists.
        :param bool by_chromosome: whether to chunk reads on a per-chromosome
            basis, implicitly imposing requirement for aligned reads.
        :param str intermediate_output_type: type of output file generated for
            each chunk of reads processed.
        :param str output_type: type of final output file generated; this is
            used by both intermediate files that are created and by the combine()
            step that creates final output.
        :raise ValueError: if given neither `outfile` path nor `action` action
            name, or if output file already exists and a new one is required.
        """

        # # Establish root logger only if client application hasn't done so.
        # # That is, create a root logger with a handler if one doesn't exist.
        # if not logging.getLogger().handlers:
        #     global _LOGGER
        #     _LOGGER = setup_logger(make_root=True, propagate=False)

        # Initial path setup and filetype handling.
        name_reads_file = os.path.basename(path_reads_file)
        readsfile_basename, _ = os.path.splitext(name_reads_file)
        self.path_reads_file = path_reads_file

        # Use given output path or derive it from processing behavior.
        if outfile:
            self.outfile = outfile
        elif action:
            self.outfile = make_outfile_name(
                    readsfile_basename, action, output_type)
        else:
            raise ValueError("Either path to output file or "
                             "name of processing action is required.")

        # Perform check after establishing the setting so that it works
        # regardless of whether the path was explicit or inferred.
        if os.path.exists(self.outfile):
            if require_new_outfile:
                raise ValueError("Output file already exists: '{}'".
                                 format(self.outfile))
            else:
                _LOGGER.warning(
                        "Output file already exists and "
                        "will be overwritten: '{}'".format(self.outfile))

        # Create temp folder that's deleted upon exit.
        if not temp_folder_parent_path:
            temp_folder_parent_path = os.path.dirname(self.outfile)

        # Handling of temporary folder.
        prefix = "tmp_{}_".format(readsfile_basename)
        if action:
            prefix += "{}_".format(action)
        tempfolder = tempfile.mkdtemp(
                prefix=prefix, dir=temp_folder_parent_path)
        self.temp_folder = tempfolder

        # Add a couple lines so that tests can execute quietly.
        # Tests handle cleanup separately, so if this existence
        # check is omitted, exceptions can be squelched, but their
        # messages can still appear due to rmtree()'s use of sys.exc_info().
        def clean():
            if os.path.exists(tempfolder):
                shutil.rmtree(tempfolder)
        
        if not retain_temp:
            atexit.register(clean)

        # Behavior/execution parameters.
        self.cores = int(cores)
        self.limit = limit
        self.require_aligned = by_chromosome or not allow_unaligned
        self.intermediate_output_type = intermediate_output_type
        self.by_chromosome = by_chromosome
        self._size_by_chromosome = None

    @abc.abstractmethod
    def __call__(self, chunk_id, reads_chunk):
        """
        Perform 'chunk'-wise processing implemented in the subclass.

        A concrete implementation operates on a chunk of sequencing reads. 
        By default, a 'chunk' simply consists of contiguous reads from the 
        input file, but a subclass that overrides the partitioning mechanism 
        can change that. In that case, the provided 'chunk_id' will likely 
        be meaningful, but otherwise it can be ignored. If processing fails 
        for the given chunk, the concrete implementation should communicate 
        this by returning None.

        :param int | str chunk_id: reads chunk identifier; this is passed in
            case it's of use in the client implementation, but for many use
            cases it can be ignored.
        :param Iterable reads_chunk: chunk of sequencing reads to process,
            likely pysam.AlignedSegment
        :return NoneType | object: null if the chunk's processing failed,
            otherwise a non-null result.
        """
        pass

    @property
    def files(self):
        """
        Refer to the pararead files mapping.
        
        :return Mapping[str, object]: pararead files mapping.
        """
        return PARA_READ_FILES

    @property
    def readsfile(self):
        """
        :return pysam.AlignmentFile | pysam.VariantFile: instance of the reads
            file abstraction appropriate for the given type of input data
            (e.g., BAM or VCF).
        :raise pararead.exceptions.CommandOrderException: if a command
            prerequisite for a parallel reads processor operation has not yet
            been performed.
        """
        return self.fetch_file(READS_FILE_KEY)


    @staticmethod
    def empty_action(read_chunk_key=None):
        """
        Action to take when processing an empty reads chunk.

        :param str read_chunk_key: key for the empty chunk of reads.
        """
        if read_chunk_key:
            _LOGGER.debug("Empty read chunk: {}".format(read_chunk_key))

    def fetch_file(self, file_key):
        """
        Retrieve one of the files registered with pararead.

        :param str file_key: which file to fetch
        :return object: likely pysam.AlignmentFile -- file ADT instance
            associated with the requested key.
        :raise pararead.exceptions.CommandOrderException: if the indicated file
            hasn't been registered.
        """
        try:
            return self.files[file_key]
        except KeyError:
            raise CommandOrderException(
                "No {} established; has {} been called?".format(
                    READS_FILE_KEY, ParaReadProcessor.register_files.__name__))

    def check_command(self, cmd):
        """
        Determine whether it appears that a command may be run.

        :param str cmd: command to check for runnability
        :return OSError: if it's possible to verify that running given command
            would fail
        """
        if not is_command_callable(cmd):
            raise OSError("{} is not callable".format(cmd))

    def get_chrom_size(self, chrom):
        """
        Determine the size of the given chromosome.
        
        :param str chrom: name of chromosome of interest.
        :return int: size of chromosome of interest.
        :raise pararead.exceptions.CommandOrderException: if there's no
            chromosome sizes map yet.
        :raise pararead.exceptions.UnknownChromosomeException: if requested
            chromosome is not in the sizes map.
        """
        if not self._size_by_chromosome:
            raise CommandOrderException(
                    "No size-by-chromosome mapping; "
                    "has a reads file been registered?")
        try:
            return self._size_by_chromosome[chrom]
        except KeyError:
            raise UnknownChromosomeException(
                    chrom, known=self._size_by_chromosome.keys())


    def register_files(self, **file_builder_kwargs):
        """
        Add to module map any large/unpicklable variables required by __call__.

        :raise pararead.exceptions.FileTypeException: if path to the reads file
            given doesn't appear to match one of the supported file types.
        """

        _LOGGER.info("Registering input file: '%s'", self.path_reads_file)
        reads_file_maker = create_reads_builder(self.path_reads_file)

        # Here, check_sq is necessary so that ParaRead can process
        # unaligned files, which is occasionally desirable.
        builder = reads_file_maker.ctor
        kwargs = file_builder_kwargs

        # Builder has minimal requirements that must be met.
        # Thus, those take precedence over user provisions.
        kwargs.update(reads_file_maker.kwargs)

        if not self.require_aligned:
            kwargs['check_sq'] = False

        readsfile = builder(self.path_reads_file, **kwargs)
        PARA_READ_FILES[READS_FILE_KEY] = readsfile

        # Cache mapping from chromosome name to size for easy access.
        self._size_by_chromosome = parse_bam_header(
                readsfile, require_aligned=self.require_aligned)

        def ensure_closed():
            if readsfile.is_open:
                readsfile.close()
        atexit.register(ensure_closed)

    def run(self, chunksize=None, interleave_chunk_sizes=False):
        """
        Do the processing defined partitioned across each unit (chromosome).

        :param int chunksize: number of reads per processing chunk; if
            unspecified, the default heuristic of size s.t. each core gets ~ 4
            chunks.
        :param bool interleave_chunk_sizes: whether to interleave reads chunk
            sizes. If off (default), just use the distribution that Python
            determines.
        :return Iterable[str]: names of chromosomes for which result is non-null.
        :raise pararead.exception.MissingHeaderException: if attempting to run
            with an unaligned reads file in the context of an aligned file
            requirement.
        """
        # Because this class is a function class (implements __call__), I can
        # call "self()" as a function, which is what runs the match function
        # on a single chromosome. By mapping self() across multiple chroms,
        # I get the parallel version.

        try:
            readsfile = PARA_READ_FILES[READS_FILE_KEY]
        except KeyError:
            _LOGGER.error(
                    "No '{}' has been established; call 'register_files' "
                    "before 'run'".format(READS_FILE_KEY))
            raise

        if not self.by_chromosome:
            read_chunk_keys = self.chunk_reads(readsfile, chunksize=chunksize)
        else:
            size_by_chromosome = parse_bam_header(
                readsfile=readsfile, chroms=self.limit,
                require_aligned=self.require_aligned)
            if self.cores == 1:
                # TODO: handle case (unaligned input) of null return.
                # TODO: pysam's fetch() may make this OK but not distribute.
                read_chunk_keys = size_by_chromosome.keys()
            else:
                if size_by_chromosome is None:
                    # Unaligned files lack any chrom header lines.
                    # If so, we'll get back a null to disambiguate
                    # empty filter case if permitting unaligned
                    # input. If requiring aligned input, the call
                    # will have already generated an exception to that effect.
                    _LOGGER.warning(
                            "Failed attempt to parse chromosomes as read "
                            "chunk keys; arbitrarily chunking reads instead.")
                    read_chunk_keys = self.chunk_reads(
                            readsfile, chunksize=chunksize)
                else:
                    if interleave_chunk_sizes:
                        # Interleave chromosomes by size so that if tasks are
                        # pre-allocated to workers, we'll get about even bins.
                        read_chunk_keys = interleave_chromosomes_by_size(
                                size_by_chromosome.items())
                    else:
                        read_chunk_keys = size_by_chromosome.keys()

        _LOGGER.info("Temporary files will be stored in: '{}'".
                     format(self.temp_folder))
        _LOGGER.info("Processing with {} cores...".format(self.cores))

        # Some implementors may have a strand mode attribute.
        # If so, log it here to avoid duplicate messaging, as it
        # will remain constant across processed chunks (chromosomes).
        try:
            _LOGGER.info("STRAND MODE: {}".format(self.use_strand))
        except AttributeError:
            pass

        # TODO: handle non-chromosome-based case.
        idxstats = readsfile.get_index_statistics()
        reads_by_chrom = {istat.contig: istat.total for istat in idxstats}
        empties, nonempties = [], []
        for c in read_chunk_keys:
            target = empties if 0 == reads_by_chrom[c] else nonempties
            target.append(c)

        # Maps for order preservation. This permits arbitrary result return,
        # i.e. something other than the chunk key itself, when the process
        # completes. It would be a bit simpler to filter on the results
        # directly, but clients may want to implement a processor that
        # has a result with meaning beyond a signal/flag that it succeeded
        # for a particular chunk ID. That is, it may produce a result with
        # downstream meaning, and not be used simply for effect on disk.
        if self.cores == 1:
            results = map(self, nonempties)
        else:
            workers = multiprocessing.Pool(self.cores)
            # The typical call to map fails to acknowledge KeyboardInterrupts.
            # This fix helps: http://stackoverflow.com/a/1408476/946721

            
            _LOGGER.debug("Cores: '{}'".format(self.cores))
            _LOGGER.debug("Nonempties: '{}'".format(nonempties))
            results = workers.map_async(self, nonempties).get(9999999)

        # TODO: note the dependence on order here.
        result_by_chunk = [(c, self.empty_action(c)) for c in empties] + \
                           list(zip(nonempties, results))
        bad_chunks, good_chunks = \
                partition_chunks_by_null_result(result_by_chunk)

        if bad_chunks:
            _LOGGER.debug("Discarding {} chunk(s) of reads: {}".
                         format(len(bad_chunks), bad_chunks))
            _LOGGER.debug("Keeping {} chunk(s) of reads: {}".
                         format(len(good_chunks), good_chunks))
        else:
            _LOGGER.info("Using all reads")

        return good_chunks

    def fetch_chunk(self, chromosome):
        """
        Pull a chunk of sequencing reads from a file.
        
        :param str chromosome: identifier for chunk of reads to select.
        :return Iterable[pysam.AlignedSegment]: collection of aligned reads
        """
        if not self.by_chromosome:
            raise NotImplementedError(
                    "Provide a fetch_chunk implementation "
                    "if not partitioning reads by chromosome.")
        readsfile = PARA_READ_FILES[READS_FILE_KEY]
        return readsfile.fetch(chromosome, multiple_iterators=True)

    def combine(self, good_chromosomes, strict=False, chrom_sep=None):
        """
        Aggregate output from independent read chunks into single output file.
        
        :param Iterable[str] good_chromosomes: identifier (e.g., chromosome)
            for each chunk of reads processed.
        :param bool strict: whether to throw an exception upon encountering a
            missing file. If not, simply log a warning message and continue the
            aggregation process that's underway, working with what is available.
        :param str chrom_sep: delimiter between output from each chromosome.
        :return Iterable[str]: path to each file successfully combined.
        :raise pararead.exceptions.MissingOutputFileException: if executing in
            strict mode, and there's a reads chunk key for which the derived
            filepath does not exist.
        :raise pararead.exceptions.IllegalChunkException: if a chunk of reads
            outside of those declared to be of interest is requested to
            participate in the combination.
        """

        if not good_chromosomes:
            _LOGGER.warning("No successful chromosomes, so no combining.")
            return

        # Check that the combination request accords with the chunks
        # declared to be of interest
        if self.limit:
            missing_chunks = set(good_chromosomes) - set(self.limit)
            if missing_chunks:
                raise IllegalChunkException(
                        requested=missing_chunks, of_interest=self.limit)

        _LOGGER.info("Merging {} files into output file: '{}'".
                     format(len(good_chromosomes), self.outfile))

        # Track what we actually combine (particularly if non-strict
        # with respect to chunk(s) for which output file is missing.
        paths_combined_files = []

        if len(good_chromosomes) == 1:
            _LOGGER.debug("Just one good chromosome; ignoring delimiter.")
            chrom_sep = None

        with open(self.outfile, 'w') as outfile:
            for chrom in good_chromosomes:
                reads_chunk_output = self._tempf(chrom)

                # Handle case in which chunk's output is missing.
                if not os.path.exists(reads_chunk_output):
                    if strict:
                        raise MissingOutputFileException(
                                reads_chunk_key=chrom,
                                filepath=reads_chunk_output)
                    else:
                        _LOGGER.warning(
                                "Missing output file for reads chunk '%s', "
                                "skipping: '%s'", chrom, reads_chunk_output)
                        continue

                # Append lines from this chunk's output.
                with open(reads_chunk_output, 'r') as tmpf:
                    for line in tmpf:
                        outfile.write(line)
                if chrom_sep:
                    outfile.write(chrom_sep)
                paths_combined_files.append(reads_chunk_output)

        return paths_combined_files

    @pending_feature
    def chunk_reads(self, readsfile, chunksize=None):
        """
        Partition sequencing reads into equally-sized 'chunks' for 
        parallel processing. This treats all reads equally, in that 
        they are grouped by contiguous occurrence within the given 
        input file. This facilitates processing of unaligned reads, 
        but it means that reads from the same chromosome will not be 
        processed together. This can be overridden if that's desired.

        :param Iterable readsfile: collection of reads to chunk, likely
            pysam.AlignmentFile pysam.VariantFile; reads to split into chunks.
        :param int chunksize: number of units (i.e., reads) per processing chunk.
            If unspecified, this is derived using the instance's 
            cores count and chunks-per-core parameter, along with 
            a count of the number of units (reads). Note that if 
            this is unspecified, there will be some additional time 
            used to count the reads to derive chunk size.
        :return Iterable[(int, itertools.groupby)]: pairs of chunk key/ID and
            chunk reads chunk itself.
        """
        if not chunksize or chunksize < 1:
            _, reads_clone = itertools.tee(readsfile)
            num_chunks = self.cores * CHUNKS_PER_CORE

            # Count the reads.
            _LOGGER.info("Deriving chunk size for %d chunks: "
                         "%d cores x %d chunks/core",
                         num_chunks, self.cores, CHUNKS_PER_CORE)
            _, reads_clone = itertools.tee(readsfile)
            num_reads = sum(1 for _ in reads_clone)
            _LOGGER.info("Reads count: %d", num_reads)

            chunksize = int(num_reads / num_chunks)

        return itertools.groupby(
            enumerate(readsfile), key=lambda ipair: int(ipair[0] / chunksize))

    def _tempf(self, chrom):
        """
        Derive name for temporary file from chromosome name.

        :param str chrom: name of chromosome (single processing partition).
        :return str: name for tempfile corresponding to given unit name.

        """
        return os.path.join(
                self.temp_folder,
                "{}.{}".format(chrom or "ALL", self.intermediate_output_type))
