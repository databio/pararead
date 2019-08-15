""" Basic tests for ParaRead """

import itertools
import os

import pytest
from pysam import AlignmentFile

from pararead.exceptions import \
    CommandOrderException, IllegalChunkException, \
    MissingHeaderException, MissingOutputFileException
from pararead.processor import ParaReadProcessor
from tests import \
    NUM_CORES_DEFAULT, NUM_READS_BY_FILE, \
    PATH_ALIGNED_FILE, PATH_UNALIGNED_FILE
from tests.helpers import IdentityProcessor, loglines


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


class ConstructorTests:
    """ Basic tests for ParaReadProcessor. """

    @pytest.mark.parametrize(
            argnames="filepath",
            argvalues=[PATH_ALIGNED_FILE, PATH_UNALIGNED_FILE])
    def test_requires_outfile_or_action(self, filepath):
        """ Explicit output file or action name to derive one is needed. """
        with pytest.raises(ValueError):
            IdentityProcessor(filepath)


class FileRegistrationTests:
    """ Tests for registration of files with the ParaReadProcessor. """

    @pytest.mark.parametrize(
            argnames="require_aligned", argvalues=[False, True])
    @pytest.mark.parametrize(
            argnames="pysam_kwargs", argvalues=[{}, {"check_sq": False}])
    def test_adds_pysam_kwargs(self, require_aligned,
                               pysam_kwargs, remove_reads_file):
        """ Unaligned input BAM needs check_sq=False to be created. """

        # Note that remove_reads_file is here to clear the module-scoped map.

        # Explicitly set by_chromosome=False to prevent it from
        # controlling the requirement regarding aligned reads.
        processor = IdentityProcessor(
                path_reads_file=PATH_UNALIGNED_FILE, action="test",
                allow_unaligned=not require_aligned, by_chromosome=False)

        if require_aligned:
            exp_error = MissingHeaderException if pysam_kwargs else ValueError
            with pytest.raises(exp_error):
                processor.register_files(**pysam_kwargs)
        else:
            # No exception --> pass (file registration is just for effect.)
            processor.register_files(**pysam_kwargs)

    @pytest.mark.parametrize(
        argnames=["path_reads_file", "require_aligned"],
        argvalues=[(PATH_ALIGNED_FILE, False), (PATH_ALIGNED_FILE, True),
                   (PATH_UNALIGNED_FILE, False)])
    def test_creates_fresh_reads_file(self, path_reads_file,
                                      require_aligned, remove_reads_file):
        """ Reads file pysam object is created by register_files(). """

        # Note that remove_reads_file is here to clear the module-scoped map.

        # Explicitly set by_chromosome=False to prevent it from
        # controlling the requirement regarding aligned reads.
        processor = IdentityProcessor(
                path_reads_file=path_reads_file, action="test",
                allow_unaligned=not require_aligned, by_chromosome=False)

        # The pysam readsfile shouldn't exist before register_files().
        with pytest.raises(CommandOrderException):
            processor.readsfile

        # Now do the registration, creating the pysam readsfile instance.
        processor.register_files()
        readsfile = processor.readsfile

        # Check out the new readsfile.
        assert isinstance(readsfile, AlignmentFile)
        num_reads = sum(1 for _ in readsfile)
        assert NUM_READS_BY_FILE[path_reads_file] == num_reads


class CombinerTests:
    """ Processor provides function to combine intermediate results. """

    # NOTE: seemingly-unused fixtures are present for effectful-ness.

    CHROMOSOME_CHUNK_KEY = "chromosome"
    ARBITRARY_CHUNK_KEY = "arbitrary"
    CHROM_NAMES = ["chr{}".format(i) for i in range(1, 23)] + \
                  ["chrX", "chrY", "chrM"]
    ARBITRARY_NAMES = ["random0", "arbitrary1", "contig2"]
    COMBO_REQUEST_NAMES = CHROM_NAMES + ARBITRARY_NAMES
    CHUNK_NAMES = {CHROMOSOME_CHUNK_KEY: CHROM_NAMES,
                   ARBITRARY_CHUNK_KEY: ARBITRARY_NAMES}

    @pytest.mark.parametrize(
            argnames="error_if_missing", argvalues=[False, True])
    @pytest.mark.skip
    def test_nothing_to_combine(self, tmpdir, path_logs_file,
                                num_cores, error_if_missing):
        """ Complete lack of output is sufficient to warrant a warning. """

        # Create the processor and do combine() step.
        path_output_file = tmpdir.join("output.txt").strpath
        processor = IdentityProcessor(
                PATH_ALIGNED_FILE, cores=num_cores, outfile=path_output_file)
        num_logs_before_combine = len(loglines(path_logs_file))
        processor.combine(good_chromosomes=[], strict=error_if_missing)
        log_records = loglines(path_logs_file)

        # The log record should be a warning, and there's only one.
        assert 1 == len(log_records) - num_logs_before_combine
        assert "WARN" in log_records[num_logs_before_combine]

    @pytest.mark.parametrize(
            argnames="which_names",
            argvalues=[CHROMOSOME_CHUNK_KEY, ARBITRARY_CHUNK_KEY])
    def test_missing_output_files(
            self, which_names, extant_files, fixed_tempfolder_processor):
        """ Missing-output chunks be skipped or exceptional. """
        with pytest.raises(MissingOutputFileException):
            fixed_tempfolder_processor.combine(self.COMBO_REQUEST_NAMES, strict=True)

    @pytest.mark.parametrize(
        argnames="which_names",
        argvalues=[CHROMOSOME_CHUNK_KEY, ARBITRARY_CHUNK_KEY])
    def test_missing_output_files_non_strict_retval(
            self, which_names, extant_files,
            fixed_tempfolder_processor, path_logs_file):
        """ Combiner returns just the paths that were used. """
        observed_combined_filepaths = fixed_tempfolder_processor.combine(
                self.COMBO_REQUEST_NAMES, strict=False)
        assert extant_files == observed_combined_filepaths

    @pytest.mark.parametrize(
        argnames="which_names",
        argvalues=[CHROMOSOME_CHUNK_KEY, ARBITRARY_CHUNK_KEY])
    @pytest.mark.skip
    def test_missing_output_files_non_strict_messaging(
            self, which_names, extant_files,
            fixed_tempfolder_processor, path_logs_file):
        """ If non-strict, combiner warns about requested-but-missing. """

        # Do the combine step and get the logged messages.
        num_logs_before_combine = len(loglines(path_logs_file))
        fixed_tempfolder_processor.combine(
                self.COMBO_REQUEST_NAMES, strict=False)
        logs_from_combine = loglines(path_logs_file)[num_logs_before_combine:]

        # As a control, check that we are in fact over-requesting in combine().
        num_extant_files = len(extant_files)
        num_requested_files = len(self.COMBO_REQUEST_NAMES)
        assert num_extant_files < num_requested_files

        # The control makes this assertion meaningful.
        num_skips_expected = num_requested_files - num_extant_files
        num_warns_observed = \
                sum(1 for msg in logs_from_combine if "WARN" in msg)
        assert num_skips_expected == num_warns_observed

    @pytest.mark.parametrize(
        argnames=["filetype", "combined_output_type"],
        argvalues=list(itertools.product(["bed", "tsv"], ["bed", "tsv"])))
    @pytest.mark.parametrize(
        argnames="which_names",
        argvalues=[CHROMOSOME_CHUNK_KEY, ARBITRARY_CHUNK_KEY])
    def test_different_format(self, tmpdir, filetype, combined_output_type,
                              which_names, extant_files, num_cores):
        """ File content is actually combined, and formats can differ. """

        # Manual creation of the processor here to control output type.
        path_output_file = tmpdir.join(
                "testfile.{}".format(combined_output_type)).strpath
        processor = IdentityProcessor(
                PATH_ALIGNED_FILE, cores=num_cores,
                outfile=path_output_file, intermediate_output_type=filetype)
        processor.temp_folder = tmpdir.strpath

        # Write to the dummy output file for each chunk.
        expected_lines = {fp: "file{}: {}\n".format(i, fp)
                          for i, fp in enumerate(extant_files)}
        for fp, line in expected_lines.items():
            with open(fp, 'w') as f:
                f.write(line)

        # For control, enforce that combined output doesn't already exist.
        assert not os.path.exists(path_output_file)
        processor.combine(self.CHUNK_NAMES[which_names], strict=True)
        assert os.path.isfile(path_output_file)

        # Check that output was combined accurately.
        with open(path_output_file, 'r') as combined:
            observed_lines = combined.readlines()
        assert set(expected_lines.values()) == set(observed_lines)

    @pytest.mark.parametrize(
        argnames="which_names",
        argvalues=[CHROMOSOME_CHUNK_KEY, ARBITRARY_CHUNK_KEY])
    def test_enforces_chunks_limit(
            self, which_names, extant_files, 
            fixed_tempfolder_processor, path_logs_file):
        """ Combination applies only to chunks of interest. """

        # Tell the processor that only certain chunks are of interest.
        extant_read_chunks = self.CHUNK_NAMES[which_names]
        fixed_tempfolder_processor.limit = extant_read_chunks

        # Request an uninteresting chunk in the combine() step.
        # This is erroneous because it would not have been processed.
        bad_chunk_name = "not-a-chunk"
        with pytest.raises(IllegalChunkException) as error:
            fixed_tempfolder_processor.combine(
                    extant_read_chunks + [bad_chunk_name])
        assert bad_chunk_name in str(error.value)

    @pytest.fixture(scope="function")
    def extant_files(self, request, tmpdir):
        """
        Ensure the existence of certain files for a test case.

        The processor's combine step aggregates results from individually 
        and independently processed reads chunks. Thus, it needs a file to 
        exist for each chunk that it attempts to include in the final 
        aggregated output. This creates empty files named according to 
        the components that a test case specifies it's interested in combining.

        Parameters
        ----------
        request : pytest.fixtures.SubRequest
            Test case requesting the parameterization.
        tmpdir : py._path.local.LocalPath
            Path to temporary folder for the test case.

        """

        # The test case may be parameterized with respect to chunk names.
        if "which_names" in request.fixturenames:
            chunk_names_key = request.getfixturevalue("which_names")
            chunk_names = self._names_from_key(chunk_names_key)
        else:
            chunk_names = self.CHROM_NAMES

        # The test case may also be parameterized with respect to file type.
        if "filetype" in request.fixturenames:
            extension = request.getfixturevalue("filetype")
        else:
            extension = "txt"

        # "Touch" each file, storing the corresponding path to 
        # communicate back to the requesting test case.
        files = []
        for chunk in chunk_names:
            path_out_file = tmpdir.join("{}.{}".format(chunk, extension))
            path_out_file.ensure(file=True)
            files.append(path_out_file.strpath)
        return files

    @pytest.fixture(scope="function")
    def fixed_tempfolder_processor(self, tmpdir, num_cores):
        """
        Create processor with known temporary folder.

        Provide some basic required values to the constructor, and critically, 
        fix the temporary folder to the location into which the individual 
        reads chunk files were ensured to exist. The processor constructor 
        provides tempfolder parent as a parameter, but then it creates a 
        temporary folder within that, used to search for the individual 
        chunk output files. Those files have already been created in a test 
        temp folder, though, so the processor needs to know about that. 
        Cleanup of both locations is handled. The processor registers its 
        own temporary folder for removal while pytest cleans up its folder.

        Parameters
        ----------
        tmpdir : py._path.local.LocalPath
            Path to where chunks' dummy output files have been placed.
        num_cores : int
            Number of cores to use for the test case, parameterized.

        Returns
        -------
        pararead.ParaReadProcessor
            New processor instance, with updated temp folder knowledge.

        """
        path_output_file = tmpdir.join("test-output.txt").strpath
        processor = IdentityProcessor(
            PATH_ALIGNED_FILE, cores=num_cores, outfile=path_output_file)
        processor.temp_folder = tmpdir.strpath
        return processor

    def _names_from_key(self, names_key):
        """
        Get chunk names based on an argument to a test case parameter.
        
        Parameters
        ----------
        names_key : str
            Argument to test case parameter indicating which set of chunk 
            identifiers will have files created.

        Returns
        -------
        Iterable of str
            Collection of chunk identifiers for which files will exist.

        """
        try:
            return self.CHUNK_NAMES[names_key]
        except KeyError:
            return self.CHROM_NAMES


@pytest.mark.skip("Not implemented")
class IntegrationTests:
    """ A couple of sample end-to-end tests through a simple processor. """
    pass


class FilesystemTests:
    """ Tests regarding interaction between Processor and filesystem """

    @pytest.mark.skip("Implement for context manager use only.")
    def test_removes_tempfolder(self):
        """ Folder for temporary files should be removed. """
        pass

    @pytest.mark.skip("implement for context manager use only.")
    def test_closes_readsfile(self):
        pass


class ArbitraryPartitionTests:
    """ Tests for processor's run() method. """

    @pytest.mark.skip("Not implemented")
    def test_cores_count(self):
        pass

    @pytest.mark.skip("Not implemented")
    def test_chunksize_inference(self):
        pass

    @pytest.mark.skip("Not implemented")
    def test_fixed_chunksize(self):
        pass
