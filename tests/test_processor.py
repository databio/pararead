""" Basic tests for ParaRead """

import pytest
from pysam import AlignmentFile
from pararead.exceptions import ExecutionOrderException
from pararead.processor import ParaReadProcessor
from tests import \
        NUM_CORES_DEFAULT, NUM_READS_BY_FILE, \
        PATH_ALIGNED_FILE, PATH_UNALIGNED_FILE
from tests.helpers import IdentityProcessor


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



class ConstructorTests:
    """
    Basic tests for ParaReadProcessor.
    """

    def test_is_abstract(self):
        """ ParaReadProcessor must be extended. """
        with pytest.raises(TypeError) as exc:
            # Provide filler arguments for ParaReadProcessor parameters
            # in an effort to ensure that the TypeError comes from the
            # requirement that the class is abstract, rather than from
            # missing arguments for required parameters.
            ParaReadProcessor(path_reads_file="dummy.bam",
                              cores=NUM_CORES_DEFAULT, outfile="dummy.txt")
        # As a fallback, check that the exception message mentions "abstract."
        assert "abstract" in exc.value.message


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

        if require_aligned and not pysam_kwargs:
            with pytest.raises(ValueError):
                processor.register_files()
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
        with pytest.raises(ExecutionOrderException):
            processor.readsfile()

        # Now do the registration, creating the pysam readsfile instance.
        processor.register_files()
        readsfile = processor.readsfile()

        # Check out the new readsfile.
        assert isinstance(readsfile, AlignmentFile)
        num_reads = sum(1 for _ in readsfile)
        assert NUM_READS_BY_FILE[path_reads_file] == num_reads



class FilesystemTests:
    """ Tests regarding interaction between Processor and filesystem """


    @pytest.mark.skip("Implement for context manager use only.")
    def test_removes_tempfolder(self):
        """ Folder for temporary files should be removed. """
        pass


    @pytest.mark.skip("implement for context manager use only.")
    def test_closes_readsfile(self):
        pass



class ExecutionTests:
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



class CombinerTests:
    """ Processor provides function to combine intermediate results. """
    pass
