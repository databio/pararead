""" Basic tests for ParaRead """

import itertools
import pytest
from pararead.processor import ParaReadProcessor
from tests import PATH_ALIGNED_FILE, PATH_UNALIGNED_FILE
from tests.helpers import IdentityProcessor


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



class ConstructorTests:
    """
    Basic tests for ParaReadProcessor.
    """

    def test_is_abstract(self):
        with pytest.raises(TypeError) as exc:
            # Provide filler arguments for ParaReadProcessor parameters
            # in an effort to ensure that the TypeError comes from the
            # requirement that the class is abstract, rather than from
            # missing arguments for required parameters.
            ParaReadProcessor(path_reads_file="dummy.bam",
                              cores=4, outfile="dummy.txt")
        # As a fallback, check that the exception message mentions "abstract."
        assert "abstract" in exc.value.message


    @pytest.mark.parametrize(
            argnames=["filepath", "num_cores"],
            argvalues=itertools.product(
                    [PATH_ALIGNED_FILE, PATH_UNALIGNED_FILE], [1, 4]),
            ids=lambda (fp, cores): "{}; cores={}".format(fp, cores))
    def test_requires_outfile_or_action(self, filepath, num_cores):
        """ Explicit output file or action name to derive one is needed. """
        with pytest.raises(ValueError):
            IdentityProcessor(filepath, cores=num_cores)



class FilesystemTests:
    """ Tests regarding interaction between Processor and filesystem """


    def test_removes_tempfolder(self, identity_processor):
        """ Folder for temporary files should be removed. """
        pass

    def test_closes_readsfile(self):
        pass




class FileRegistrationTests:
    """ Tests for registration of files with the ParaReadProcessor. """


    def test_filetype_inference(self):
        pass


    def test_filetype_exception(self):
        pass



class ExecutionTests:
    """ Tests for processor's run() method. """


    def test_KeyError_if_unregistered(self):
        pass


    def test_strand_mode_messaging(self):
        """ Child class with strand mode provides a message. """
        pass


    def test_cores_count(self):
        pass


    def test_chunksize_inference(self):
        pass


    def test_fixed_chunksize(self):
        pass



class CombinerTests:
    """ Processor provides function to combine intermediate results. """
    pass
