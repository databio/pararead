""" Basic tests for ParaRead """

import pytest
from pararead.processor import ParaReadProcessor


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
            ParaReadProcessor(path_reads_file="dummy.bam", chunksize=1000,
                              cores=4, outfile="dummy.txt")
        # As a fallback, check that the exception message mentions "abstract."
        assert "abstract" in exc.value.message

    def test_requires_outfile_or_action(self):
        pass


    def test_removes_tempfolder(self):
        pass



class FileRegistrationTests:
    """ Tests for registration of files with the ParaReadProcessor. """


    def test_closes_readsfile(self):
        pass


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
    