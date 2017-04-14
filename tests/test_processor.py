""" Basic tests for ParaRead """

import pytest
from pararead.processor import ParaReadProcessor


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



class ParaReadProcessorTests:
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
