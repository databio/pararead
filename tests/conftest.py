""" Configuration of tests """

import pytest
import pysam

from pararead import processor
from tests import \
    IS_ALIGNED_PARAM_NAME, PATH_ALIGNED_FILE, PATH_UNALIGNED_FILE
from tests.helpers import IdentityProcessor, ReadsfileWrapper


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



@pytest.fixture(scope="function")
def aligned_reads_file():
    """
    Create aligned reads file.
    
    Returns
    -------
    pysam.AlignmentFile
        File with pysam.AlignedSegments for reads.

    """
    return pysam.AlignmentFile(PATH_ALIGNED_FILE)



@pytest.fixture(scope="function")
def unaligned_reads_file():
    """
    Create unaligned reads file

    Returns
    -------
    pysam.AlignmentFile
        File with pysam.AlignedSegments for reads.

    """
    return ReadsfileWrapper(PATH_UNALIGNED_FILE)



@pytest.fixture(scope="function", params=[1, 2, 4])
def num_cores(request):
    """ For quick tests, validate across several cores count values. """
    return request.param



@pytest.fixture(scope="function")
def identity_processor(request, num_cores, tmpdir):
    """
    Provide a basic processor for a fast test of some behavior.
    
    Parameters
    ----------
    request : pytest.fixtures.SubRequest
        Test case requesting the fixture parameterization.

    Returns
    -------
    pararead.ParaReadProcessor
        A very basic processor, returning elements with no 
        or very trivial modification(s) for speed.

    """

    if IS_ALIGNED_PARAM_NAME in request.fixturenames \
            and not request.getfixturevalue(IS_ALIGNED_PARAM_NAME):
        path_reads_file = PATH_UNALIGNED_FILE
    else:
        path_reads_file = PATH_ALIGNED_FILE

    path_output_file = tmpdir.join("placeholder-testfile.txt").strpath

    return IdentityProcessor(
            path_reads_file, cores=num_cores, outfile=path_output_file)



@pytest.fixture(scope="function")
def remove_reads_file(request):
    """ Clear out the reads file mapping after executing a test case. """
    def clear_pararead():
        processor.PARA_READ_FILES = {}
    request.addfinalizer(clear_pararead)
