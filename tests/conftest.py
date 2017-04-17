""" Configuration of tests """

import pytest
import pysam

from pararead import processor, setup_pararead_logger
from pararead.logs import DEV_LOGGING_FMT
from tests import \
    IS_ALIGNED_PARAM_NAME, NAME_TEST_LOGFILE, \
    PATH_ALIGNED_FILE, PATH_UNALIGNED_FILE
from tests.helpers import IdentityProcessor, ReadsfileWrapper


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



def pytest_generate_tests(metafunc):
    """ Additional runtime parameterization of test cases. """
    if "num_cores" in metafunc.fixturenames:
        metafunc.parametrize(argnames="num_cores", argvalues=[1, 2, 4])



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



@pytest.fixture(scope="function")
def path_logs_file(request, tmpdir):

    logfile = tmpdir.join(NAME_TEST_LOGFILE).strpath
    logger = setup_pararead_logger(
            logfile=logfile, stream_format=DEV_LOGGING_FMT)

    def clear_handlers():
        logger.handlers = []

    request.addfinalizer(clear_handlers)
    return logfile
