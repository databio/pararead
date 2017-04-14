""" Configuration of tests """

import pytest
import pysam

from tests import \
    CORES_PARAM_NAME, IS_ALIGNED_PARAM_NAME, \
    PATH_ALIGNED_FILE, PATH_UNALIGNED_FILE
from tests.helpers import IdentityProcessor, ReadsfileWrapper

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


"""
def pytest_generate_tests(metafunc):
    if "require_aligned" in metafunc.fixturenames:
        metafunc.parametrize(
            argnames=["is_aligned", "path_reads_file"],
            argvalues=[(False, PATH_UNALIGNED_FILE),
                       (True, PATH_ALIGNED_FILE)],
            ids=lambda (is_aln, rfile): "{} ({})".format(
                    rfile.filename, "aligned" if is_aln else "unaligned"))
"""


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
def identity_processor(request):

    if CORES_PARAM_NAME not in request.fixturenames:
        raise Exception("To create a processor, test case must "
                        "specify '{}' as a fixture".format(CORES_PARAM_NAME))
    num_cores = request.getfixturevalue(CORES_PARAM_NAME)

    if IS_ALIGNED_PARAM_NAME in request.fixturenames \
            and not request.getfixturevalue(IS_ALIGNED_PARAM_NAME):
        path_reads_file = PATH_UNALIGNED_FILE
    else:
        path_reads_file = PATH_ALIGNED_FILE

    return IdentityProcessor(path_reads_file, cores=num_cores)
