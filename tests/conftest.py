""" Configuration of tests """

import os
import pytest
import pysam

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


NAME_ALIGNED_FILE = "meth_spikein_k1_k3.bam"
NAME_UNALIGNED_FILE = "rrbs1.bam"


@pytest.fixture(scope="session")
def path_test_data():
    """
    Provide shared path to data for tests.
    
    Returns
    -------
    str
        Path to data for tests.

    """
    return os.path.join(os.path.dirname(__file__), "data")


@pytest.fixture(scope="function")
def aligned_reads_file(path_test_data):
    """
    Create path to aligned reads file.
    
    Parameters
    ----------
    path_test_data : str
        Path to tests' data folder.

    Returns
    -------
    pysam.AlignmentFile
        File with pysam.AlignedSegments for reads.

    """
    path_reads_file = os.path.join(path_test_data, NAME_ALIGNED_FILE)
    return pysam.AlignmentFile(path_reads_file, mode='rb', check_sq=False)



@pytest.fixture(scope="function")
def unaligned_reads_file(path_test_data):
    """
    
    Parameters
    ----------
    path_test_data : str
        Path to folder containing test data.

    Returns
    -------
    pysam.AlignmentFile
        File with pysam.AlignedSegments for reads.

    """
    path_reads_file = os.path.join(path_test_data, NAME_UNALIGNED_FILE)
    return pysam.AlignmentFile(path_reads_file, mode='rb', check_sq=False)
