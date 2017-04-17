""" Tests' shared constants """

import os

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


DATA_PATH = os.path.join(os.path.dirname(__file__), "data")
PATH_ALIGNED_FILE = os.path.join(DATA_PATH, "meth_spikein_k1_k3.bam")
PATH_UNALIGNED_FILE = os.path.join(DATA_PATH, "rrbs1.bam")
IS_ALIGNED_PARAM_NAME = "is_aligned"
NUM_CORES_DEFAULT = 4
NAME_TEST_LOGFILE = "test.log"

NUM_READS_BY_FILE = {
    PATH_ALIGNED_FILE: 123,
    PATH_UNALIGNED_FILE: 1895
}
