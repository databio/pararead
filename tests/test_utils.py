""" Test the utility functions. """

from functools import partial
import itertools
import os

import pytest
import pysam

from pararead.exceptions import MissingHeaderException
from pararead.utils import \
    chromosomes_from_bam_header, partition_chromosomes_by_null_result
from conftest import NAME_ALIGNED_FILE, NAME_UNALIGNED_FILE

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



class ChromosomesFromBamHeaderTests:
    """ Tests for fetching chromosome names from BAM file header. """

    CHROMOSOME_NAMES = ["K1_unmethylated", "K3_methylated"]
    FILE_BY_ALIGNMENT = {True: NAME_ALIGNED_FILE, False: NAME_UNALIGNED_FILE}
    CHROMOSOMES_BY_FILE = {NAME_ALIGNED_FILE: CHROMOSOME_NAMES}


    @pytest.mark.parametrize(
            argnames="chromosomes",
            argvalues=[{"K1_unmethylated"}, ("K3_methylated", )])
    def test_specificity(self, chromosomes, aligned_reads_file):
        observed = chromosomes_from_bam_header(
                readsfile=aligned_reads_file, chroms=chromosomes)
        assert list(chromosomes) == list(observed)


    @pytest.mark.parametrize(
            argnames=["require_aligned", "is_aligned"],
            argvalues=itertools.product([False, True], [False, True]))
    def test_alignment_requirement(self, require_aligned,
                                   is_aligned, path_test_data):

        name_reads_file = self.FILE_BY_ALIGNMENT[is_aligned]
        path_reads_file = os.path.join(path_test_data, name_reads_file)
        readsfile = pysam.AlignmentFile(
            path_reads_file, mode='rb', check_sq=False)

        func = partial(chromosomes_from_bam_header,
                       readsfile=readsfile, require_aligned=require_aligned)
        if not is_aligned and require_aligned:
            with pytest.raises(MissingHeaderException):
                func()
        else:
            expected = self.CHROMOSOME_NAMES if is_aligned else None
            assert expected == func()


    def retains_no_chromosomes(self):
        pass


class PartitionChromosomesByResultTests:
    """ Tests for filtering a collection of results """
    pass
