""" Test the utility functions. """

from functools import partial
import itertools
import os

import pytest
import pysam

from pararead.exceptions import MissingHeaderException
from pararead.utils import \
    chromosomes_from_bam_header, partition_chunks_by_null_result
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
    def test_filters_chromosomes(self, chromosomes, aligned_reads_file):
        """ Chromosome name-from-BAM-header fetch can filter. """
        observed = chromosomes_from_bam_header(
                readsfile=aligned_reads_file, chroms=chromosomes)
        assert list(chromosomes) == list(observed)


    @pytest.mark.parametrize(
            argnames=["require_aligned", "is_aligned"],
            argvalues=itertools.product([False, True], [False, True]))
    def test_allows_alignment_requirement_flexibility(
            self, require_aligned, is_aligned, path_test_data):
        """ The chromosome fetcher affords option to require aligned input. """

        # Create the reads file.
        name_reads_file = self.FILE_BY_ALIGNMENT[is_aligned]
        path_reads_file = os.path.join(path_test_data, name_reads_file)
        readsfile = pysam.AlignmentFile(
            path_reads_file, mode='rb', check_sq=False)

        # Parameterize the function call
        func = partial(chromosomes_from_bam_header,
                       readsfile=readsfile, require_aligned=require_aligned)

        # Determine expected behavior.
        if not is_aligned and require_aligned:
            # Only 1 of the 4 cases is exceptional.
            with pytest.raises(MissingHeaderException):
                func()
        else:
            # Unaligned input returns null chromosome names collection;
            # aligned input retains all no chromosomes with no filtration.
            expected = self.CHROMOSOME_NAMES if is_aligned else None
            assert expected == func()


    def test_retains_no_chromosomes(self, aligned_reads_file):
        """ Set of chromosomes disjoint from those present retains nothing. """
        retain_chroms = ["not-present-1", "chrNull"]
        assert [] == chromosomes_from_bam_header(
                readsfile=aligned_reads_file, chroms=retain_chroms)



class PartitionChromosomesByResultTests:
    """ Tests for filtering a collection of results """


    def test_discards_nulls_keeps_non_nulls(self):
        result_by_chromosome = {
                "chr1": None, "chr2": "non-null",
                "chrX": ["totally", "arbitrary"], "chrY": None}
        expected_scraps = ["chr1", "chrY"]
        expected_keeps = ["chr2", "chrX"]
        scraps, keeps = \
                partition_chunks_by_null_result(result_by_chromosome)
        # Sort is for comparison in case of Mapping rather than Sequence.
        assert expected_scraps == scraps
        assert expected_keeps == keeps
