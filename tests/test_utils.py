""" Test the utility functions. """

from functools import partial
import itertools

import pytest
from pysam import AlignmentFile, VariantFile

from pararead.exceptions import FileTypeException, MissingHeaderException
from pararead.utils import (
    create_reads_builder,
    parse_bam_header,
    partition_chunks_by_null_result,
    READS_FILE_MAKER,
)
from tests import PATH_ALIGNED_FILE, PATH_UNALIGNED_FILE
from tests.helpers import ReadsfileWrapper


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


class CreateReadsBuilderTests:
    """Tests for creation of reads file factory builder."""

    @pytest.mark.parametrize(
        argnames="unsupported_filetype",
        argvalues=list(
            itertools.chain(
                *[(ft.lower(), ft.upper()) for ft in ["BED", "TXT", "GTF", "TSV"]]
            )
        ),
    )
    def test_filetype_exception(self, unsupported_filetype):
        """Unsupported input filetype is exceptional."""
        filename = "testfile.{}".format(unsupported_filetype)
        with pytest.raises(FileTypeException):
            create_reads_builder(filename)

    @pytest.mark.parametrize(
        argnames="filetype",
        argvalues=list(
            itertools.chain(
                *[(ft.lower(), ft.upper()) for ft in READS_FILE_MAKER.keys()]
            )
        ),
    )
    def test_infers_filetype(self, filetype):
        """We infer pysam type from file extension."""
        ft_key = filetype.upper()
        filename = "testfile.{}".format(filetype)
        readsfile_builder = create_reads_builder(filename)
        inferred_constructor = readsfile_builder.ctor
        if ft_key in ["BAM", "SAM", "CRAM"]:
            assert AlignmentFile == inferred_constructor
        elif ft_key in ["BCF", "VCF"]:
            assert VariantFile == inferred_constructor
        else:
            pytest.fail("Add test case for filetype: '{}'").format(filetype)


class ParseBamHeaderTests:
    """Tests for fetching chromosome names from BAM file header."""

    CHROMOSOME_NAMES = ["K1_unmethylated", "K3_methylated"]

    @pytest.mark.parametrize(
        argnames="chromosomes", argvalues=[{"K1_unmethylated"}, ("K3_methylated",)]
    )
    def test_keeps_specific_chromosomes(self, chromosomes, aligned_reads_file):
        """Chromosome name-from-BAM-header fetch grabs requested chroms."""
        observed = parse_bam_header(
            readsfile=aligned_reads_file, chroms=chromosomes
        ).keys()
        assert set(chromosomes) == set(observed)

    @pytest.mark.parametrize(
        argnames=["require_aligned", "is_aligned"],
        argvalues=itertools.product([False, True], [False, True]),
    )
    def test_allows_alignment_requirement_flexibility(
        self, require_aligned, is_aligned
    ):
        """The chromosome fetcher affords option to require aligned input."""

        path_reads_file = PATH_ALIGNED_FILE if is_aligned else PATH_UNALIGNED_FILE
        with ReadsfileWrapper(path_reads_file, check_sq=is_aligned) as readsfile:

            # Parameterize the function call
            func = partial(
                parse_bam_header, readsfile=readsfile, require_aligned=require_aligned
            )

            # Determine expected behavior.
            if not is_aligned and require_aligned:
                # Only 1 of the 4 cases is exceptional.
                with pytest.raises(MissingHeaderException):
                    func()
            elif not is_aligned:
                assert func() is None
            else:
                # Unaligned input returns null chromosome names collection;
                # aligned input retains all no chromosomes with no filtration.
                assert set(self.CHROMOSOME_NAMES) == set(func())

    @pytest.mark.parametrize(argnames="is_aligned", argvalues=[False, True])
    def test_no_chromosome_intersection(self, is_aligned):
        """Set of chromosomes disjoint from those present retains nothing."""

        retain_chroms = ["not-present-1", "chrNull"]
        path_reads_file = PATH_ALIGNED_FILE if is_aligned else PATH_UNALIGNED_FILE

        with ReadsfileWrapper(path_reads_file, check_sq=is_aligned) as readsfile:
            observed = parse_bam_header(
                readsfile=readsfile, chroms=retain_chroms, require_aligned=False
            )

        expected = {} if is_aligned else None
        assert expected == observed


class PartitionChromosomesByResultTests:
    """Tests for filtering a collection of results"""

    def test_discards_nulls_keeps_non_nulls(self):
        result_by_chromosome = {
            "chr1": None,
            "chr2": "non-null",
            "chrX": ["totally", "arbitrary"],
            "chrY": None,
        }
        expected_scraps = ["chr1", "chrY"]
        expected_keeps = ["chr2", "chrX"]
        scraps, keeps = partition_chunks_by_null_result(result_by_chromosome)
        # Sort is for comparison in case of Mapping rather than Sequence.
        assert expected_scraps == scraps
        assert expected_keeps == keeps
