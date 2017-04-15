""" Parallel reads processor utilities. """

import itertools
import operator as op
import sys
if sys.version_info < (3, 3):
    from collections import Mapping, Sequence
else:
    from collections.abc import Mapping, Sequence
from .exceptions import MissingHeaderException


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


__all__ = ["interleave_chromosomes_by_size",
           "make_outfile_name", "parse_bam_header",
           "partition_chunks_by_null_result",
           "pending_feature", "unbuffered_write"]



def interleave_chromosomes_by_size(size_by_chromosome):
    """
    Arrange chromosome names to facilitate even binning.
    
    Intersperse/interleave chromosomes such that ones at opposite ends of 
    the spectrum of sizes are adjacent.
    
    Parameters
    ----------
    size_by_chromosome : Iterable of (str, int) or Mapping[str, int]

    Returns
    -------
    Iterable of str
        Names of chromosomes

    """

    if not size_by_chromosome:
        return []
    if isinstance(size_by_chromosome, Mapping):
        size_by_chromosome = size_by_chromosome.items()

    ordered_chromosomes = zip(*sorted(size_by_chromosome,
                                      key=op.itemgetter(1)))[0]
    num_chromosomes = len(ordered_chromosomes)
    meridian = int(num_chromosomes / 2)
    first_half, second_half = \
            ordered_chromosomes[:meridian], ordered_chromosomes[meridian:]

    interleaved = list(itertools.chain(*zip(first_half, second_half[::-1])))

    num_interleaved = len(interleaved)
    if num_interleaved != num_chromosomes:
        # Account for odd number of chromosomes; zip will truncate.
        assert num_interleaved == num_chromosomes - 1
        maybe_missing_chrom = ordered_chromosomes[-1]
        assert maybe_missing_chrom not in interleaved
        interleaved.append(maybe_missing_chrom)

    return interleaved



def make_outfile_name(readsfile_basename, processing_action, output_type):
    """
    Create a name for an output file based on action performed.

    Parameters
    ----------
    readsfile_basename : str
        Path-less and extension-less version of name of file of reads.
    processing_action : str
        Name for the processing action being performed by a
        parallel processor of sequencing reads (i.e., a class
        derived from ParaReadProcessor).
    output_type : str
        Type of output file for which name is being created.

    Returns
    -------
    str
        (Fallback) name for output file, used by the
        ParaReadProcessor constructor if a null or
        empty output filename is provided at creation.

    """
    return "{}_{}.{}".format(readsfile_basename,
                             processing_action, output_type)



def parse_bam_header(readsfile, chroms=None, require_aligned=False):
    """
    Get a list of chromosomes (and lengths) in this readsfile from header.

    Parameters
    ----------
    readsfile : pysam.libcalignmentfile.AlignmentFile
        File with aligned sequencing reads datasets
    chroms : collections.Iterable of str, optional
        Chromosomes of interest. If empty/null/False,
        assume that all chromosomes with read(s) are of interest.
    require_aligned : bool, default False
        Whether to throw an exception if given unaligned input

    Returns
    -------
    None or Iterable of str, int
        Null if no chromosomes are in the header

    """

    try:
        all_sizes_by_chrom = {headline['SN']: headline['LN']
                              for headline in readsfile.header['SQ']}
    except KeyError:
        all_sizes_by_chrom = {}

    if not all_sizes_by_chrom:
        if require_aligned:
            raise MissingHeaderException(readsfile.filename)
        else:
            return None

    if not chroms:
        # No filtration --> retain all chromosomes.
        return all_sizes_by_chrom

    return {c: s for c, s in all_sizes_by_chrom.items() if c in set(chroms)}



def partition_chunks_by_null_result(result_by_chromosome):
    """
    Bin chromosome name by whether processing result was null.

    Parameters
    ----------
    result_by_chromosome : Sequence of (str, object) or Mapping[str, object]
        Mapping from chromosome name to processing result

    Returns
    -------
    (list of str, list of str)
        Sequence of names of chromosomes for which result was null,
        and an analogous sequence for those with a null result.

    """
    # Ideally, the filtration strategy would be an argument.
    # Due to cPickle's disdain for anonymous functions and general
    # inability to serialize a callable, though, the filtration
    # strategy employed is fixed to be whether the result is null.
    bad_chroms, good_chroms = [], []
    if isinstance(result_by_chromosome, Mapping):
        res_by_chr = result_by_chromosome.items()
        if not isinstance(result_by_chromosome, Sequence):
            # Non-mapping input suggests pre-sorting, so sort Mapping to
            # facilitate concordance between results for varied input types.
            # Don't destroy order of an already-ordered Mapping, though.
            res_by_chr = sorted(res_by_chr)
    else:
        res_by_chr = result_by_chromosome
    for c, r in res_by_chr:
        chroms = bad_chroms if r is None else good_chroms
        chroms.append(c)
    return bad_chroms, good_chroms



def pending_feature(not_yet_implemented):
    def raise_error(*args, **kwargs):
        raise NotImplementedError("{} is not fully implemented".
                                  format(not_yet_implemented.__name__))
    return raise_error



def unbuffered_write(txt):
    """ Writes unbuffered output by flushing after each stdout.write call """
    sys.stdout.write(txt)
    sys.stdout.flush()
