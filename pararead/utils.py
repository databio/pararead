""" Parallel reads processor utilities. """

import itertools
import sys
if sys.version_info < (3, 3):
    from collections import Mapping, Sequence
else:
    from collections.abc import Mapping, Sequence
from .exceptions import MissingHeaderException

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


__all__ = ["chromosomes_from_bam_header", "make_outfile_name",
           "partition_chunks_by_null_result", "unbuffered_write"]



def chromosomes_from_bam_header(readsfile, chroms=None, require_aligned=False):
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
    None or list of str
        Names of chromosomes present in the input file that are of interest 
        as specified by the 

    """
    chroms_in_header = [headline['SN'] for headline in readsfile.header['SQ']]
    if not chroms_in_header:
        if require_aligned:
            raise MissingHeaderException(readsfile.filename)
        else:
            return None
    if not chroms:
        return chroms_in_header
    return [c for c in chroms_in_header if c in set(chroms)]



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



def make_outfile_name(readsfile_basename, processing_action):
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

    Returns
    -------
    str
        (Fallback) name for output file, used by the
        ParaReadProcessor constructor if a null or
        empty output filename is provided at creation.

    """
    return "{}_{}.txt".format(readsfile_basename, processing_action)



def unbuffered_write(txt):
    """ Writes unbuffered output by flushing after each stdout.write call """
    sys.stdout.write(txt)
    sys.stdout.flush()
