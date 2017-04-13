""" Parallel reads processor utilities. """


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


__all__ = ["chromosomes_from_bam_header",
           "filter_chromosomes", "make_outfile_name"]



def chromosomes_from_bam_header(chroms, readsfile):
    """
    Get a list of chromosomes (and lengths) in this readsfile from header.

    Parameters
    ----------
    chroms : collections.Iterable of str
        Chromosomes of interest. If empty/null/False,
        assume that all chromosomes with read(s) are of interest.
    readsfile : pysam.libcalignmentfile.AlignmentFile
        File with aligned sequencing reads datasets

    Returns
    -------
    list of str
        Mapping from reference sequence name to
        reference sequence length for each reference
        sequence that intersects with those of interest

    """
    chroms = set(chroms or [])
    chrlist = [item['SN'] for item in readsfile.header['SQ']
               if not chroms or item['SN'] in chroms]
    return chrlist



def filter_chromosomes(result_by_chromosome):
    """
    Bin chromosome name by whether processing result was null.

    Parameters
    ----------
    result_by_chromosome : Sequence of (str, object)
        Mapping from chromosome name to processing result

    Returns
    -------
    (list of str, list of str)
        Sequence of names of chromosomes for which result was null,
        and an analogous sequence for those with a null result.

    """
    bad_chroms, good_chroms = [], []
    for c, r in result_by_chromosome:
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
