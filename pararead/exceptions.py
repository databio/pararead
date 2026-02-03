"""Specific exception types."""

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


class CommandOrderException(Exception):
    """The parallel reads processor needs certain method call sequence."""

    def __init__(self, reason=""):
        super().__init__(reason)


class FileTypeException(Exception):
    """Extension not matching any of those of supported file types."""

    def __init__(self, got, known):
        """
        Declare filetype received and those supported.

        :param str got: file type, name, or path of offending file.
        :param str | Iterable[str] known: supported filetype(s).
        """
        reason = f"'{got}' is not among supported filetypes: {known}"
        super().__init__(reason)


class IllegalChunkException(Exception):
    """Illegal reads chunk ID."""

    def __init__(self, requested, of_interest):
        reason = f"Requested {requested} but processing was restricted to: {of_interest}"
        super().__init__(reason)


class MissingHeaderException(Exception):
    """A reads file undergoing processing must have a header."""

    def __init__(self, filepath=""):
        reason = (
            "No chromosomes in header; this file is "
            "empty or unaligned. Aligned reads are required{}".format(
                f": '{filepath}'" if filepath else "."
            )
        )
        super().__init__(reason)


class MissingOutputFileException(Exception):
    """
    Filepath for particular chunk output doesn't exist.

    Based on its internal settings, the processor's combine() step derives
    a filepath to the output file for each reads chunk indicated by key
    in the argument that it receives. If one of those is missing, it may
    be considered an exceptional case.
    """

    def __init__(self, reads_chunk_key, filepath):
        reason = f"Path to output file for reads chunk '{reads_chunk_key}' does not exist: '{filepath}'"
        super().__init__(reason)


class UnknownChromosomeException(Exception):
    """Represent case in which data about a chromosome is not available."""

    def __init__(self, requested, known=None):
        reason = requested
        if known:
            reason += f"; known: {known}"
        super().__init__(reason)
