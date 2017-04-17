""" Specific exception types. """

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



class ExecutionOrderException(Exception):
    """ The parallel reads processor needs certain method call sequence. """
    def __init__(self, reason=""):
        super(ExecutionOrderException, self).__init__(reason)



class FileTypeException(Exception):
    """ Extension not matching any of those of supported file types. """

    def __init__(self, got, known):
        """
        Declare filetype received and those supported.

        Parameters
        ----------
        got : str
            File type, name, or path of offending file.
        known : str | collections.Iterable of str
            Supported filetype(s).
        """
        reason = "'{}' is not among supported filetypes: {}".format(got, known)
        super(FileTypeException, self).__init__(reason)



class MissingHeaderException(Exception):
    """ A reads file undergoing processing must have a header. """
    def __init__(self, filepath=""):
        reason = "No chromosomes in header; this file is " \
                 "empty or unaligned. Aligned reads are required{}".\
                 format(": '{}'".format(filepath) if filepath else ".")
        super(MissingHeaderException, self).__init__(reason)
