""" Test helpers types and functions. """

from pysam import AlignmentFile
from pararead import ParaReadProcessor
from pararead.processor import CORES_PARAM_NAME
from tests import NUM_CORES_DEFAULT

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



class IdentityProcessor(ParaReadProcessor):
    """ Essentially a mock for test cases, simply echoing input. """


    def __init__(self, *args, **kwargs):
        if CORES_PARAM_NAME not in kwargs:
            kwargs[CORES_PARAM_NAME] = NUM_CORES_DEFAULT
        super(IdentityProcessor, self).__init__(*args, **kwargs)

    def __call__(self, chunk, index):
        """
        Simply echo the input, reversing the order for simpler indexing.

        Parameters
        ----------
        chunk : Iterable of object, likely pysam.AlignedSegment
            Chunk of input elements to process.
        index : int or str
            Identifier for the given input chunk.

        Returns
        -------
        (int or str, Iterable of object)
            Same as input pair, but with order reversed.

        """
        return index, chunk



class ReadsfileWrapper(object):
    """ Wrap a pysam reads file for context management. """

    def __init__(self, path_reads_file, reads_ctor=AlignmentFile,
                 **builder_kwargs):
        self.path_reads_file = path_reads_file
        self.reads_ctor = reads_ctor
        self.builder_kwargs = builder_kwargs

    def __enter__(self):
        self.readsfile = self.reads_ctor(self.path_reads_file,
                                         **self.builder_kwargs)
        return self.readsfile

    def __exit__(self, *args):
        try:
            self.readsfile.close()
        except Exception:
            pass