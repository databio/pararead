#!/usr/bin/env python
""" Counting reads, as an example/template for implementing a processor. """

import argparse
import sys

from pararead import add_logging_options, ParaReadProcessor

__author__ = "Vince Reuter"
__email__ = "vince.reuter@gmail.com"


def _parse_cmdl(cmdl):
    """ Define and parse command-line interface. """

    parser = argparse.ArgumentParser(
        description="Read count as template for ParaReadProcessor "
                    "implementation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "readsfile", help="Path to sequencing reads file.")

    parser.add_argument(
        "-O", "--outfile", required=True, help="Path to output file.")
    parser.add_argument(
        "-C", "--cores", required=False, default=1, help="Number of cores.")

    parser = add_logging_options(parser)
    return parser.parse_args(cmdl)


class ReadCounter(ParaReadProcessor):
    """ Sequencing reads counter. """
    def __call__(self, chromosome, _=None):
        n_reads = self.readsfile.count(chromosome)
        with open(self._tempf(chromosome), 'w') as f:
            f.write("{}\t{}".format(chromosome, n_reads))
        return chromosome


def main(cmdl):
    """ Run the script. """
    from pararead import logger_via_cli

    opts = _parse_cmdl(cmdl)
    logger = logger_via_cli(opts)

    logger.debug("Creating counter")
    counter = ReadCounter(opts.readsfile, cores=opts.cores,
                          outfile=opts.outfile, action="CountReads")
    logger.debug("Registering files")
    counter.register_files()

    logger.info("Counting reads: {}".format(opts.readsfile))
    good_chromosomes = counter.run()
    logger.info("Collecting read counts: {}".format(opts.outfile))
    counter.combine(good_chromosomes, chrom_sep="\n")


if __name__ == "__main__":
    main(sys.argv[1:])
