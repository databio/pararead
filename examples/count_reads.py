#!/usr/bin/env python
"""Counting reads, as an example/template for implementing a processor."""

import argparse
import logging
import sys

from pararead import ParaReadProcessor

__author__ = "Vince Reuter"
__email__ = "vince.reuter@gmail.com"

_LOGGER = logging.getLogger(__name__)


def _parse_cmdl(cmdl):
    """Define and parse command-line interface."""

    parser = argparse.ArgumentParser(
        description="Read count as template for ParaReadProcessor implementation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("readsfile", help="Path to sequencing reads file.")

    parser.add_argument("-O", "--outfile", required=True, help="Path to output file.")

    parser.add_argument("-C", "--cores", required=False, default=1, help="Number of cores.")

    parser.add_argument(
        "-t",
        "--limit",
        dest="limit",
        help="Limit to these chromosomes",
        nargs="+",
        default=None,
    )

    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")

    return parser.parse_args(cmdl)


class ReadCounter(ParaReadProcessor):
    """Sequencing reads counter."""

    def __call__(self, chromosome, _=None):
        # doesn't work because count can't handle multiple iterators
        # n_reads = self.readsfile.count(chromosome)

        # Here we want to go through each read and do what we need.
        # In this simple example, all we're doing is counting them.
        n_reads = 0
        reads = self.fetch_chunk(chromosome)
        for read in reads:
            n_reads += 1

        _LOGGER.debug(f"Chromosome: '{chromosome}'; n_reads {n_reads}")
        with open(self._tempf(chromosome), "w") as f:
            f.write(f"{chromosome}\t{n_reads}")
        return chromosome


def main(cmdl):
    """Run the script."""

    args = _parse_cmdl(cmdl)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    _LOGGER.debug("Creating counter")
    counter = ReadCounter(
        args.readsfile,
        cores=args.cores,
        outfile=args.outfile,
        action="CountReads",
        limit=args.limit,
    )
    _LOGGER.debug("Registering files")
    counter.register_files()

    _LOGGER.info(f"Counting reads: {args.readsfile}")
    good_chromosomes = counter.run()
    _LOGGER.info(f"Collecting read counts: {args.outfile}")
    counter.combine(good_chromosomes, chrom_sep="\n")


if __name__ == "__main__":
    main(sys.argv[1:])
