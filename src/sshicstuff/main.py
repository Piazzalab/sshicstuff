#! /usr/bin/env python
# Based on Rémy Greinhofer (rgreinho) tutorial on subcommands in docopt
# https://github.com/rgreinho/docopt-subcommands-example


"""
Single Stranded DNA Hi-C pipeline for generating oligo 4-C profiles and aggregated contact matrices.

usage:
    sshicstuff [-hv] <command> [<args>...]

options:
    -h, --help                  shows the help
    -v, --version               shows the version

The subcommands are:
    aggregate           Aggregate all 4C-like profiles on centromeric or telomeric regions.

    associate           Associate oligo/probe name to fragment/read ID that contains it.

    compare             Compare the capture efficiency of a sample with that of a wild type
                        It may be another sample.

    coverage            Calculate the coverage per fragment and save the result to a bedgraph.

    dsdnaconly          Keep only Hi-C (dsdna sites) reads from a sparse matrix file (i.e., remove all ssDNA reads).
                        Generate a new sparse matrix file with only dsDNA reads.

    filter              Filter reads from a sparse matrix and keep only pairs of reads that contain at least one
                        oligo/probe (ssdna reads vs whole genome).

    genomaker           Create a chromosome artificial that is the concatenation of the
                        annealing oligos and the enzyme sequence.

    pipeline            Run the entire pipeline from filtering to aggregation.

    plot                Plot a 4C-like profile.

    profile             Generate a 4C-like profile for each ssDNA oligo.

    rebin               Rebin change binning resolution of a 4C-like profile

    ssdnaconly          Keep only ssDNA reads from a sparse matrix file (i.e., remove all dsdna reads).
                        Generate a new sparse matrix file with only ssDNA reads.

    stats               Generate statistics and normalization for contacts made by each probe.

    subsample           Subsample and compress FASTQ file using seqtk.

    view                Open a graphical user interface to visualize 4-C like profile


"""

from docopt import docopt
from docopt import DocoptExit
import sshicstuff.commands as commands

import importlib.metadata
__version__ = importlib.metadata.version("sshicstuff")


def main():
    args = docopt(__doc__, version=__version__, options_first=True)
    # Retrieve the command to execute.
    command_name = args.pop("<command>").capitalize()

    # Retrieve the command arguments.
    command_args = args.pop("<args>")
    if command_args is None:
        command_args = {}
    # After 'popping' '<command>' and '<args>', what is left in the
    # args dictionary are the global arguments.

    # Retrieve the class from the 'commands' module.
    try:
        command_class = getattr(commands, command_name)
    except AttributeError:
        print("Unknown command.")
        raise DocoptExit()
    # Create an instance of the command.
    command = command_class(command_args, args)
    # Execute the command.
    command.execute()


if __name__ == "__main__":
    main()
