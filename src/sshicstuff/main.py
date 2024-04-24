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
    gui                 Open a graphical user interface for the pipeline to visualize profiles.

    subsample           Subsample and compress FASTQ file using seqtk.

    genomaker           Create a chromosome artificial that is the concatenation of the
                        annealing oligos and the enzyme sequence.

    associate           Associate oligo/probe name to fragment/read ID that contains it.

    hiconly             Keep only Hi-C reads from a sparse matrix file (i.e., remove all ssDNA reads).

    filter              Filter reads from a sparse matrix and keep only pairs of reads that contain at least one
                        oligo/probe.

    coverage            Calculate the coverage per fragment and save the result to a bedgraph.

    profile             Generate a 4C-like profile for each ssDNA oligo.

    rebin               Rebin change binning resolution of a 4C-like profile

    stats               Generate statistics and normalization for contacts made by each probe.

    aggregate           Aggregate all 4C-like profiles on centromeric or telomeric regions.
"""

from docopt import docopt
from docopt import DocoptExit

import sshicstuff.commands as commands
from sshicstuff.version import __version__


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
