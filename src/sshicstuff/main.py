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
                        This will copy the oligo_capture.csv file and add a new columns with the fragment ID, start and end.

    compare             Compare the capture efficiency of a sample with that of a wild type
                        It may be another sample.

    coverage            Calculate the coverage per fragment and save the result to a bedgraph.
                        Coverage can be normalized and binned at different resolutions.

    dsdnaconly          Keep only Hi-C (dsdna sites) reads from a sparse matrix file (i.e., remove all ssDNA reads).
                        Generate a new sparse matrix file with only dsDNA reads.

    filter              Filter reads from a sparse matrix and keep only pairs of reads that contain at least one
                        oligo/probe (ssdna reads vs whole genome).

    genomaker           Create a chromosome artificial that is the concatenation of the
                        annealing oligos and the enzyme sequence.

    merge               Merge multiple sparse matrix files into a single one.

    pipeline            Run the entire pipeline.
                        It contains the following steps: 
                            - associate
                            - dsdnaconly
                            - ssdnaconly
                            - coverage of dsdna and ssdna reads separately
                            - filter
                            - coverage of all reads at multiple resolutions
                            - profile (4C-like)
                            - stats
                            - rebin
                            - aggregate on centromeric and telomeric regions

    oligo4sshic         generate oligonucleotides for single-strand Hi-C experiments (RUST based sub-module).

    plot4c              Plot a 4C-like profile. Similar graph as those got with the 'view' interactive command (plotly).

    plotmatrix          Plot a contact matrix of contacts made between all the probes. (matplotlib)

    profile             Generate a 4C-like profile for each ssDNA oligo.

    rebin               Rebin change binning resolution of a 4C-like profile

    ssdnaconly          Keep only ssDNA reads from a sparse matrix file (i.e., remove all dsdna reads).
                        Generate a new sparse matrix file with only ssDNA reads.

    stats               Generate statistics and normalization for contacts made by each probe.

    subsample           Subsample and compress FASTQ file using seqtk.

    view                Open a graphical user interface to visualize 4-C like profile (flask + dash + plotly).
                        This will open a web browser with the 4C-like profile u created with the 'profile' command.


"""

import importlib.metadata

from docopt import DocoptExit, docopt

import sshicstuff.commands as commands

__version__ = importlib.metadata.version("sshicstuff")


def main():
    """Main entry point for the sshicstuff CLI."""

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
    except AttributeError as exc:
        print("Unknown command.")
        raise DocoptExit() from exc
    # Create an instance of the command.
    command = command_class(command_args, args)
    # Execute the command.
    command.execute()


if __name__ == "__main__":
    main()
