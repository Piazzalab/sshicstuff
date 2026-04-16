#! /usr/bin/env python
# Based on Rémy Greinhofer (rgreinho) tutorial on subcommands in docopt
# https://github.com/rgreinho/docopt-subcommands-example


"""
Single-stranded DNA Hi-C (ssHi-C) analysis toolkit for generating probe-centered contact profiles,
probe-to-probe interaction matrices, and genome-wide contact summaries.

usage:
    sshicstuff [-hv] <command> [<args>...]

options:
    -h, --help                  Show this help message
    -v, --version               Show version

The available subcommands are:

    aggregate
        Aggregate 4C-like profiles around centromeric or telomeric regions to
        produce meta-contact profiles reflecting large-scale chromosomal organization.

    associate
        Map each oligo/probe to its corresponding restriction fragment.
        Adds fragment ID and genomic coordinates to the oligo capture table.

    compare
        Compare probe capture efficiency between a sample and a reference.
        Outputs probe-wise efficiency ratios.

    coverage
        Compute contact coverage per fragment or genomic bin from a sparse matrix.
        Supports normalization and multi-resolution binning.

    dsdnaonly
        Extract dsDNA–dsDNA contacts from a sparse matrix.
        Removes all interactions involving ssDNA-associated fragments.

    ssdnaonly
        Extract ssDNA–ssDNA contacts from a sparse matrix.
        Retains only interactions between probe-associated fragments.

    filter
        Retain only contacts involving at least one probe-associated fragment.
        Enriches capture signal over genome-wide background.

    merge
        Merge multiple sparse matrices by summing contact counts.
        Requires consistent fragment definitions across inputs.

    pipeline
        Run the full ssHi-C processing workflow, including:
            - probe-to-fragment association
            - separation of dsDNA and ssDNA contacts
            - coverage computation
            - probe-based filtering
            - 4C-like profile generation
            - statistical analysis and normalization
            - rebinning
            - aggregation on centromeres and telomeres

    design
        Design oligonucleotides for ssHi-C experiments and generate modified
        reference genomes, including annealing and capture probe tables.

    plot4c
        Generate static visualizations of 4C-like contact profiles.

    profile
        Generate genome-wide 4C-like contact profiles from filtered interactions.

    probe2probe
        Build a probe-by-probe interaction matrix from filtered contacts.
        Optionally export to Cooler format.

    rebin
        Change the resolution of 4C-like profiles by aggregating contacts into
        larger genomic bins.

    stats
        Compute probe-level statistics, including capture efficiency and
        cis/trans contact distributions.

    subsample
        Subsample FASTQ reads using seqtk for dataset size reduction.

    view
        Launch an interactive web interface for exploring 4C-like profiles.
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
