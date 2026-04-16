#! /usr/bin/env python
"""
Single-stranded DNA Hi-C (ssHi-C) analysis toolkit.

Generates probe-centered contact profiles, probe-to-probe interaction
matrices, and genome-wide contact summaries from ssHi-C sparse matrices.

usage:
    sshicstuff [-hv] <command> [<args>...]

options:
    -h, --help     Show this help message.
    -v, --version  Show version.

Available sub-commands
----------------------
    aggregate    Aggregate 4C-like profiles around centromeres or telomeres.
    associate    Map each oligo/probe to its restriction fragment.
    compare      Compare probe capture efficiency against a reference.
    coverage     Compute fragment or bin-level contact coverage.
    design       Design oligos and produce the modified reference genome.
    dsdnaonly    Extract dsDNA-only contacts from a sparse matrix.
    filter       Retain only probe-associated contacts from a sparse matrix.
    merge        Merge multiple sparse matrices by summing contacts.
    pipeline     Run the complete ssHi-C processing workflow.
    plot4c       Generate static 4C-like profile visualizations.
    probe2probe  Build a probe × probe contact matrix.
    profile      Build genome-wide 4C-like probe profiles.
    rebin        Aggregate a 4C-like profile to a larger bin size.
    stats        Compute per-probe contact statistics.
    ssdnaonly    Extract ssDNA-only contacts from a sparse matrix.
    subsample    Subsample a FASTQ file with seqtk.
    view         Launch the interactive profile viewer.
"""

import importlib.metadata

from docopt import DocoptExit, docopt

import sshicstuff.commands as commands

__version__ = importlib.metadata.version("sshicstuff")


def main():
    """Main entry point for the sshicstuff CLI."""
    args = docopt(__doc__, version=__version__, options_first=True)

    command_name = args.pop("<command>").capitalize()
    command_args = args.pop("<args>") or {}

    try:
        command_class = getattr(commands, command_name)
    except AttributeError as exc:
        print(f"Unknown command: {command_name!r}")
        raise DocoptExit() from exc

    command_class(command_args, args).execute()


if __name__ == "__main__":
    main()