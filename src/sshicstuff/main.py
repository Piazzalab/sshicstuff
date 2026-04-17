#! /usr/bin/env python
"""
Single-stranded DNA Hi-C (ssHi-C) analysis toolkit.

The toolkit is now cool-first: it consumes fragment-level ``.cool``
files for all 2-D matrix operations, while still accepting the legacy
GRAAL triplet in the ``pipeline`` command for backward compatibility.

usage:
    sshicstuff [-hv] <command> [<args>...]

options:
    -h, --help     Show this help message.
    -v, --version  Show version.

Available sub-commands
----------------------
    aggregate    Aggregate 4C-like profiles around centromeres or telomeres.
    associate    Map each oligo/probe to its restriction fragment.
    balance      ICE-balance a fragment-level cool file.
    compare      Compare probe capture efficiency against a reference.
    coverage     Compute fragment or bin-level contact coverage from a cool.
    design       Design oligos and produce the modified reference genome.
    dsdnaonly    Extract dsDNA-only contacts from a cool.
    filter       Retain only probe-associated contacts from a cool.
    graal2cool   Convert a legacy graal sparse matrix to fragment-level cool.
    merge        Merge multiple fragment-aligned cool files.
    pipeline     Run the complete ssHi-C processing workflow.
    plot4c       Generate static 4C-like profile visualizations.
    probe2probe  Build a probe × probe contact matrix.
    profile      Build genome-wide 4C-like probe profiles.
    rebin        Aggregate a 4C-like profile to a larger bin size.
    stats        Compute per-probe contact statistics.
    ssdnaonly    Extract ssDNA-only contacts from a cool.
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