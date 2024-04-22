import os
import pytest
from os.path import join, dirname
import sshicstuff.commands as shcmd


CWD = os.getcwd()
TESTDIR = join(dirname(CWD), "test_data")
CAPTURE = join(TESTDIR, "capture_oligo_positions.csv")
ANNEALING = join(TESTDIR, "annealing_oligo_positions.csv")
COORDS = join(TESTDIR, "chr_coords.csv")
GENOME = join(TESTDIR, "S288c_DSB_LY_Capture.fa")
DPNII = "gatc"
SEED = 42


def test_subsample():
    args = (
        "-i {0} -c -F -n 10000 -s {1}"
    ).format(join(TESTDIR, "AD162/AD162_sub100K.end1.fastq.gz"), SEED)
    proc = shcmd.Subsample(args.split(" "), {})
    proc.execute()


def test_genomaker():
    add_chr_donor = join(TESTDIR, "chr_artificial_donor.fa")
    args = (
        "-e {0} -g {1} -o {2} -a {3} -s N -l 80 -f 150"
    ).format(DPNII, GENOME, ANNEALING, add_chr_donor)
    proc = shcmd.Genomaker(args.split(" "), {})
    proc.execute()


def test_associate():
    args = (
        "-f {0} -o {1} -F"
    ).format(join(TESTDIR, "AD162/fragments_list.txt"), CAPTURE)
    proc = shcmd.Associate(args.split(" "), {})
    proc.execute()


def test_hiconly():
    args = (
        "-c {0} -m {1} -n 2 -F"
    ).format(CAPTURE, join(TESTDIR, "AD162/AD162_pcrdupkept.txt"))
    proc = shcmd.Hiconly(args.split(" "), {})
    proc.execute()


def test_filter():
    args = (
        "-c {0} -f {1} -m {2} -F"
    ).format(
        CAPTURE,
        join(TESTDIR, "AD162/fragments_list.txt"),
        join(TESTDIR, "AD162/AD162_pcrdupkept.txt")
        )
    proc = shcmd.Filter(args.split(" "), {})
    proc.execute()


def test_coverage():
    args = (
        "-f {0} -m {1} -F -N"
    ).format(
        join(TESTDIR, "AD162/fragments_list.txt"),
        join(TESTDIR, "AD162/AD162_pcrdupkept.txt")
    )
    proc = shcmd.Coverage(args.split(" "), {})
    proc.execute()


def test_profile():
    pass


def test_rebin():
    pass


def test_stats():
    pass


def test_aggregate():
    pass



