import os
from os.path import join, dirname
import sshicstuff.commands as shcmd


CWD = os.getcwd()
TESTDIR = join(dirname(CWD), "test_data")
CAPTURE = join(TESTDIR, "capture_oligo_positions.csv")
CAPTURE_FRAGMENTS = CAPTURE.replace(".csv", "_fragments_associated.csv")
ANNEALING = join(TESTDIR, "annealing_oligo_positions.csv")
COORDS = join(TESTDIR, "chr_coords.tsv")
GENOME = join(TESTDIR, "S288c_DSB_LY_Capture.fa")
DPNII = "gatc"
SEED = 42


def test_subsample():
    args = (
        "-i {0} -c -F -n 10000 -s {1}"
    ).format(join(TESTDIR, "AD162/AD162_sub100K.end1.fastq.gz"), SEED)
    proc = shcmd.Subsample(args.split(" "), {})
    proc.execute()


def test_associate():
    args = (
        "-f {0} -o {1} -F"
    ).format(join(TESTDIR, "AD162/fragments_list.txt"), CAPTURE)
    proc = shcmd.Associate(args.split(" "), {})
    proc.execute()


def test_hiconly():
    args = (
        "-c {0} -m {1} -n 3 -F"
    ).format(CAPTURE_FRAGMENTS, join(TESTDIR, "AD162/AD162_pcrdupkept.txt"))
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
    args = (
        "-c {0} -C {1} -f {2} -a {3} -F -N"
    ).format(
        CAPTURE,
        COORDS,
        join(TESTDIR, "AD162/AD162_pcrdupkept_filtered.tsv"),
        join(TESTDIR, "additional_probe_groups.tsv")
    )
    proc = shcmd.Profile(args.split(" "), {})
    proc.execute()


def test_rebin():
    args = (
        "-b 10000 -c {0} -p {1} -F"
    ).format(
        COORDS,
        join(TESTDIR, "AD162/AD162_pcrdupkept_0kb_profile_frequencies.tsv"),
    )
    proc = shcmd.Rebin(args.split(" "), {})
    proc.execute()


def test_stats():
    args = (
        "-c {0} -C {1} -m {2} -p {3} -F -r 50000"
    ).format(
        CAPTURE,
        COORDS,
        join(TESTDIR, "AD162/AD162_pcrdupkept.txt"),
        join(TESTDIR, "AD162/AD162_pcrdupkept_0kb_profile_frequencies.tsv")
    )

    proc = shcmd.Stats(args.split(" "), {})
    proc.execute()


def test_aggregate():
    args = (
        "-c {0} -h {1} -p {2} -T {3} -w 150000 -N -I -L"
    ).format(
        CAPTURE,
        COORDS,
        join(TESTDIR, "AD162/AD162_pcrdupkept_10kb_profile_frequencies.tsv"),
        "-E chr3 -E chr2 -E 2_micron -E mitochondrion -E chr_artificial_donor -E chr_artificial_ssDNA"
    )
    proc = shcmd.Aggregate(args.split(" "), {})
    proc.execute()


