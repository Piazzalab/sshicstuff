import os
import re
from os.path import join
import subprocess
import pandas as pd


def check_nan(str_):
    return str_ != str_


def run_pipeline(sample, reference):
    subprocess.run([
        "python3", script,
        "-s", sample,
        "-f", fragments,
        "-o", oligos,
        "-c", centromeres,
        "-r", reference,
        "-b", *binning,
        "-a", additional,
        "--window-size-centros", str(ws_centros),
        "--window-size-telos", str(ws_telos),
        "--excluded-chr", *excluded_chr,
        "--exclude-probe-chr"],
        check=True)


if __name__ == "__main__":
    pcr_type = "pcrfree"

    base_dir = "/home/nicolas/Documents/Projects/ssdna-hic"
    script = join(base_dir, "sshic", "pipeline.py")
    data_dir = join(base_dir, "data")
    inputs_dir = join(data_dir,  "inputs")
    refs_dir = join(inputs_dir, "refs")
    samples_dir = join(data_dir, "samples", pcr_type.lower())
    samples_only = []

    fragments = join(inputs_dir,  "fragments_list_S288c_DSB_LY_Capture_artificial_DpnIIHinfI.txt")
    oligos = join(inputs_dir, "capture_oligo_positions.csv")
    additional = join(inputs_dir, "additional_probe_groups.tsv")
    centromeres = join(inputs_dir, "S288c_chr_centro_coordinates.tsv")
    binning = ["1000", "2000", "3000", "5000", "10000", "20000", "40000", "50000", "80000", "100000"]
    ws_centros = 150000
    ws_telos = 150000
    excluded_chr = ["chr2", "chr3", "2_micron", "mitochondrion", "chr_artificial"]

    df_samp2ref: pd.DataFrame = pd.read_csv(join(inputs_dir, f"sample_vs_ref_ponder.tsv"), sep="\t")

    sparse_list = sorted([
            file for file in os.listdir(samples_dir) if not os.path.isdir(os.path.join(samples_dir, file))],
        key=lambda x: int(re.search(r'AD(\d+)', x).group(1))
    )

    for samp in sparse_list:
        samp_name = samp.split("_")[0]
        if samp_name not in samples_only and len(samples_only) > 0:
            continue
        samp_path = join(samples_dir, samp)
        ref1 = df_samp2ref.loc[df_samp2ref["sample"] == samp_name, "reference1"].tolist()[0]
        ref2 = df_samp2ref.loc[df_samp2ref["sample"] == samp_name, "reference2"].tolist()[0]

        ref1_path = join(refs_dir, ref1) + '.tsv'
        run_pipeline(samp_path, ref1_path)
        if not check_nan(ref2):
            ref2_path = join(refs_dir, ref2) + '.tsv'
            run_pipeline(samp_path, ref2_path)
