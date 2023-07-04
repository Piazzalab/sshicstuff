import os
from os.path import join
import subprocess
import numpy as np
import pandas as pd


base_dir = "/home/nicolas/Documents/Projects/ssHiC"
script = join(base_dir, "hic_ssdna", "sshic", "pipeline.py")
data_dir = join(base_dir, "data", "samples")
inputs_dir = join(data_dir,  "inputs")
refs_dir = join(inputs_dir, "refs")


fragments = join(inputs_dir,  "fragments_list_S288c_DSB_LY_Capture_artificial_DpnIIHinfI.txt")
oligos = join(inputs_dir, "capture_oligo_positions.csv")
additional = join(inputs_dir, "additional_probe_groups.tsv")
centromeres = join(inputs_dir, "S288c_chr_centro_coordinates.tsv")
binning = ["1000", "2000", "3000", "5000", "10000", "20000", "40000", "50000", "80000", "100000"]
ws_centros = 150000
ws_telos = 150000
excluded_chr = ["chr2", "chr3", "2_micron", "mitochondrion", "chr_artificial"]


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


for pcr_type in ["PCRfree", "PCRdupkept"]:
    samples_dir = join(data_dir, pcr_type)
    df_samp2ref: pd.DataFrame = pd.read_csv(join(inputs_dir, f"sample_vs_ref_ponder_{pcr_type}.tsv"), sep="\t")
    for samp in os.listdir(samples_dir):
        samp_name = samp.split(".")[0]
        samp_path = join(samples_dir, samp)
        ref1 = join(refs_dir, df_samp2ref.loc[df_samp2ref["sample"] == samp_name, "reference1"].tolist()[0])
        ref2 = df_samp2ref.loc[df_samp2ref["sample"] == samp_name, "reference2"].tolist()[0]

        ref1_path = join(refs_dir, ref1) + '.tsv'
        run_pipeline(samp_path, ref1_path)
        if ~np.isnan(ref2):
            ref2_path = join(refs_dir, ref2) + '.tsv'
            run_pipeline(samp_path, ref2_path)
