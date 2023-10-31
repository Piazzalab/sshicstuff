import os
import re
from os.path import join, dirname
import subprocess
import pandas as pd


def check_nan(str_):
    return str_ != str_


def run_pipeline(sample, reference):
    command = [
        "python3", pipeline,
        "-s", sample,
        "-f", fragments,
        "-o", oligos,
        "-c", centromeres,
        "-b", *binning,
        "-a", additional,
        "--window-size-centros", str(ws_centros),
        "--window-size-telos", str(ws_telos),
        "--excluded-chr", *excluded_chr,
        "--exclude-probe-chr"
    ]

    if reference is not None:
        command.extend(["-r", reference])

    subprocess.run(command, check=True)


if __name__ == "__main__":
    cwd = os.getcwd()
    pipeline = join(cwd, "pipeline.py")

    base_dir = dirname(cwd)
    data_dir = join(base_dir, "data")
    inputs_dir = join(data_dir,  "inputs")
    smaplesheet = join(inputs_dir, "samplesheet.csv")
    refs_dir = join(data_dir, "references")
    samples_dir = join(data_dir, "samples")

    # parameters for pipeline
    fragments = join(inputs_dir,  "fragments_list_S288c_DSB_LY_Capture_artificial_DpnIIHinfI.txt")
    oligos = join(inputs_dir, "capture_oligo_positions.csv")
    additional = join(inputs_dir, "additional_probe_groups.tsv")
    centromeres = join(inputs_dir, "S288c_chr_centro_coordinates.tsv")
    binning = ["1000", "2000", "3000", "5000", "10000", "20000", "40000", "50000", "80000", "100000"]
    ws_centros = 150000
    ws_telos = 150000
    excluded_chr = ["chr2", "chr3", "2_micron", "mitochondrion", "chr_artificial"]

    df_samplesheet: pd.DataFrame = pd.read_csv(smaplesheet, sep=",")
    samples = {}
    for _, row in df_samplesheet.iterrows():
        samples[row.loc["name"]] = []
        if len(row) > 1:
            for i in range(1, len(row)):
                if not check_nan(row.iloc[i]):
                    samples[row.loc["name"]].append(row.iloc[i])

    sparse_list = sorted([
            file for file in os.listdir(samples_dir) if not os.path.isdir(os.path.join(samples_dir, file))],
        key=lambda x: int(re.search(r'AD(\d+)', x).group(1))
    )

    for samp in sparse_list:
        for name in samples:
            if name in samp:
                samp_name = name
                break
        else:
            raise ValueError(f"Sample {samp} not found in samplesheet")

        samp_path = join(samples_dir, samp)
        if samples[samp_name]:
            for ref_name in samples[samp_name]:
                reference_path = join(refs_dir, ref_name, '.tsv')
                run_pipeline(samp_path, reference_path)
        else:
            run_pipeline(samp_path, None)
