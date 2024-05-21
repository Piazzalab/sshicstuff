import os
from os.path import join
import pandas as pd
import numpy as np
import base64

from sshicstuff import utils


CWD = os.getcwd()
TEMPORARY_DIRECTORY = join(CWD, "__cache__")


colors_rgba = [
    'rgba(0, 0, 255, 0.8)',  # blue
    'rgba(255, 0, 0, 0.8)',  # red
    'rgba(249, 172, 37, 0.8)',  # yellow
    'rgba(245, 0, 87, 0.8)',  # pink
    'rgba(29, 233, 182, 0.8)',  # green
    'rgba(255, 234, 0, 0.8)',  # yellow 2
    'rgba(255, 11, 0, 0.8)',  # orange
    'rgba(141, 110, 99, 0.8)',  # brown
    'rgba(255, 64, 129, 0.8)',  # pink 2
    'rgba(120, 144, 156, 0.8)',  # blue grey
    'rgba(0, 131, 143, 0.8)',  # cyan
    'rgba(171, 71, 188, 0.8)',  # purple
    'rgba(255, 152, 0, 0.8)',  # amber
    'rgba(0, 150, 136, 0.8)',  # teal
    'rgba(0, 184, 212, 0.8)',  # cyan 2
    'rgba(0, 200, 83, 0.8)',  # green 2
    'rgba(229, 115, 115, 0.8)',  # red 2
    'rgba(255, 167, 38, 0.8)',  # orange 2
    'rgba(61, 90, 254, 0.8)',  # indigo
    'rgba(68, 138, 255, 0.8)',  # blue 2
    'rgba(121, 134, 203, 0.8)',  # deep purple
    'rgba(170, 102, 68, 0.8)',  # deep orange
    'rgba(255, 171, 145, 0.8)',  # pink 3
    'rgba(255, 209, 128, 0.8)'  # amber 2
]

colors_hex = ['#000000', '#0c090a', '#2c3e50', '#34495e', '#7f8c8d', '#8e44ad', '#2ecc71', '#2980b9',
              '#f1c40f', '#d35400', '#e74c3c', '#c0392b', '#1abc9c', '#16a085', '#bdc3c7', '#2c3e50',
              '#7f8c8d', '#f39c12', '#27ae60', '#9b59b6', '#3498db', '#e67e22', '#95a5a6', '#d35400',
              '#f1c40f', '#2980b9', '#e74c3c', '#2ecc71', '#8e44ad', '#34495e', '#1abc9c', '#c0392b',
              '#16a085', '#27ae60', '#7f8c8d', '#f39c12', '#bdc3c7', '#000000', '#0c090a', '#2c3e50',
              '#34495e', '#7f8c8d', '#8e44ad', '#2ecc71', '#2980b9', '#f1c40f', '#d35400', '#e74c3c',
              '#c0392b', '#1abc9c', '#16a085', '#bdc3c7', '#2c3e50', '#7f8c8d', '#f39c12', '#27ae60',
              '#9b59b6', '#3498db', '#e67e22', '#95a5a6', '#d35400', '#f1c40f', '#2980b9', '#e74c3c',
              '#2ecc71', '#8e44ad', '#34495e', '#1abc9c', '#c0392b', '#16a085', '#27ae60', '#7f8c8d',
              '#f39c12', '#bdc3c7', '#000000', '#0c090a', '#2c3e50', '#34495e', '#7f8c8d', '#8e44ad',
              '#2ecc71', '#2980b9', '#f1c40f', '#d35400', '#e74c3c', '#c0392b', '#1abc9c', '#16a085',
              '#bdc3c7', '#2c3e50', '#7f8c8d', '#f39c12', '#27ae60', '#9b59b6', '#3498db', '#e67e22',
              '#95a5a6', '#d35400', '#f1c40f', '#2980b9', '#e74c3c', '#2ecc71', '#8e44ad', '#34495e',
              '#1abc9c', '#c0392b', '#16a085', '#27ae60', '#7f8c8d', '#f39c12', '#bdc3c7']


def save_file(name, content):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(join(TEMPORARY_DIRECTORY, name), "wb") as fp:
        fp.write(base64.decodebytes(data))


def uploaded_files():
    """List the files in the upload directory."""
    files = []
    for filename in os.listdir(TEMPORARY_DIRECTORY):
        path = join(TEMPORARY_DIRECTORY, filename)
        if os.path.isfile(path):
            files.append(filename)
    return files


def rebin_live(
        df: pd.DataFrame,
        bin_size: int,
        df_coords: pd.DataFrame,
):
    """
    Rebin function for the GUI to change resolution of contacts in live mode.
    """

    chr_sizes = dict(zip(df_coords.chr, df_coords.length))
    chr_list, chr_bins = [], []

    for c, l in chr_sizes.items():
        chr_list.append([c] * (l // bin_size + 1))
        chr_bins.append(np.arange(0, (l // bin_size + 1) * bin_size, bin_size))

    chr_list = np.concatenate(chr_list)
    chr_bins = np.concatenate(chr_bins)

    df_template = pd.DataFrame({
        'chr': chr_list,
        'chr_bins': chr_bins,
        'genome_bins': np.arange(0, len(chr_bins)*bin_size, bin_size)
    })

    df["end"] = df["start"] + df["sizes"]
    df["start_bin"] = df["start"] // bin_size * bin_size
    df["end_bin"] = df["end"] // bin_size * bin_size
    df.drop(columns=["genome_start"], inplace=True)

    # chr_filter = df["chr"].unique()
    # if len(chr_filter) == 1:
    #     unique_chr = chr_filter[0]
    #     unique_start = df["start_bin"].min()
    #     unique_end = df["end_bin"].max()
    #     df_template = df_template[
    #         (df_template["chr"] == unique_chr) &
    #         (df_template["chr_bins"] >= unique_start) &
    #         (df_template["chr_bins"] <= unique_end)
    #         ]

    df_cross_bins = df[df["start_bin"] != df["end_bin"]].copy()
    df_in_bin = df.drop(df_cross_bins.index)
    df_in_bin["chr_bins"] = df_in_bin["start_bin"]

    df_cross_bins_a = df_cross_bins.copy()
    df_cross_bins_b = df_cross_bins.copy()
    df_cross_bins_a["chr_bins"] = df_cross_bins["start_bin"]
    df_cross_bins_b["chr_bins"] = df_cross_bins["end_bin"]

    fragments_columns = df.filter(regex='^\d+$').columns.to_list()

    correction_factors = (df_cross_bins_b["end"] - df_cross_bins_b["chr_bins"]) / df_cross_bins_b["sizes"]
    for c in fragments_columns:
        df_cross_bins_a[c] *= (1 - correction_factors)
        df_cross_bins_b[c] *= correction_factors

    df_binned = pd.concat([df_cross_bins_a, df_cross_bins_b, df_in_bin])
    df_binned.drop(columns=["start_bin", "end_bin"], inplace=True)

    df_binned = df_binned.groupby(["chr", "chr_bins"]).sum().reset_index()
    df_binned = utils.sort_by_chr(df_binned, chr_list, 'chr_bins')
    df_binned = pd.merge(df_template, df_binned,  on=['chr', 'chr_bins'], how='left')
    df_binned.drop(columns=["start", "end", "sizes"], inplace=True)
    df_binned.fillna(0, inplace=True)

    return df_binned
