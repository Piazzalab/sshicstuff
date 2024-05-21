import os
import re
from os.path import join
import pandas as pd
import numpy as np
import base64

# plotly
from plotly import graph_objs as go


"""
###################
VARIABLES TO IMPORT
###################
"""

CWD = os.getcwd()
TEMPORARY_DIRECTORY = join(CWD, "__cache__")

chr_colors = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
    '#bcbd22', '#17becf', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94',
    '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5', '#393b79', '#5254a3', '#6b6ecf', '#9c9ede',
]

chr_to_exclude = ["chr_artificial_donor", "chr_artificial_ssDNA"]

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

empty_figure = go.Figure(
    layout=go.Layout(
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        annotations=[
            dict(
                x=0.5,
                y=0.5,
                text="No data available",
                showarrow=False,
                font=dict(size=28)
            )
        ],
        hovermode='closest',
        plot_bgcolor='white',
        paper_bgcolor='white'
    )
)


"""
#################
METHODS TO IMPORT
#################
"""


def build_chr_bins_lists(df_coords: pd.DataFrame, bin_size: int):
    chr_sizes = dict(zip(df_coords.chr, df_coords.length))
    chr_list, chr_bins = [], []

    for c, l in chr_sizes.items():
        chr_list.append([c] * (l // bin_size + 1))
        chr_bins.append(np.arange(0, (l // bin_size + 1) * bin_size, bin_size))

    chr_list = np.concatenate(chr_list)
    chr_bins = np.concatenate(chr_bins)
    return chr_list, chr_bins


def make_colorbar(df_coords: pd.DataFrame, bins_size: int):
    x_colors = []
    chr_ticks = []
    chr_ticks_pos = []
    prev_chr = None
    chr_start = 0

    chr_list, chr_bins = build_chr_bins_lists(df_coords, bins_size)
    chr_2_num = {c: i for i, c in enumerate(pd.unique(df_coords.chr))}
    n_bins = len(chr_bins)

    full_chr_bins = []
    for ii_, chr_ in enumerate(chr_list):
        chr_num = chr_2_num[chr_]
        full_chr_bins.append(f"{chr_num}:{chr_bins[ii_]}")
        x_colors.append(chr_colors[chr_num % len(colors_hex)])
        if chr_num != prev_chr:
            if prev_chr is not None:
                chr_ticks_pos.append((ii_ + chr_start) / 2)
            chr_ticks.append(chr_num)
            chr_start = ii_
            prev_chr = chr_num
    chr_ticks_pos.append((n_bins + chr_start) / 2)

    color_bar = go.Bar(
        x=full_chr_bins,
        y=[1] * n_bins,
        marker=dict(color=x_colors, line=dict(width=0)),
        showlegend=False,
        hoverinfo='none'
    )

    return color_bar, chr_ticks_pos


def rebin_live(df: pd.DataFrame, bin_size: int, df_coords: pd.DataFrame):
    """
    Rebin function for the GUI to change resolution of contacts in live mode.
    """

    chr_list, chr_bins = build_chr_bins_lists(df_coords, bin_size)

    df_template = pd.DataFrame({
        'chr': chr_list,
        'chr_bins': chr_bins,
        'genome_bins': np.arange(0, len(chr_bins)*bin_size, bin_size)
    })

    df["end"] = df["start"] + df["sizes"]
    df["start_bin"] = df["start"] // bin_size * bin_size
    df["end_bin"] = df["end"] // bin_size * bin_size
    df.drop(columns=["genome_start"], inplace=True)

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
    df_binned = sort_by_chr(df_binned, chr_list, 'chr_bins')
    df_binned = pd.merge(df_template, df_binned,  on=['chr', 'chr_bins'], how='left')
    df_binned.drop(columns=["start", "end", "sizes"], inplace=True)
    df_binned.fillna(0, inplace=True)

    return df_binned


def save_file(name, content):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(join(TEMPORARY_DIRECTORY, name), "wb") as fp:
        fp.write(base64.decodebytes(data))


def sort_by_chr(df: pd.DataFrame, chr_list: list[str], *args: str):
    chr_list = np.unique(chr_list)
    chr_with_number = [c for c in chr_list if re.match(r'chr\d+', c)]
    chr_with_number.sort(key=lambda x: int(x[3:]))
    chr_without_number = [c for c in chr_list if c not in chr_with_number]

    order = chr_with_number + chr_without_number
    df['chr'] = df['chr'].apply(lambda x: order.index(x) if x in order else len(order))

    if args:
        df = df.sort_values(by=['chr', *args])
    else:
        df = df.sort_values(by=['chr'])

    df['chr'] = df['chr'].map(lambda x: order[x])
    df.index = range(len(df))
    return df


def transform_data(data: np.array, y_max: float, user_y_max: float, y_min: float, re_scale: bool):
    re_scale_output = ""
    if re_scale:
        if y_max <= 1.:
            # squared root transformation
            new_data = np.sqrt(data + 1e-8)
            y_max = np.sqrt(y_max) if not user_y_max else user_y_max
            y_min = np.sqrt(y_min) if y_min > 0 else 0
            re_scale_output = "sqrt"
        else:
            # log transformation
            new_data = np.log(data + 1)
            y_max = np.log(y_max) if not user_y_max else user_y_max
            y_min = 0
            re_scale_output = "log"
    else:
        new_data = data

    return new_data, y_max, y_min, re_scale_output


def uploaded_files():
    """List the files in the upload directory."""
    files = []
    for filename in os.listdir(TEMPORARY_DIRECTORY):
        path = join(TEMPORARY_DIRECTORY, filename)
        if os.path.isfile(path):
            files.append(filename)
    return files
