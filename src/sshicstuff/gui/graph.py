import pandas as pd
import numpy as np

import plotly.graph_objs as go
from plotly.subplots import make_subplots

from sshicstuff.gui.common import sort_by_chr
from sshicstuff.gui.common import chr_to_exclude
from sshicstuff.gui.colors import colors_hex, colors_rgba, chr_colors


"""
###################
     METHODS
###################
"""


def transform_data(data: np.array, y_max: float, user_y_max: float, y_min: float):
    if data.max().max() <= 1.:
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

    return new_data, y_max, y_min, re_scale_output


def build_bins_template(df_coords: pd.DataFrame, bin_size: int) -> pd.DataFrame:
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
    return df_template


def rebin_live(df: pd.DataFrame, df_template: pd.DataFrame, bin_size: int):
    """
    Rebin function for the GUI to change resolution of contacts in live mode.
    """

    chr_list = df_template["chr"]
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

    fragments_columns = df.filter(regex='^\d+$|^\$').columns.to_list()

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


def colorbar_maker(df_bins: pd.DataFrame):
    x_colors = []
    chr_ticks = []
    chr_ticks_pos = []
    prev_chr = None
    chr_start = 0

    chr_list = df_bins.chr.values
    chr_bins = df_bins.chr_bins.values
    chr_2_num = {c: i for i, c in enumerate(pd.unique(df_bins.chr))}
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


def figure_maker(
        relayout_data: dict,
        binsize: int,
        df_coords: pd.DataFrame,
        df: pd.DataFrame,
        sample_name: str,
        probes: list,
        chr_region: str,
        rescale: bool,
        user_x_min: str,
        user_x_max: str,
        user_y_min: str,
        user_y_max: str,
        width: int,
        height: int,
):

    chr_list = df.chr.values
    chr_list_unique = pd.unique(df.chr)
    n_chr = len(chr_list_unique)

    # Set coordinates preference
    full_genome_size = df_coords.loc[n_chr-1, 'cumu_start'] + df_coords.loc[n_chr-1, 'length']
    x_min = 0
    x_max = full_genome_size

    coord_mode = ["genomic", "unbinned"]
    x_col = "genome_start" if binsize == 0 else "genome_bins"
    if chr_region:
        coord_mode[0] = "chromosomal"
        max_chr_region_len = df_coords.loc[df_coords.chr == chr_region, 'length'].values[0]
        x_min = float(user_x_min) if user_x_min else 0
        x_max = float(user_x_max) if user_x_max else max_chr_region_len
        df = df[df['chr'] == chr_region]
        df = df[(df['start'] >= x_min) & (df['start'] <= x_max)]
        x_col = "start" if binsize == 0 else "chr_bins"

    if binsize > 0:
        df_bins = build_bins_template(df_coords, binsize)
        chr_bins = df_bins.chr_bins.values
        genome_bins = df_bins.genome_bins.values
        n_bins = len(chr_bins)
        if chr_region:
            df_bins = df_bins[df_bins['chr'] == chr_region]

        coord_mode[1] = "binned"
        x_min = x_min // binsize * binsize
        x_max = x_max // binsize * binsize + binsize
        df_bins = df_bins[(df_bins['chr_bins'] >= x_min) & (df_bins['chr_bins'] <= x_max)]
        df = rebin_live(df, df_bins, binsize)

    y_min = float(user_y_min) if user_y_min else 0.
    y_max = float(user_y_max) if user_y_max else df[probes].max().max()

    rescale_output = ""
    if rescale:
        data = df[probes].values
        new_data, y_max, y_min, rescale_output = transform_data(data, y_max, user_y_max, y_min)
        df[probes] = new_data

    # Making the figure(s)
    if chr_region:
        fig = go.Figure()
        for j in range(len(probes)):
            frag = probes[j]
            trace_id = j
            fig.add_trace(
                go.Scattergl(
                    x=df[x_col],
                    y=df[frag],
                    name=frag,
                    mode='lines',
                    line=dict(width=1, color=colors_rgba[trace_id]),
                    marker=dict(size=4)
                )
            )
        fig.update_layout(
            title=f"{sample_name}",
            xaxis=dict(
                title=f"{chr_region} coordinates",
                range=[x_min, x_max]
            ),
            yaxis=dict(
                title="Contact frequency",
                range=[y_min, y_max]
            ),
            xaxis_type='linear',
            xaxis_tickformat="d",
            plot_bgcolor='white',
            paper_bgcolor='white',
            width=width,
            height=height,
        )

    else:
        # Make the genome color bar per chromosome
        df_10kb = build_bins_template(df_coords, 10000)
        colorbar, chr_ticks_pos = colorbar_maker(df_10kb)

        fig = make_subplots(
            rows=2, cols=1, row_heights=[0.94, 0.06], vertical_spacing=0.06,
            specs=[[{'type': 'scatter'}], [{'type': 'bar'}]]
        )

        for j in range(len(probes)):
            frag = probes[j]
            trace_id = j
            fig.add_trace(
                go.Scattergl(
                    x=df[x_col],
                    y=df[frag],
                    name=frag,
                    mode='lines',
                    line=dict(width=1, color=colors_rgba[trace_id]),
                    marker=dict(size=4)
                ),
                row=1, col=1
            )
            fig.add_trace(colorbar, row=2, col=1)

            fig.update_layout(
                title=f"{sample_name}",
                xaxis=dict(
                    title=dict(text="Genomic coordinates", standoff=80),
                    range=[x_min, x_max]
                ),
                xaxis2=dict(
                    tickmode='array',
                    tickvals=chr_ticks_pos,
                    ticktext=df['chr'].unique(),
                    tickfont=dict(size=12),
                ),
                yaxis=dict(
                    title="Contact frequency"
                ),
                yaxis2=dict(
                    showticklabels=False,
                ),

                xaxis_showgrid=False,
                yaxis_showgrid=False,
                xaxis_type='linear',
                xaxis_tickformat="d",
                xaxis_range=[x_min, x_max],
                yaxis_range=[y_min, y_max],
                hovermode='closest',
                plot_bgcolor='white',
                paper_bgcolor='white',
                width=width,
                height=height,
            )

    return fig, rescale_output
