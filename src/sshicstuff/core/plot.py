"""
This module contains the functions to build the figures for the GUI and 
the pipeline plots.
"""

import os
import re

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.io as pio
from plotly.subplots import make_subplots
from matplotlib import pyplot as plt

import sshicstuff.core.methods as methods

# for the colorbar of the plots (hexadecimal colors)
chr_colorbar = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
    '#bcbd22', '#17becf', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94',
    '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5', '#393b79', '#5254a3', '#6b6ecf', '#9c9ede',
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



def build_bins_template(df_coords: pd.DataFrame, bin_size: int) -> pd.DataFrame:
    """
    Build a dataframe template of bins for the whole genome.
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

    fragments_columns = df.filter(regex=r'^\d+$|^\$').columns.to_list()

    correction_factors = (df_cross_bins_b["end"] - df_cross_bins_b["chr_bins"]) / df_cross_bins_b["sizes"]
    for c in fragments_columns:
        df_cross_bins_a[c] *= (1 - correction_factors)
        df_cross_bins_b[c] *= correction_factors

    df_binned = pd.concat([df_cross_bins_a, df_cross_bins_b, df_in_bin])
    df_binned.drop(columns=["start_bin", "end_bin"], inplace=True)

    df_binned = df_binned.groupby(["chr", "chr_bins"]).sum().reset_index()
    df_binned = methods.sort_by_chr(df_binned, chr_list, 'chr_bins')
    df_binned = pd.merge(df_template, df_binned,  on=['chr', 'chr_bins'], how='left')
    df_binned.drop(columns=["start", "end", "sizes"], inplace=True)
    df_binned.fillna(0, inplace=True)

    return df_binned


def colorbar_maker(df_bins: pd.DataFrame):
    """
    Build the color bar for the genome plot.
    """

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
    chr_colors_list = chr_colorbar
    for ii_, chr_ in enumerate(chr_list):
        chr_num = chr_2_num[chr_]
        full_chr_bins.append(f"{chr_num}:{chr_bins[ii_]}")
        x_colors.append(chr_colors_list[chr_num % len(chr_colors_list)])
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
        binsize: int,
        rolling_window: int,
        df_coords: pd.DataFrame,
        df: pd.DataFrame,
        sample_name: str,
        probes: list,
        chr_region: str,
        log_scale: bool,
        user_x_min: str,
        user_x_max: str,
        user_y_min: str,
        user_y_max: str,
        width: int,
        height: int,
):
    
    """
    Build the figure for the GUI-browser.
    """

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
        if chr_region:
            df_bins = df_bins[df_bins['chr'] == chr_region]

        coord_mode[1] = "binned"
        x_min = x_min // binsize * binsize
        x_max = x_max // binsize * binsize + binsize
        df_bins = df_bins[(df_bins['chr_bins'] >= x_min) & (df_bins['chr_bins'] <= x_max)]
        df = rebin_live(df, df_bins, binsize)

        if rolling_window > 1:
            for chr_ in chr_list_unique:
                df.loc[df['chr'] == chr_, probes] = (
                    df.loc[df['chr'] == chr_, probes].rolling(window=rolling_window, min_periods=1).mean())

    y_min = float(user_y_min) if user_y_min else 0.
    y_max = float(user_y_max) if user_y_max else df[probes].max().max()

    if log_scale:
        data = df[probes].values
        data[data == 0] = np.nan
        new_data = np.log10(data)
        y_min = np.nanmin(new_data) if not user_y_max else float(user_y_max)
        y_max = np.nanmax(new_data) if not user_y_min else float(user_y_min)
        df[probes] = new_data

    y_ticks = np.linspace(y_min, y_max, 5)
    y_tick_text = [f"{tick:.3f}" for tick in y_ticks]

    colors_rgba = methods.generate_colors('rgba', len(probes))

    # Making the figure(s)
    if chr_region:
        fig = go.Figure()
        for trace_id, frag in enumerate(probes):
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
                tickformat='d',
                range=[x_min, x_max],
                showgrid=False,
            ),
            yaxis=dict(
                title="Contact frequency",
                tickvals=y_ticks,
                ticktext=y_tick_text,
                range=[y_min, y_max],
                showgrid=False,
            ),
            xaxis_type='linear',
            xaxis_tickformat="d",
            yaxis_tickformat='%.4e' if log_scale else '%.4d',
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

        for trace_id, frag in enumerate(probes):
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
            tickformat='d',
            range=[x_min, x_max],
            showgrid=False,
            ),
            xaxis2=dict(
            tickmode='array',
            tickvals=chr_ticks_pos,
            ticktext=df['chr'].unique(),
            tickfont=dict(size=12),
            ),
            yaxis=dict(
            title="Contact frequency",
            tickvals=y_ticks,
            ticktext=y_tick_text,
            range=[y_min, y_max],
            showgrid=False,
            ),
            yaxis2=dict(
            showticklabels=False,
            ),
            xaxis_type='linear',
            hovermode='closest',
            plot_bgcolor='white',
            paper_bgcolor='white',
            width=width,
            height=height,
        )

    return fig

def plot_profiles(
    profile_contacts_path: str,
    oligo_capture_path: str,
    chr_coord_path: str,
    output_dir: str = None,
    extension: str = "pdf",
    rolling_window: int = 1,
    region: str = None,
    log_scale: bool = False,
    user_y_min: float = None,
    user_y_max: float = None,
    width: int = 1200,
    height: int = 600,
):
    """
    Plot 4C-like profiles.
    """

    profile_type = "contacts"
    if "frequencies" in profile_contacts_path:
        profile_type = "frequencies"

    df: pd.DataFrame = pd.read_csv(profile_contacts_path, sep="\t")
    frags_col = df.filter(regex=r"^\d+$|^\$").columns.to_list()
    df_oligo: pd.DataFrame = pd.read_csv(oligo_capture_path, sep=",")
    probes_to_frag = dict(
        zip(df_oligo["fragment"].astype(str), df_oligo["name"].astype(str))
    )
    df_coords = pd.read_csv(chr_coord_path, sep="\t")

    full_genome_size = df_coords.loc[:, "length"].sum()
    x_min = 0
    x_max = full_genome_size

    b = re.search(r"_(\d+)kb_profile_", profile_contacts_path).group(1)
    if b == 0:
        binsize = 0
        bin_suffix = "0kb"
    else:
        binsize = int(b) * 1000  # kb to bp
        bin_suffix = f"{b}kb"

    sample_name = os.path.basename(profile_contacts_path).split(f"_{bin_suffix}_")[0]

    # Create the output directory
    if not output_dir:
        output_dir = os.path.dirname(profile_contacts_path)

    output_dir = os.path.join(output_dir, "plots")
    output_dir = os.path.join(output_dir, bin_suffix)
    if log_scale:
        output_dir = os.path.join(output_dir, "log")
    else:
        output_dir = os.path.join(output_dir, "raw")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    chr_list_unique = pd.unique(df.chr)

    # Set coordinates preference
    coord_mode = ["genomic", "unbinned"]
    x_col = "genome_start" if binsize == 0 else "genome_bins"
    if region:
        coord_mode[0] = "chromosomal"
        region_list = region.split("-")
        if len(region_list) == 1:
            chr_, start_, end_ = region_list[0], "", ""
        else:
            chr_, start_, end_ = region_list
        max_chr_region_len = df_coords.loc[df_coords.chr == chr_, "length"].values[0]
        x_col = "start" if binsize == 0 else "chr_bins"

        if not start_ and not end_:
            df = df[df["chr"] == chr_]
            df_coords = df_coords[df_coords["chr"] == chr_]
            x_min = 0
            x_max = max_chr_region_len
        else:
            start_, end_ = int(start_), int(end_)
            x_min = start_
            x_max = end_

        df = df[(df["chr"] == chr_) & (df[x_col] >= x_min) & (df[x_col] <= x_max)]
        x_col = "start" if binsize == 0 else "chr_bins"

    if binsize > 0:
        df_bins = build_bins_template(df_coords, binsize)
        if region:
            df_bins = df_bins[df_bins["chr"] == chr_]

        coord_mode[1] = "binned"
        x_min = x_min // binsize * binsize
        x_max = x_max // binsize * binsize + binsize
        df_bins = df_bins[
            (df_bins["chr_bins"] >= x_min) & (df_bins["chr_bins"] <= x_max)
        ]

        if rolling_window > 1:
            for chr_ in chr_list_unique:
                df.loc[df["chr"] == chr_, frags_col] = (
                    df.loc[df["chr"] == chr_, frags_col]
                    .rolling(window=rolling_window, min_periods=1)
                    .mean()
                )

    y_min = float(user_y_min) if user_y_min else 0.0
    y_max = float(user_y_max) if user_y_max else df[frags_col].max().max()

    log_status = ""
    if log_scale:
        log_status = "log_"
        data = df[frags_col].values
        data[data == 0] = np.nan
        new_data = np.log10(data)
        y_min = np.nanmin(new_data) if not user_y_max else float(user_y_max)
        y_max = np.nanmax(new_data) if not user_y_min else float(user_y_min)
        df[frags_col] = new_data

    y_ticks = np.linspace(y_min, y_max, 5)
    y_tick_text = [f"{tick:.3f}" for tick in y_ticks]

    colors_rgba = methods.generate_colors("rgba", len(frags_col))

    if region:
        for ii_f, frag in enumerate(frags_col):
            fig = go.Figure()
            probe = probes_to_frag.get(frag, "")
            output_path = os.path.join(
                output_dir,
                f"{sample_name}_{frag}_{probe}_{profile_type}_{bin_suffix}_{log_status}{region}.{extension}",
            )

            fig.add_trace(
                go.Scattergl(
                    x=df[x_col],
                    y=df[frag],
                    name=frag,
                    mode="lines",
                    line=dict(width=1, color=colors_rgba[ii_f]),
                    marker=dict(size=4),
                    showlegend=False,
                )
            )

            fig.update_layout(
                title=f"{sample_name}",
                xaxis=dict(
                    title=dict(text=f"{region} coordinates", standoff=80),
                    tickformat="d",
                    range=[x_min, x_max],
                    showgrid=False,
                ),
                yaxis=dict(
                    title=f"{profile_type.capitalize()}{' - log' if log_scale else ''}",
                    tickvals=y_ticks,
                    ticktext=y_tick_text,
                    range=[y_min, y_max],
                    showgrid=False,
                ),
                xaxis_type="linear",
                xaxis_tickformat="d",
                yaxis_tickformat="%.4e" if log_scale else "%.4d",
                plot_bgcolor="white",
                paper_bgcolor="white",
                width=width,
                height=height,
            )

            pio.write_image(fig, output_path, engine="kaleido")

    else:
        df_10kb_tmp = build_bins_template(df_coords=df_coords, bin_size=10000)
        colorbar, chr_ticks_pos = colorbar_maker(df_10kb_tmp)

        # plot for each prob or group of probes
        for ii_f, frag in enumerate(frags_col):
            probe = probes_to_frag.get(frag, "")
            output_path = os.path.join(
                output_dir,
                f"{sample_name}_{frag}_{probe}_{profile_type}_{bin_suffix}_{log_status}.{extension}",
            )

            fig = make_subplots(
                rows=2,
                cols=1,
                row_heights=[0.94, 0.06],
                vertical_spacing=0.06,
                specs=[[{"type": "scatter"}], [{"type": "bar"}]],
            )

            fig.add_trace(
                go.Scatter(
                    x=df[x_col],
                    y=df[frag],
                    mode="lines",
                    line=dict(width=1, color=colors_rgba[ii_f]),
                    marker=dict(size=4),
                    showlegend=False,
                ),
                row=1,
                col=1,
            )

            fig.add_trace(colorbar, row=2, col=1)
            title = f" {sample_name} <br> {frag} {'- ' + probe}"
            fig.update_layout(
                title=title,
                xaxis=dict(
                    title=dict(text="Genomic coordinates", standoff=80),
                    tickformat="d",
                    range=[x_min, x_max],
                    showgrid=False,
                ),
                xaxis2=dict(
                    tickmode="array",
                    tickvals=chr_ticks_pos,
                    ticktext=df["chr"].unique(),
                    tickfont=dict(size=12),
                ),
                yaxis=dict(
                    title=f"{profile_type.capitalize()}{' - log' if log_scale else ''}",
                    tickvals=y_ticks,
                    ticktext=y_tick_text,
                    range=[y_min, y_max],
                    showgrid=False,
                ),
                yaxis2=dict(
                    showticklabels=False,
                ),
                xaxis_showgrid=False,
                yaxis_showgrid=False,
                xaxis_type="linear",
                xaxis_tickformat="d",
                yaxis_tickformat="%.4e" if log_scale else "%.4d",
                xaxis_range=[x_min, x_max],
                yaxis_range=[y_min, y_max],
                hovermode="closest",
                plot_bgcolor="white",
                paper_bgcolor="white",
                width=width,
                height=height,
            )

            pio.write_image(fig, output_path, engine="kaleido")



def plot_probes_matrix(
        probes_matrix_path: str,
        probes_a: list = None,
        probes_b: list= None,
        extension: str = 'png',
        output_dir: str = None,
        title: str = None,
        logscale: bool = False,
        cmap: str = 'Reds',
        vmin: int = 0,
        vmax: int = None,
        dpi: int = 300,
        fsize: int = 14
) -> None:
    
    fontsize = 12
    if not title:
        title = os.path.basename(probes_matrix_path)

    # Check if the probes matrix file exists
    methods.check_if_exists(probes_matrix_path)

    # Set output directory to the input file's directory if not provided
    if not output_dir:
        output_dir = os.path.dirname(probes_matrix_path)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_path = os.path.join(output_dir, f"probes_matrix_heatmap.{extension}")

    # Read the probes matrix into a DataFrame (first column as index)
    df = pd.read_csv(probes_matrix_path, sep="\t", index_col=0)

    # Optionally filter the DataFrame by specified probe lists
    if probes_a:
        df = df.loc[probes_a, :]
    if probes_b:
        df = df.loc[:, probes_b]

    # Get the data array from the DataFrame
    array = df.values
    if logscale:
        array = np.log10(array + 1)

    # If vmax is not provided, set it to the 99th percentile of the data
    if vmax is None:
        vmax = np.percentile(array, 99)

    # Set image display parameters
    im_kwargs = {'vmin': vmin, 'vmax': vmax, 'cmap': cmap, 'interpolation': "none"}

    # Create the heatmap
    plt.figure()
    # figsize must be squared if n row == n col
    if df.shape[0] == df.shape[1]:
        plt.figure(figsize=(fsize, fsize))
    # else must me proportional to the number of rows and columns
    elif df.shape[0] > df.shape[1]:
     r = df.shape[0] / df.shape[1]
     plt.figure(figsize=(fsize, fsize * r))
    else:
        r = df.shape[1] / df.shape[0]
        plt.figure(figsize=(fsize * r, fsize))

    
    img = plt.imshow(array, **im_kwargs)
    cbar = plt.colorbar(img)
    if logscale:
        cbar_label = "Log10(Frequencies + 1)"
    else:
        cbar_label = "Frequencies"
    cbar.set_label(cbar_label, fontsize=fontsize+2)
    cbar.ax.tick_params(labelsize=fontsize)

    # Set tick labels with rotation and font size adjustments
    plt.xticks(ticks=np.arange(len(df.columns)), labels=df.columns, rotation=65, fontsize=fontsize, ha='right')
    plt.yticks(ticks=np.arange(len(df.index)), labels=df.index, fontsize=fontsize)

    # Add axis labels
    plt.xlabel("Probes", fontsize=fontsize + 2)
    plt.ylabel("Probes", fontsize=fontsize + 2)

    # Set the title if provided
    if title:
        plt.title(title, fontsize=fontsize + 4)

    # Save the figure and close
    plt.savefig(output_path, bbox_inches="tight", pad_inches=0.1, dpi=dpi)
    plt.close()