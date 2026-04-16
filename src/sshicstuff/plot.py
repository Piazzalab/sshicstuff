"""
Plotting functions for ssHi-C profiles, probe matrices, and GUI figures.

This module is intentionally isolated from data-production logic: every
function here takes a file path or a DataFrame as input and produces a
figure or writes an image.  No analytical computation (rebinning,
normalisation, aggregation) should live here.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.io as pio
from matplotlib import pyplot as plt
from plotly.subplots import make_subplots

from sshicstuff.core import schemas
from sshicstuff.core.io import detect_delimiter, require_exists

logger = logging.getLogger(__name__)

pio.kaleido.scope.mathjax = None  # faster Kaleido export

# Chromosome colour palette (hex, 24 colours)
CHR_COLOUR_PALETTE: list[str] = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
    "#393b79", "#5254a3", "#6b6ecf", "#9c9ede",
]

EMPTY_FIGURE = go.Figure(
    layout=go.Layout(
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        annotations=[dict(
            x=0.5, y=0.5, text="No data available",
            showarrow=False, font=dict(size=28),
        )],
        hovermode="closest",
        plot_bgcolor="white",
        paper_bgcolor="white",
    )
)


# ---------------------------------------------------------------------------
# Colour helpers
# ---------------------------------------------------------------------------

def generate_colours(colour_type: str, n: int, alpha: float = 0.8, seed: int = 42) -> list[str]:
    """Return *n* colours in ``'hex'`` or ``'rgba'`` format."""
    import random
    rng = random.Random(seed)

    if colour_type == "hex":
        colours = [f"#{rng.randint(0, 0xFFFFFF):06x}" for _ in range(n)]
        fallback = "#0080FF"
        return [fallback if c == "#ffffff" else c for c in colours]

    if colour_type == "rgba":
        colours = [
            f"rgba({rng.randint(0, 255)}, {rng.randint(0, 255)}, {rng.randint(0, 255)}, {alpha})"
            for _ in range(n)
        ]
        fallback = f"rgba(0, 128, 255, {alpha})"
        white = f"rgba(255, 255, 255, {alpha})"
        return [fallback if c == white else c for c in colours]

    raise ValueError(f"colour_type must be 'hex' or 'rgba', got '{colour_type}'.")


# ---------------------------------------------------------------------------
# Chromosome colour bar
# ---------------------------------------------------------------------------

def build_chr_colourbar(df_bins: pd.DataFrame) -> tuple[go.Bar, list[float]]:
    """Build a Plotly chromosome-coloured bar trace for the genome view.

    Parameters
    ----------
    df_bins:
        Binned template DataFrame with ``chr`` and ``chr_bins`` columns.

    Returns
    -------
    (bar_trace, chr_tick_positions)
    """
    chr_list = df_bins[schemas.COL_CHR].values
    chr_bins = df_bins[schemas.COL_CHR_BINS].values
    chrom_to_idx = {c: i for i, c in enumerate(pd.unique(df_bins[schemas.COL_CHR]))}
    n_bins = len(chr_bins)

    full_bins, x_colours = [], []
    prev_idx = None
    chr_start = 0
    chr_tick_pos = []

    for ii, chrom in enumerate(chr_list):
        idx = chrom_to_idx[chrom]
        full_bins.append(f"{idx}:{chr_bins[ii]}")
        x_colours.append(CHR_COLOUR_PALETTE[idx % len(CHR_COLOUR_PALETTE)])
        if idx != prev_idx:
            if prev_idx is not None:
                chr_tick_pos.append((ii + chr_start) / 2)
            chr_start = ii
            prev_idx = idx

    chr_tick_pos.append((n_bins + chr_start) / 2)

    bar = go.Bar(
        x=full_bins,
        y=[1] * n_bins,
        marker=dict(color=x_colours, line=dict(width=0)),
        showlegend=False,
        hoverinfo="none",
    )
    return bar, chr_tick_pos


def _build_bins_template(df_coords: pd.DataFrame, bin_size: int) -> pd.DataFrame:
    """Build a complete bin template DataFrame for all chromosomes."""
    parts = []
    for _, row in df_coords.iterrows():
        n = row[schemas.COL_LENGTH] // bin_size + 1
        bins = np.arange(0, n * bin_size, bin_size)
        parts.append(pd.DataFrame({
            schemas.COL_CHR: row[schemas.COL_CHR],
            schemas.COL_CHR_BINS: bins,
            schemas.COL_GENOME_BINS: np.arange(0, len(bins) * bin_size, bin_size),
        }))
    return pd.concat(parts, ignore_index=True)


# ---------------------------------------------------------------------------
# 4C profile plots (static export)
# ---------------------------------------------------------------------------

def plot_profiles(
    profile_path: str | Path,
    oligo_capture_path: str | Path,
    chr_coord_path: str | Path,
    output_dir: str | Path | None = None,
    extension: str = "pdf",
    region: str | None = None,
    rolling_window: int = 1,
    log_scale: bool = False,
    user_y_min: float | None = None,
    user_y_max: float | None = None,
    width: int = 1200,
    height: int = 600,
) -> None:
    """Write one image per probe/group from a 4C-like profile.

    Parameters
    ----------
    profile_path:
        Binned or unbinned probe profile TSV.
    oligo_capture_path:
        Oligo capture table (to map fragment IDs → probe names).
    chr_coord_path:
        Chromosome coordinate file.
    output_dir:
        Destination directory.
    extension:
        Output image format (``"pdf"``, ``"png"``, ``"svg"``).
    region:
        Genomic region string ``"chrN:start-end"`` or ``"chrN:"``.
    rolling_window:
        Smoothing window size in bins.
    log_scale:
        Apply log10 scaling.
    user_y_min / user_y_max:
        Manual y-axis limits.
    width / height:
        Figure dimensions in pixels.
    """
    profile_path = Path(profile_path)
    oligo_capture_path = Path(oligo_capture_path)
    chr_coord_path = Path(chr_coord_path)

    require_exists(profile_path)
    require_exists(oligo_capture_path)
    require_exists(chr_coord_path)

    measure = "frequencies" if "frequencies" in profile_path.name else "contacts"

    df = pd.read_csv(str(profile_path), sep="\t")
    frag_cols = [c for c in df.columns if re.match(r"^\d+$|^\$", str(c))]

    df_oligo = pd.read_csv(str(oligo_capture_path), sep=detect_delimiter(oligo_capture_path))
    frag_to_probe = dict(
        zip(df_oligo[schemas.COL_FRAGMENT].astype(str), df_oligo[schemas.COL_NAME].astype(str))
    )

    df_coords = pd.read_csv(str(chr_coord_path), sep=detect_delimiter(chr_coord_path))
    df_coords.columns = [c.lower() for c in df_coords.columns]

    # Determine bin size from filename
    m = re.search(r"_(\d+)kb_profile_", profile_path.name)
    if m and m.group(1) != "0":
        binsize = int(m.group(1)) * 1000
        bin_suffix = f"{m.group(1)}kb"
    else:
        binsize = 0
        bin_suffix = "0kb"

    sample_name = profile_path.name.split(f"_{bin_suffix}_")[0]

    if output_dir is None:
        output_dir = profile_path.parent
    scale_tag = "log" if log_scale else "linear"
    out_dir = Path(output_dir) / "plots" / bin_suffix / scale_tag
    out_dir.mkdir(parents=True, exist_ok=True)

    # Region filtering
    x_min, x_max = 0, int(df_coords[schemas.COL_LENGTH].sum())
    chr_region = None
    if region:
        parts = region.split(":")
        chr_region = parts[0]
        coord_str = parts[1] if len(parts) > 1 else ""
        max_len = int(df_coords.loc[df_coords[schemas.COL_CHR] == chr_region, schemas.COL_LENGTH].values[0])
        if coord_str and "-" in coord_str:
            x_min, x_max = map(int, coord_str.split("-"))
        else:
            x_min, x_max = 0, max_len
        df = df[(df[schemas.COL_CHR] == chr_region) & (df.get("start", df.get(schemas.COL_CHR_BINS, 0)).between(x_min, x_max))]

    # Smoothing
    if rolling_window > 1 and binsize > 0:
        for chrom in df[schemas.COL_CHR].unique():
            mask = df[schemas.COL_CHR] == chrom
            df.loc[mask, frag_cols] = (
                df.loc[mask, frag_cols]
                .rolling(window=rolling_window, min_periods=1)
                .mean()
            )

    # Axis config
    y_min = float(user_y_min) if user_y_min is not None else 0.0
    y_max = float(user_y_max) if user_y_max is not None else float(df[frag_cols].max().max())
    x_col = "start" if binsize == 0 else schemas.COL_CHR_BINS
    if not chr_region:
        x_col = schemas.COL_GENOME_START if binsize == 0 else schemas.COL_GENOME_BINS

    log_tag = "log_" if log_scale else ""
    colours = generate_colours("rgba", len(frag_cols))

    if chr_region:
        for ii, frag in enumerate(frag_cols):
            probe = frag_to_probe.get(frag, "")
            y = df[frag].values.copy().astype(float)
            if log_scale:
                y[y == 0] = np.nan
                y = np.log10(y)

            fig = go.Figure(go.Scattergl(
                x=df[x_col], y=y, mode="lines",
                line=dict(width=1, color=colours[ii]),
                showlegend=False,
            ))
            fig.update_layout(
                title=f"{sample_name}",
                xaxis=dict(title=f"{chr_region} coordinates", tickformat="d",
                           range=[x_min, x_max], showgrid=False),
                yaxis=dict(title=f"{measure.capitalize()}{' (log)' if log_scale else ''}",
                           range=[y_min, y_max], showgrid=False),
                plot_bgcolor="white", paper_bgcolor="white",
                width=width, height=height,
            )
            out_path = out_dir / f"{sample_name}_{frag}_{probe}_{measure}_{bin_suffix}_{log_tag}{region}.{extension}"
            pio.write_image(fig, str(out_path), engine="kaleido")

    else:
        df_10kb = _build_bins_template(df_coords, 10_000)
        colourbar, chr_tick_pos = build_chr_colourbar(df_10kb)

        for ii, frag in enumerate(frag_cols):
            probe = frag_to_probe.get(frag, "")
            y = df[frag].values.copy().astype(float)
            if log_scale:
                y[y == 0] = np.nan
                y = np.log10(y)

            fig = make_subplots(
                rows=2, cols=1, row_heights=[0.94, 0.06], vertical_spacing=0.06,
                specs=[[{"type": "scatter"}], [{"type": "bar"}]],
            )
            fig.add_trace(
                go.Scatter(x=df[x_col], y=y, mode="lines",
                           line=dict(width=1, color=colours[ii]),
                           showlegend=False),
                row=1, col=1,
            )
            fig.add_trace(colourbar, row=2, col=1)

            chrom_labels = df[schemas.COL_CHR].unique().tolist()
            fig.update_layout(
                title=f"{sample_name} — {frag} {('— ' + probe) if probe else ''}",
                xaxis=dict(title=dict(text="Genomic coordinates", standoff=80),
                           tickformat="d", range=[x_min, x_max], showgrid=False),
                xaxis2=dict(tickmode="array", tickvals=chr_tick_pos,
                            ticktext=chrom_labels, tickfont=dict(size=12)),
                yaxis=dict(title=f"{measure.capitalize()}{' (log)' if log_scale else ''}",
                           range=[y_min, y_max], showgrid=False),
                yaxis2=dict(showticklabels=False),
                hovermode="closest",
                plot_bgcolor="white", paper_bgcolor="white",
                width=width, height=height,
            )
            out_path = out_dir / f"{sample_name}_{frag}_{probe}_{measure}_{bin_suffix}_{log_tag}.{extension}"
            pio.write_image(fig, str(out_path), engine="kaleido")

    logger.info("[Plot] Profiles written to %s", out_dir)


# ---------------------------------------------------------------------------
# Probe × probe heatmap (matplotlib)
# ---------------------------------------------------------------------------

def plot_probe_matrix(
    matrix_path: str | Path,
    output_path: str | Path | None = None,
    probes_rows: list[str] | None = None,
    probes_cols: list[str] | None = None,
    title: str | None = None,
    log_scale: bool = False,
    cmap: str = "Reds",
    vmin: float = 0,
    vmax: float | None = None,
    dpi: int = 300,
    fig_size: int = 14,
) -> None:
    """Render a probe × probe matrix as a matplotlib heatmap.

    Parameters
    ----------
    matrix_path:
        Probe matrix TSV (first column is the index).
    output_path:
        Destination image path.  Defaults to a PNG next to *matrix_path*.
    probes_rows / probes_cols:
        Subset the matrix rows / columns.
    title:
        Figure title.
    log_scale:
        Apply log10(x + 1) transformation.
    cmap:
        Matplotlib colourmap.
    vmin / vmax:
        Colour scale limits.
    dpi / fig_size:
        Figure resolution and size in inches.
    """
    matrix_path = Path(matrix_path)
    require_exists(matrix_path)

    if output_path is None:
        output_path = matrix_path.with_suffix(".png")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(str(matrix_path), sep="\t", index_col=0)
    if probes_rows:
        df = df.loc[probes_rows, :]
    if probes_cols:
        df = df.loc[:, probes_cols]

    arr = df.values.astype(float)
    if log_scale:
        arr = np.log10(arr + 1)

    if vmax is None:
        vmax = float(np.percentile(arr, 99))

    n_rows, n_cols = df.shape
    if n_rows == n_cols:
        fig_w = fig_h = fig_size
    elif n_rows > n_cols:
        fig_w, fig_h = fig_size, fig_size * n_rows / n_cols
    else:
        fig_w, fig_h = fig_size * n_cols / n_rows, fig_size

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    img = ax.imshow(arr, vmin=vmin, vmax=vmax, cmap=cmap, interpolation="none", aspect="auto")

    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label(
        "log10(contacts + 1)" if log_scale else "Contacts",
        fontsize=12,
    )
    cbar.ax.tick_params(labelsize=11)

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(df.columns, rotation=65, ha="right", fontsize=11)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(df.index, fontsize=11)
    ax.set_xlabel("Probes", fontsize=12)
    ax.set_ylabel("Probes", fontsize=12)

    if title:
        ax.set_title(title, fontsize=14)
    elif not title:
        ax.set_title(matrix_path.stem, fontsize=14)

    fig.savefig(str(output_path), bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    logger.info("[Plot] Probe matrix heatmap → %s", output_path.name)


# ---------------------------------------------------------------------------
# GUI interactive figure builder
# ---------------------------------------------------------------------------

def build_interactive_figure(
    df: pd.DataFrame,
    df_coords: pd.DataFrame,
    sample_name: str,
    probes: list[str],
    binsize: int,
    rolling_window: int,
    chr_region: str | None,
    log_scale: bool,
    user_x_min: str | None,
    user_x_max: str | None,
    user_y_min: str | None,
    user_y_max: str | None,
    width: int = 1200,
    height: int = 600,
) -> go.Figure:
    """Build a Plotly figure for the interactive profile viewer (GUI/Dash).

    Parameters
    ----------
    df:
        Profile DataFrame (binned or unbinned).
    df_coords:
        Chromosome coordinates.
    sample_name:
        Sample identifier for the figure title.
    probes:
        Column names in *df* to display.
    binsize:
        Active bin size (0 = fragment-level).
    rolling_window:
        Smoothing window (1 = no smoothing).
    chr_region:
        If set, restrict the view to this chromosome.
    log_scale / user_x_min / user_x_max / user_y_min / user_y_max:
        Display options.
    width / height:
        Figure dimensions in pixels.
    """
    n_chr = len(df_coords)
    full_genome_size = int(df_coords[schemas.COL_LENGTH].sum())
    x_min = float(user_x_min) if user_x_min else 0
    x_max = float(user_x_max) if user_x_max else full_genome_size
    x_col = schemas.COL_GENOME_START if binsize == 0 else schemas.COL_GENOME_BINS

    if chr_region:
        max_len = int(
            df_coords.loc[df_coords[schemas.COL_CHR] == chr_region, schemas.COL_LENGTH].values[0]
        )
        x_min = float(user_x_min) if user_x_min else 0
        x_max = float(user_x_max) if user_x_max else max_len
        df = df[df[schemas.COL_CHR] == chr_region]
        x_col = schemas.COL_START if binsize == 0 else schemas.COL_CHR_BINS

    if binsize > 0 and rolling_window > 1:
        for chrom in df[schemas.COL_CHR].unique():
            mask = df[schemas.COL_CHR] == chrom
            df.loc[mask, probes] = (
                df.loc[mask, probes]
                .rolling(window=rolling_window, min_periods=1)
                .mean()
            )

    data = df[probes].values.astype(float)
    if log_scale:
        data[data == 0] = np.nan
        data = np.log10(data)
        df = df.copy()
        df[probes] = data

    y_min = float(user_y_min) if user_y_min else float(np.nanmin(data)) if log_scale else 0.0
    y_max = float(user_y_max) if user_y_max else float(np.nanmax(data)) if log_scale else float(df[probes].max().max())
    y_ticks = np.linspace(y_min, y_max, 5)
    colours = generate_colours("rgba", len(probes))

    if chr_region:
        fig = go.Figure()
        for ii, probe in enumerate(probes):
            fig.add_trace(go.Scattergl(
                x=df[x_col], y=df[probe], name=probe, mode="lines",
                line=dict(width=1, color=colours[ii]),
            ))
        fig.update_layout(
            title=sample_name,
            xaxis=dict(title=f"{chr_region} coordinates", tickformat="d",
                       range=[x_min, x_max], showgrid=False),
            yaxis=dict(title="Contact frequency", tickvals=y_ticks,
                       range=[y_min, y_max], showgrid=False),
            plot_bgcolor="white", paper_bgcolor="white",
            width=width, height=height,
        )
    else:
        df_10kb = _build_bins_template(df_coords, 10_000)
        colourbar, chr_tick_pos = build_chr_colourbar(df_10kb)
        chrom_labels = df[schemas.COL_CHR].unique().tolist()

        fig = make_subplots(
            rows=2, cols=1, row_heights=[0.94, 0.06], vertical_spacing=0.06,
            specs=[[{"type": "scatter"}], [{"type": "bar"}]],
        )
        for ii, probe in enumerate(probes):
            fig.add_trace(
                go.Scattergl(x=df[x_col], y=df[probe], name=probe, mode="lines",
                             line=dict(width=1, color=colours[ii])),
                row=1, col=1,
            )
        fig.add_trace(colourbar, row=2, col=1)
        fig.update_layout(
            title=sample_name,
            xaxis=dict(title=dict(text="Genomic coordinates", standoff=80),
                       tickformat="d", range=[x_min, x_max], showgrid=False, fixedrange=True),
            xaxis2=dict(tickmode="array", tickvals=chr_tick_pos, ticktext=chrom_labels,
                        tickfont=dict(size=12), fixedrange=True),
            yaxis=dict(title="Contact frequency", tickvals=y_ticks, range=[y_min, y_max],
                       showgrid=False, fixedrange=True),
            yaxis2=dict(showticklabels=False, fixedrange=True),
            hovermode="closest",
            plot_bgcolor="white", paper_bgcolor="white",
            width=width, height=height,
        )

    return fig