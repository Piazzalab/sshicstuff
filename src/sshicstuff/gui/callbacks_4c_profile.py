"""
Dash callbacks for the ssHiC Browser (4C-like profile viewer) tab.
"""

from __future__ import annotations

import os
import re
from uuid import uuid4

import pandas as pd
import plotly.io as pio
import dash_bootstrap_components as dbc

from dash import callback, dcc, html, no_update
from dash.dependencies import Input, Output, State
from flask import session

from sshicstuff.gui.cache import (
    __CACHE_DIR__ as BASE_CACHE_DIR,
    APP_INSTANCE_ID,
    save_file_cache,
    uploaded_files_cache,
)
from sshicstuff.plot import EMPTY_FIGURE, build_interactive_figure
from sshicstuff.core import schemas

CHR_EXCLUDED = list(schemas.ARTIFICIAL_CHROMOSOMES) + ["chr_artificial_donor"]
SESSION_KEY = "sshicstuff_session_id"


# ---------------------------------------------------------------------------
# Session cache
# ---------------------------------------------------------------------------

def get_user_cache_dir() -> str:
    sid = session.get(SESSION_KEY)
    if sid is None:
        sid = str(uuid4())
        session[SESSION_KEY] = sid
    cache_dir = os.path.join(str(BASE_CACHE_DIR), APP_INSTANCE_ID, sid)
    os.makedirs(cache_dir, exist_ok=True)
    return cache_dir


# ---------------------------------------------------------------------------
# Slider labels
# ---------------------------------------------------------------------------

@callback(
    Output("binning-slider-output-container", "children"),
    Input("binning-slider", "value"),
)
def update_binning_label(value):
    return [
        "Binning resolution: ",
        html.Span(f"{value} kb", style={"fontWeight": "700", "color": "#2245b7"}),
    ]


@callback(
    Output("window-slider-output-container", "children"),
    Input("window-slider", "value"),
)
def update_smoothing_label(value):
    return [
        "Smoothing window: ",
        html.Span(str(value), style={"fontWeight": "700", "color": "#2245b7"}),
    ]


# ---------------------------------------------------------------------------
# Advanced settings collapse toggle
# ---------------------------------------------------------------------------

@callback(
    Output("advanced-collapse", "is_open"),
    Input("advanced-toggle", "n_clicks"),
    State("advanced-collapse", "is_open"),
    prevent_initial_call=True,
)
def toggle_advanced(n_clicks, is_open):
    return not is_open


# ---------------------------------------------------------------------------
# File management
# ---------------------------------------------------------------------------

@callback(
    [
        Output("oligo-dropdown",    "options"),
        Output("coord-dropdown",    "options"),
        Output("samples-dropdown",  "options"),
        Output("clear-list-4c",     "n_clicks"),
        Output("alert-clean-cache-4c", "children"),
        Output("alert-upload-4c",   "children"),
    ],
    [
        Input("upload-files-4c",    "filename"),
        Input("upload-files-4c",    "contents"),
        Input("clear-list-4c",      "n_clicks"),
    ],
)
def update_file_list(uploaded_filenames, uploaded_file_contents, n_clicks):
    cache_dir = get_user_cache_dir()
    upload_alert = None

    if uploaded_filenames and uploaded_file_contents:
        for name, data in zip(uploaded_filenames, uploaded_file_contents):
            save_file_cache(name, data, cache_dir)
        upload_alert = dbc.Alert(
            f"Uploaded {len(uploaded_filenames)} file(s): {', '.join(uploaded_filenames)}",
            color="success", dismissable=True, className="py-2",
        )

    files = uploaded_files_cache(cache_dir)
    clear_alert = None

    if n_clicks and n_clicks > 0:
        for f in files:
            os.remove(os.path.join(cache_dir, f))
        files = []
        clear_alert = dbc.Alert("Cache cleared.", color="info",
                                dismissable=True, className="py-2")

    if not files:
        return [], [], [], 0, clear_alert, upload_alert

    inputs, samples = [], []
    for f in files:
        full_path = os.path.join(cache_dir, f)
        if "profile" in f:
            samples.append({"label": f, "value": full_path})
        else:
            inputs.append({"label": f, "value": full_path})

    return inputs, inputs, samples, 0, clear_alert, upload_alert


# ---------------------------------------------------------------------------
# Loaded-file badges
# ---------------------------------------------------------------------------

def _loaded_badge(value):
    if not value:
        return None
    name = os.path.basename(value)
    return html.Span(f"✓ {name}", className="file-loaded-badge")


@callback(Output("oligo-loaded-badge",  "children"), Input("oligo-dropdown",   "value"))
def badge_oligo(v):  return _loaded_badge(v)

@callback(Output("coord-loaded-badge",  "children"), Input("coord-dropdown",   "value"))
def badge_coord(v):  return _loaded_badge(v)

@callback(Output("sample-loaded-badge", "children"), Input("samples-dropdown", "value"))
def badge_sample(v): return _loaded_badge(v)


# ---------------------------------------------------------------------------
# Chromosome region dropdown
# ---------------------------------------------------------------------------

@callback(
    Output("region-dropdown", "options"),
    Input("coord-dropdown", "value"),
)
def update_region_dropdown(coord_value):
    if not coord_value:
        return []
    df = pd.read_csv(coord_value, sep="\t")
    df = df[~df["chr"].isin(CHR_EXCLUDED)]
    chrs = df["chr"].unique().tolist()
    return [{"label": c, "value": c} for c in chrs]


@callback(
    Output("chr-length-info", "children"),
    [Input("region-dropdown", "value"), Input("coord-dropdown", "value")],
)
def update_chr_length(chr_region, coord_value):
    if not chr_region or not coord_value:
        return ""
    df = pd.read_csv(coord_value, sep="\t")
    row = df[df["chr"] == chr_region]
    if row.empty:
        return ""
    length = int(row["length"].values[0])
    return f"Length: {length:,} bp"


# ---------------------------------------------------------------------------
# Probes dropdown
# ---------------------------------------------------------------------------

@callback(
    Output("probes-dropdown", "options"),
    [Input("oligo-dropdown", "value"), Input("samples-dropdown", "value")],
)
def update_probes_dropdown(oligo_value, sample_value):
    if not sample_value:
        return []
    df = pd.read_csv(sample_value, sep="\t")
    frag_cols = [c for c in df.columns if re.match(r"^\d+$|^\$", str(c))]

    frag_to_probe: dict[str, str] = {}
    if oligo_value:
        df_oligo = pd.read_csv(oligo_value)
        frag_to_probe = dict(
            zip(df_oligo["fragment"].astype(str), df_oligo["name"].astype(str))
        )

    return [
        {"label": f"{c} — {frag_to_probe[c]}" if c in frag_to_probe else c, "value": c}
        for c in frag_cols
    ]


# ---------------------------------------------------------------------------
# Main plot
# ---------------------------------------------------------------------------

@callback(
    Output("graph", "figure"),
    Input("plot-button", "n_clicks"),
    [
        State("binning-slider",      "value"),
        State("window-slider",       "value"),
        State("normalization-radio", "value"),
        State("coord-dropdown",      "value"),
        State("samples-dropdown",    "value"),
        State("probes-dropdown",     "value"),
        State("region-dropdown",     "value"),
        State("start-pos",           "value"),
        State("end-pos",             "value"),
        State("y-min",               "value"),
        State("y-max",               "value"),
        State("log-scale-switch",    "on"),
        State("height",              "value"),
        State("width",               "value"),
    ],
)
def update_graph(
    n_clicks,
    binning_value,
    rolling_value,
    normalization,
    coords_value,
    samples_value,
    probes_value,
    region_value,
    user_x_min,
    user_x_max,
    user_y_min,
    user_y_max,
    log_scale,
    height,
    width,
):
    if not n_clicks or not samples_value or not probes_value:
        return EMPTY_FIGURE

    df_coords = pd.read_csv(coords_value, sep="\t")
    df_coords = df_coords[~df_coords["chr"].isin(CHR_EXCLUDED)]
    df_coords = df_coords[["chr", "length"]].copy()

    sample_name = os.path.basename(samples_value).split(".")[0]
    df_full = pd.read_csv(samples_value, sep="\t")

    coord_cols = [
        schemas.COL_CHR, schemas.COL_START, schemas.COL_SIZES, schemas.COL_GENOME_START
    ]
    df = df_full[[c for c in coord_cols if c in df_full.columns] + probes_value]
    df = df[~df[schemas.COL_CHR].isin(CHR_EXCLUDED)]

    binsize = (binning_value or 0) * 1000

    # Apply normalisation before plotting
    df = _apply_normalization(df, probes_value, normalization or "raw")

    return build_interactive_figure(
        df=df,
        df_coords=df_coords,
        sample_name=sample_name,
        probes=probes_value,
        binsize=binsize,
        rolling_window=rolling_value or 1,
        chr_region=region_value,
        log_scale=bool(log_scale),
        user_x_min=str(user_x_min) if user_x_min is not None else None,
        user_x_max=str(user_x_max) if user_x_max is not None else None,
        user_y_min=str(user_y_min) if user_y_min is not None else None,
        user_y_max=str(user_y_max) if user_y_max is not None else None,
        width=width or 1400,
        height=height or 600,
    )


def _apply_normalization(df: pd.DataFrame, probes: list[str], mode: str) -> pd.DataFrame:
    """Normalise probe columns in-place (returns a copy)."""
    if mode == "raw":
        return df
    df = df.copy()
    if mode == "fraction_viewpoint":
        for probe in probes:
            total = df[probe].sum()
            if total > 0:
                df[probe] = df[probe] / total
    elif mode == "fraction_global":
        grand_total = df[probes].values.sum()
        if grand_total > 0:
            df[probes] = df[probes] / grand_total
    return df


# ---------------------------------------------------------------------------
# Exports
# ---------------------------------------------------------------------------

@callback(
    Output("download-figure-pdf", "data"),
    Input("btn-figure-pdf", "n_clicks"),
    State("graph", "figure"),
    prevent_initial_call=True,
)
def export_pdf(n_clicks, figure):
    if not n_clicks:
        return no_update
    tmp = "/tmp/sshicstuff_profile.pdf"
    pio.write_image(figure, tmp, format="pdf")
    return dcc.send_file(tmp)


@callback(
    Output("download-figure-svg", "data"),
    Input("btn-figure-svg", "n_clicks"),
    State("graph", "figure"),
    prevent_initial_call=True,
)
def export_svg(n_clicks, figure):
    if not n_clicks:
        return no_update
    tmp = "/tmp/sshicstuff_profile.svg"
    pio.write_image(figure, tmp, format="svg")
    return dcc.send_file(tmp)