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

from dash import callback, dcc
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

# Chromosomes excluded from region selection (artificial / donor)
CHR_EXCLUDED = list(schemas.ARTIFICIAL_CHROMOSOMES) + ["chr_artificial_donor"]

SESSION_KEY = "sshicstuff_session_id"


# ---------------------------------------------------------------------------
# Per-session cache helpers
# ---------------------------------------------------------------------------

def get_user_cache_dir() -> str:
    """Return a per-session cache directory namespaced by the app instance.

    Structure: ``BASE_CACHE_DIR / APP_INSTANCE_ID / SESSION_ID``
    """
    sid = session.get(SESSION_KEY)
    if sid is None:
        sid = str(uuid4())
        session[SESSION_KEY] = sid

    cache_dir = os.path.join(str(BASE_CACHE_DIR), APP_INSTANCE_ID, sid)
    os.makedirs(cache_dir, exist_ok=True)
    return cache_dir


# ---------------------------------------------------------------------------
# Callbacks
# ---------------------------------------------------------------------------

@callback(
    Output("binning-slider-output-container", "children"),
    Input("binning-slider", "value"),
)
def update_binning_label(value):
    return f"Binning resolution: {value} kb"


@callback(
    Output("window-slider-output-container", "children"),
    Input("window-slider", "value"),
)
def update_smoothing_label(value):
    return f"Smoothing window: {value}"


@callback(
    Output("alert-upload-4c", "children"),
    Input("upload-files-4c", "filename"),
    prevent_initial_call=True,
)
def show_upload_alert(filenames):
    if filenames:
        return dbc.Alert(
            f"Uploaded {len(filenames)} file(s): {', '.join(filenames)}",
            color="success",
            dismissable=True,
        )
    return None


@callback(
    [
        Output("oligo-dropdown", "options"),
        Output("coord-dropdown", "options"),
        Output("clear-list-4c", "n_clicks"),
        Output("alert-clean-cache-4c", "children"),
        Output("samples-dropdown", "options"),
    ],
    [
        Input("upload-files-4c", "filename"),
        Input("upload-files-4c", "contents"),
        Input("clear-list-4c", "n_clicks"),
    ],
)
def update_file_list(uploaded_filenames, uploaded_file_contents, n_clicks):
    cache_dir = get_user_cache_dir()

    if uploaded_filenames is not None and uploaded_file_contents is not None:
        for name, data in zip(uploaded_filenames, uploaded_file_contents):
            save_file_cache(name, data, cache_dir)

    files = uploaded_files_cache(cache_dir)
    clear_alert = None

    if n_clicks and n_clicks > 0:
        for filename in files:
            os.remove(os.path.join(cache_dir, filename))
        files = []
        clear_alert = dbc.Alert("Cache cleared.", color="info", dismissable=True)

    if not files:
        return [], [], 0, clear_alert, []

    inputs, samples = [], []
    for f in files:
        full_path = os.path.join(cache_dir, f)
        if "profile" in f:
            samples.append({"label": f, "value": full_path})
        else:
            inputs.append({"label": f, "value": full_path})

    return inputs, inputs, 0, clear_alert, samples


@callback(
    Output("region-dropdown", "options"),
    Input("coord-dropdown", "value"),
)
def update_region_dropdown(coord_value):
    if coord_value is None:
        return []

    df = pd.read_csv(coord_value, sep="\t")
    df = df[~df["chr"].isin(CHR_EXCLUDED)]
    chrs = df["chr"].unique().tolist()
    return [{"label": c, "value": c} for c in chrs]


@callback(
    Output("probes-dropdown", "options"),
    [
        Input("oligo-dropdown", "value"),
        Input("samples-dropdown", "value"),
    ],
)
def update_probes_dropdown(oligo_value, sample_value):
    if sample_value is None:
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


@callback(
    Output("graph", "figure"),
    Input("plot-button", "n_clicks"),
    [
        State("binning-slider", "value"),
        State("window-slider", "value"),
        State("coord-dropdown", "value"),
        State("samples-dropdown", "value"),
        State("probes-dropdown", "value"),
        State("region-dropdown", "value"),
        State("start-pos", "value"),
        State("end-pos", "value"),
        State("y-min", "value"),
        State("y-max", "value"),
        State("log-scale-switch", "on"),
        State("height", "value"),
        State("width", "value"),
    ],
)
def update_graph(
    n_clicks,
    binning_value,
    rolling_value,
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

    # Chromosome coordinates (cumulative positions)
    df_coords = pd.read_csv(coords_value, sep="\t")
    df_coords = df_coords[~df_coords["chr"].isin(CHR_EXCLUDED)]
    df_coords = df_coords[["chr", "length"]].copy()
    df_coords["chr_start"] = df_coords["length"].shift().fillna(0).astype("int64")
    df_coords["cumu_start"] = df_coords["chr_start"].cumsum()

    sample_name = os.path.basename(samples_value).split(".")[0]
    df_full = pd.read_csv(samples_value, sep="\t")

    coord_cols = [
        schemas.COL_CHR, schemas.COL_START, schemas.COL_SIZES, schemas.COL_GENOME_START
    ]
    df = df_full[[c for c in coord_cols if c in df_full.columns] + probes_value]
    df = df[~df[schemas.COL_CHR].isin(CHR_EXCLUDED)]

    binsize = (binning_value or 0) * 1000

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
        width=width or 1600,
        height=height or 600,
    )


@callback(
    Output("download-figure-pdf", "data"),
    Input("btn-figure-pdf", "n_clicks"),
    State("graph", "figure"),
)
def export_figure_pdf(n_clicks, figure):
    if not n_clicks:
        return None
    tmp_pdf = "/tmp/sshicstuff_figure_export.pdf"
    pio.write_image(figure, tmp_pdf, format="pdf")
    return dcc.send_file(tmp_pdf)