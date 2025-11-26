"""
Browser callbacks
"""
import os
import re
from uuid import uuid4

import pandas as pd
import plotly.io as pio
import dash_bootstrap_components as dbc

from dash import callback, dcc
from dash.dependencies import Input, Output, State
from flask import session

# Interpret __CACHE_DIR__ as *base* directory for all sessions
from sshicstuff.core.methods import (
    uploaded_files_cache,
    save_file_cache,
    __CACHE_DIR__ as BASE_CACHE_DIR,
    APP_INSTANCE_ID,
)
from sshicstuff.core.plot import empty_figure, figure_maker

CHR_ARTIFICIAL_EXCLUSION = ["chr_artificial_donor", "chr_artificial_ssDNA"]

SESSION_KEY = "sshicstuff_session_id"

def get_user_cache_dir() -> str:
    """
    Return a per-session cache directory, namespaced by the app instance.

    Directory structure:
        BASE_CACHE_DIR / APP_INSTANCE_ID / SESSION_ID
    """
    sid = session.get(SESSION_KEY)
    if sid is None:
        sid = str(uuid4())
        session[SESSION_KEY] = sid

    cache_dir = os.path.join(str(BASE_CACHE_DIR), APP_INSTANCE_ID, sid)
    os.makedirs(cache_dir, exist_ok=True)
    return cache_dir

@callback(
    Output('binning-slider-output-container', 'children'),
    [Input('binning-slider', 'value')])
def update_binning_output(value):
    return f'Binning resolution : {value} kb'


@callback(
    Output('window-slider-output-container', 'children'),
    [Input('window-slider', 'value')])
def update_smoothing_output(value):
    return f'Smoothing window : {value}'


@callback(
    Output("alert-upload-4c", "children"),
    Input("upload-files-4c", "filename"),
    prevent_initial_call=True
)
def show_upload_alert(filenames):
    if filenames:
        return dbc.Alert(f"Uploaded {len(filenames)} files : {', '.join(filenames)}", color="success", dismissable=True)
    return None



@callback(
    [Output("oligo-dropdown", "options"),
     Output("coord-dropdown", "options"),
     Output("clear-list-4c", "n_clicks"),
     Output("alert-clean-cache-4c", "children"),
     Output("samples-dropdown", "options"),],
    [Input("upload-files-4c", "filename"),
     Input("upload-files-4c", "contents"),
     Input("clear-list-4c", "n_clicks")],
)
def update_file_list(uploaded_filenames, uploaded_file_contents, n_clicks):
    cache_dir = get_user_cache_dir()

    if uploaded_filenames is not None and uploaded_file_contents is not None:
        for name, data in zip(uploaded_filenames, uploaded_file_contents):
            save_file_cache(name, data, cache_dir)

    files = uploaded_files_cache(cache_dir)
    clear_alert = None
    if n_clicks is not None:
        if n_clicks > 0:
            for filename in files:
                os.remove(os.path.join(cache_dir, filename))
            files = []
            clear_alert = dbc.Alert("Cache cleared", color="info", dismissable=True)

    n_clicks = 0
    if len(files) == 0:
        return files, files, n_clicks, clear_alert, files

    else:
        inputs = []
        samples = []
        for f in files:
            full_path = os.path.join(cache_dir, f)
            if "profile" in f:
                samples.append({'label': f, 'value': full_path})
            else:
                inputs.append({'label': f, 'value': full_path})
        return inputs, inputs, n_clicks, clear_alert, samples


@callback(
    Output("region-dropdown", "options"),
    [Input("coord-dropdown", "value")],
)
def update_region_dropdown(coord_value):
    if coord_value is None:
        return []

    df = pd.read_csv(coord_value, sep='\t')
    df = df[~df['chr'].isin(CHR_ARTIFICIAL_EXCLUSION)]
    chr_list = df['chr'].unique()
    chr_list = [f"{c}" for c in chr_list]

    return [{'label': c, 'value': c} for c in chr_list]


@callback(
    Output("probes-dropdown", "options"),
    [Input("oligo-dropdown", "value"),
     Input("samples-dropdown", "value")],
)
def update_probes_dropdown(oligo_value, sample_value):
    if sample_value is None:
        return []

    df = pd.read_csv(sample_value, sep='\t')
    col_of_interest = [c for c in df.columns if re.match(r'^\d+$|^\$', c)]

    probes_options = []
    probes_to_frag = {}

    if oligo_value:
        df2 = pd.read_csv(oligo_value)
        probes_to_frag = dict(zip(df2['fragment'].astype(str), df2['name'].astype(str)))

    for c in col_of_interest:
        label = f"{c} - {probes_to_frag[c]}" if c in probes_to_frag else c
        probes_options.append({'label': label, 'value': c})

    return probes_options


@callback(
    Output('graph', 'figure'),
    Input('plot-button', 'n_clicks'),
    [State('binning-slider', 'value'),
     State('window-slider', 'value'),
     State('coord-dropdown', 'value'),
     State('samples-dropdown', 'value'),
     State('probes-dropdown', 'value'),
     State('region-dropdown', 'value'),
     State('start-pos', 'value'),
     State('end-pos', 'value'),
     State('y-min', 'value'),
     State('y-max', 'value'),
     State('log-scale-switch', 'on'),
     State('height', 'value'),
     State('width', 'value'),]
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
        width
):
    
    # If the button has not been clicked, return an empty figure
    # adn if the samples or probes are not selected, return an empty figure as well
    if n_clicks is None or n_clicks == 0:
        return empty_figure
    if not samples_value or not probes_value:
        return empty_figure

    # coordinates & genomic (cumulative) positions stuff
    df_coords = pd.read_csv(coords_value, sep='\t')
    df_coords = df_coords[~df_coords['chr'].isin(CHR_ARTIFICIAL_EXCLUSION)]
    df_coords = df_coords[["chr", "length"]]
    df_coords["chr_start"] = df_coords["length"].shift().fillna(0).astype("int64")
    df_coords["cumu_start"] = df_coords["chr_start"].cumsum()

    # sample reading
    sample_name = samples_value.split('/')[-1].split('.')[0]
    df_samples = pd.read_csv(samples_value, sep='\t')
    df = df_samples[["chr", "start", "sizes", "genome_start"] + probes_value]
    df = df[~df['chr'].isin(CHR_ARTIFICIAL_EXCLUSION)]

    binsize = 0
    if binning_value:
        binsize = binning_value * 1000

    figure = figure_maker(
        binsize=binsize,
        rolling_window=rolling_value,
        df_coords=df_coords,
        df=df,
        sample_name=sample_name,
        probes=probes_value,
        chr_region=region_value,
        log_scale=log_scale,
        user_x_min=user_x_min,
        user_x_max=user_x_max,
        user_y_min=user_y_min,
        user_y_max=user_y_max,
        width=width,
        height=height,
    )

    return figure


@callback(
    Output("download-figure-pdf", "data"),
    Input("btn-figure-pdf", "n_clicks"),
    [State("graph", "figure")],
)
def export_figure(n_clicks, figure):
    if n_clicks is None or n_clicks == 0:
        return

    tmp_pdf = "figure.pdf"
    pio.write_image(figure, tmp_pdf, format="pdf")
    return dcc.send_file(tmp_pdf)

