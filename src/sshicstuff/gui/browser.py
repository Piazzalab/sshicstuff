import re
import os

import pandas as pd

# dash
from dash import callback
from dash.dependencies import Input, Output, State

# common.py
from sshicstuff.gui.common import TEMPORARY_DIRECTORY
from sshicstuff.gui.common import empty_figure
from sshicstuff.gui.common import uploaded_files
from sshicstuff.gui.common import save_file
from sshicstuff.gui.common import chr_to_exclude

# layout.py
from sshicstuff.gui.layout import layout
import sshicstuff.gui.graph as graph

if not os.path.exists(TEMPORARY_DIRECTORY):
    os.makedirs(TEMPORARY_DIRECTORY)
    

@callback(
    Output('slider-output-container', 'children'),
    [Input('binning-slider', 'value')])
def update_output(value):
    return f'Binning resolution : {value} kb'


@callback(
    [Output("oligo-dropdown", "options"),
     Output("coord-dropdown", "options"),
     Output("clear-list", "n_clicks"),
     Output("samples-dropdown", "options"),],
    [Input("upload-files", "filename"),
     Input("upload-files", "contents"),
     Input("clear-list", "n_clicks")],
)
def update_file_list(uploaded_filenames, uploaded_file_contents, n_clicks):
    if uploaded_filenames is not None and uploaded_file_contents is not None:
        for name, data in zip(uploaded_filenames, uploaded_file_contents):
            save_file(name, data)

    files = uploaded_files()
    if n_clicks is not None:
        if n_clicks > 0:
            for filename in files:
                os.remove(os.path.join(TEMPORARY_DIRECTORY, filename))
            files = []

    n_clicks = 0
    if len(files) == 0:
        return files, files, n_clicks, files

    else:
        inputs = []
        samples = []
        for f in files:
            if "profile" in f:
                samples.append({'label': f, 'value': os.path.join(TEMPORARY_DIRECTORY, f)})
            else:
                inputs.append({'label': f, 'value': os.path.join(TEMPORARY_DIRECTORY, f)})
        return inputs, inputs, n_clicks, samples


@callback(
    Output("region-dropdown", "options"),
    [Input("coord-dropdown", "value")],
)
def update_region_dropdown(coord_value):
    if coord_value is None:
        return []

    df = pd.read_csv(coord_value, sep='\t')
    df = df[~df['chr'].isin(chr_to_exclude)]
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
    df_coords = df_coords[~df_coords['chr'].isin(chr_to_exclude)]
    df_coords = df_coords[["chr", "length"]]
    df_coords["chr_start"] = df_coords["length"].shift().fillna(0).astype("int64")
    df_coords["cumu_start"] = df_coords["chr_start"].cumsum()

    # sample reading
    sample_name = samples_value.split('/')[-1].split('.')[0]
    df_samples = pd.read_csv(samples_value, sep='\t')
    df = df_samples[["chr", "start", "sizes", "genome_start"] + probes_value]
    df = df[~df['chr'].isin(chr_to_exclude)]

    binsize = 0
    if binning_value:
        binsize = binning_value * 1000

    figure = graph.figure_maker(
        binsize=binsize,
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

