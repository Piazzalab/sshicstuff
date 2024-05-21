import re
import os
import numpy as np
import pandas as pd

# plotly
import plotly.graph_objs as go
from plotly.subplots import make_subplots

# dash
from dash import callback
from dash.dependencies import Input, Output, State

# common.py
from sshicstuff.gui.common import TEMPORARY_DIRECTORY
from sshicstuff.gui.common import colors_rgba
from sshicstuff.gui.common import empty_figure
from sshicstuff.gui.common import uploaded_files
from sshicstuff.gui.common import save_file
from sshicstuff.gui.common import rebin_live
from sshicstuff.gui.common import transform_data
from sshicstuff.gui.common import make_colorbar
from sshicstuff.gui.common import chr_to_exclude

# layout.py
from sshicstuff.gui.layout import layout

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
    Output('re-scale-output', 'children'),
    Input('plot-button', 'n_clicks'),
    Input('graph', 'relayoutData'),
    [State('binning-slider', 'value'),
     State('coord-dropdown', 'value'),
     State('samples-dropdown', 'value'),
     State('probes-dropdown', 'value'),
     State('region-dropdown', 'value'),
     State('start-pos', 'value'),
     State('end-pos', 'value'),
     State('y-min', 'value'),
     State('y-max', 'value'),
     State('re-scale-switch', 'on'),
     State('height', 'value'),
     State('width', 'value'),]
)
def update_graph(
        n_clicks,
        relayout_data,
        binning_value,
        coords_value,
        samples_value,
        probes_value,
        region_value,
        user_x_min,
        user_x_max,
        user_y_min,
        user_y_max,
        re_scale,
        height,
        width
):
    re_scale_output = ""

    if n_clicks is None or n_clicks == 0:
        return empty_figure, re_scale_output
    if not samples_value or not probes_value:
        return empty_figure, re_scale_output

    # coordinates & genomic (cumulative) positions stuff
    df_coords = pd.read_csv(coords_value, sep='\t')
    df_coords = df_coords[~df_coords['chr'].isin(chr_to_exclude)]
    df_chr_len = df_coords[["chr", "length"]]
    df_chr_len["chr_start"] = df_chr_len["length"].shift().fillna(0).astype("int64")
    df_chr_len["cumu_start"] = df_chr_len["chr_start"].cumsum()

    df_samples = pd.read_csv(samples_value, sep='\t')
    sample_name = samples_value.split('/')[-1].split('.')[0]

    df = df_samples[["chr", "start", "sizes", "genome_start"] + probes_value]
    df = df[~df['chr'].isin(chr_to_exclude)]

    if binning_value > 0:
        binning_value *= 1000  # kbp convert to bp
        df = rebin_live(df, binning_value, df_coords)

    if region_value:
        df = df[df["chr"] == region_value]
        x_max_basal = df_chr_len.loc[df_chr_len["chr"] == region_value]["length"].tolist()[0]
        x_min = int(user_x_min) if user_x_min else 0
        x_max = int(user_x_max) if user_x_max else x_max_basal
        if x_max > x_max_basal:
            x_max = x_max_basal
        x_label = f"{region_value} position (bp)"
        x_col = "chr_bins" if binning_value else "start"

    else:
        x_min = 0
        x_max = df_chr_len["cumu_start"].max()
        x_label = "Genomic position (bp)"
        x_col = "genome_bins" if binning_value else "genome_start"

    if binning_value:
        x_min = x_min // binning_value * binning_value
        x_max = x_max // binning_value * binning_value
        df = df[(df["chr_bins"] >= x_min) & (df["chr_bins"] <= x_max)]
    else:
        df = df[(df["start"] >= x_min) & (df["start"] <= x_max)]

    y_min = float(user_y_min) if user_y_min else 0
    y_max = float(user_y_max) if user_y_max else df[probes_value].max().max()

    if re_scale:
        data = df[probes_value]
        new_data, y_max, y_min, re_scale_output = transform_data(data, y_max, user_y_max, y_min, re_scale)
        df[probes_value] = new_data

    fig = make_subplots(
        rows=2, cols=1, row_heights=[0.94, 0.06], vertical_spacing=0.06,
        specs=[[{'type': 'scatter'}], [{'type': 'bar'}]]
    )

    for j in range(len(probes_value)):
        frag = probes_value[j]
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

    color_bar, chr_tick_pos = make_colorbar(df_coords, binning_value)
    unique_chr_list = df_coords['chr'].unique()

    fig.add_trace(color_bar, row=2, col=1)
    fig.update_layout(
        title=f"{sample_name}",
        xaxis=dict(
            # title=x_label
        ),
        xaxis2=dict(
            tickmode='array',
            tickvals=chr_tick_pos,
            ticktext=unique_chr_list,
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

    if relayout_data:
        if 'xaxis.range[0]' in relayout_data:
            new_x_min = relayout_data['xaxis.range[0]']
            new_x_max = relayout_data['xaxis.range[1]']
            if user_y_max:
                y_max = user_y_max
            else:
                y_max = df[(df[x_col] >= new_x_min) & (df[x_col] <= new_x_max)][probes_value].max().max()
            fig.update_layout(xaxis_range=[new_x_min, new_x_max], yaxis_range=[y_min, y_max])

    return fig, re_scale_output
