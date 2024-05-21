import re

import numpy as np
from dash import html, dcc, callback
import dash_bootstrap_components as dbc
import dash_daq as daq
from dash.dependencies import Input, Output, State

from sshicstuff.gui.common import *


if not os.path.exists(TEMPORARY_DIRECTORY):
    os.makedirs(TEMPORARY_DIRECTORY)
    

layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H6('Upload files'),
            dcc.Upload(
                id="upload-files",
                children=html.Div(
                    ["Drag and drop or click to select a file to upload."]
                ),
                style={
                    "width": "100%",
                    "height": "80px",
                    "lineHeight": "80px",
                    "borderWidth": "2px",
                    "borderStyle": "dashed",
                    "borderRadius": "20px",
                    "textAlign": "center",
                },
                multiple=True,
            ),
        ], width=5, style={'margin-top': '0px', 'margin-bottom': '25px'}),

        dbc.Col([
            html.Button(
                id="clear-list",
                className="btn btn-danger",
                children="Clear list",
            )
        ], width=2, style={'margin-top': '50px', 'margin-bottom': '20px', 'margin-left': '20px'}),
    ]),

    dbc.Row([
        # Input files dropdown COLUMN
        dbc.Col([
            html.H6('Input files dropdown'),
            dbc.Row([
                dcc.Dropdown(
                    id='oligo-dropdown',
                    placeholder="Select capture oligos file",
                    multi=False),
            ], style={'margin-top': '10px', 'margin-bottom': '20px'}),

            dbc.Row([
                dcc.Dropdown(
                    id='coord-dropdown',
                    placeholder="Select chr. coordinates file",
                    multi=False),
            ], style={'margin-top': '10px', 'margin-bottom': '20px'}),

            dbc.Row([
                dcc.Dropdown(
                    id='samples-dropdown',
                    placeholder="Select sample file",
                    multi=False
                ),
            ], style={'margin-top': '10px', 'margin-bottom': '20px'}),

            dbc.Row([
                dcc.Dropdown(
                    id='probes-dropdown',
                    placeholder="Select probe(s) or group of probes",
                    multi=True
                ),
            ], style={'margin-top': '10px', 'margin-bottom': '20px'}),
        ], width=6),

        # Region & binning settings COLUMN
        dbc.Col([
            dbc.Row([
                html.Div(
                    id='slider-output-container', style={'font-size': '14px', 'margin-bottom': '10px'}),
                dcc.Slider(
                    id='binning-slider',
                    min=0,
                    max=100,
                    step=1,
                    value=10,
                    marks={i: str(i) for i in range(0, 101, 10)},
                    included=False,
                ),
            ]),

            dbc.Row([
                dbc.Label("Region", style={'font-size': '14px', 'margin-top': '20px'}),
            ]),

            dbc.Row([
                dbc.Col([
                    dcc.Dropdown(
                        id='region-dropdown',
                        value=None,
                        placeholder="Select chromosome",
                        multi=False),

                ], width=6),

                dbc.Col([
                    dcc.Input(
                        id='start-pos',
                        placeholder="Start",
                        value=None,
                        type='text',
                        className="custom-input"
                    ),
                ], width=3),

                dbc.Col([
                    dcc.Input(
                        id='end-pos',
                        placeholder="End",
                        value=None,
                        className="custom-input"
                    ),
                ], width=3),

            ]),
        ], width=6),
    ]),

    dbc.Row([
        dbc.Col([
            html.Button(
                id="plot-button", className="plot-button", children="Plot",
                style={'margin-top': '25px'}),
        ], width=2),

        dbc.Col([
            dbc.Label("Y min", style={'font-size': '14px', 'margin-top': '10px'}),
            dcc.Input(
                id='y-min', type='number', value=None,
                placeholder='Y min', className="custom-input"),
        ], width=1),

        dbc.Col([
            dbc.Label("Y max", style={'font-size': '14px', 'margin-top': '10px'}),
            dcc.Input(
                id='y-max', type='number', value=None,
                placeholder='Y max', className="custom-input"),
        ], width=1),

        dbc.Col([
            dbc.Label("Height", style={'font-size': '14px', 'margin-top': '10px'}),
            dcc.Input(
                id='height', type='number', value=600, step=20,
                placeholder='Height', className="custom-input"),
        ], width=1),

        dbc.Col([
            dbc.Label("Width", style={'font-size': '14px', 'margin-top': '10px'}),
            dcc.Input(
                id='width', type='number', value=1500, step=20,
                placeholder='Width', className="custom-input"),
        ], width=1),

        dbc.Col([
            daq.BooleanSwitch(
                id='re-scale-switch',
                on=False,
                label='Re-scale'
            ),
            html.Div(id='re-scale-output',
                     style={'margin-top': '10px', 'margin-left': '30px', 'font-size': '12px'}),

        ], width=1, style={'margin-top': '10px', 'margin-bottom': '10px', 'margin-left': '20px'}),

        dbc.Col([
            daq.BooleanSwitch(
                id='delimiter-switch',
                on=False,
                label='Delimiter'
            ),
        ], width=1, style={'margin-top': '10px', 'margin-bottom': '10px', 'margin-left': '20px'}),
    ]),

    dbc.Row([
        # html.Div(id='graphs', children=[], style={'margin-top': '20px', 'margin-bottom': '20px'}),
        dcc.Graph(
            id='graph',
            config={
                'displayModeBar': True,
                'scrollZoom': True,
                'doubleClick': 'reset',
                'autosizable': False
            },
            style={'height': '100%', 'width': '100%'},
            figure=empty_figure
        )
    ])
])


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
     State('delimiter-switch', 'on'),
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
        delimit,
        height,
        width
):

    re_scale_output = ""
    if n_clicks is None or n_clicks == 0:
        return empty_figure, re_scale_output
    if not samples_value or not probes_value:
        return empty_figure, re_scale_output

    fig = go.Figure()

    # coordinates & genomic (cumulative) positions stuff
    df_coords = pd.read_csv(coords_value, sep='\t')
    df_chr_len = df_coords[["chr", "length"]]
    df_chr_len["chr_start"] = df_chr_len["length"].shift().fillna(0).astype("int64")
    df_chr_len["cumu_start"] = df_chr_len["chr_start"].cumsum()

    df_samples = pd.read_csv(samples_value, sep='\t')
    sample_name = samples_value.split('/')[-1].split('.')[0]

    df = df_samples[["chr", "start", "sizes", "genome_start"] + probes_value]

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
        data = df[probes_value].values
        if y_max <= 1.:
            # sqrt transformation
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

        df[probes_value] = new_data

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
            )
        )

    fig.update_layout(
        width=width,
        height=height,
        title=f"{sample_name}",
        xaxis=dict(domain=[0.0, 0.9], title=x_label),
        yaxis=dict(title="Contact frequency"),
        xaxis_type='linear',
        xaxis_tickformat="d",
        xaxis_range=[x_min, x_max],
        yaxis_range=[y_min, y_max],
        hovermode='closest',
        plot_bgcolor='white',
        paper_bgcolor='white'
    )

    if delimit and not region_value:
        for xi, x_pos in enumerate(df_chr_len["cumu_start"].to_list()):
            name_pos = x_pos + 100
            fig.add_shape(
                type='line',
                yref='paper',
                xref='x',
                x0=x_pos, x1=x_pos,
                y0=0, y1=1,
                line=dict(color='gray', width=1, dash='dot')
            )

            fig.add_annotation(
                go.layout.Annotation(
                    x=name_pos,
                    y=1.07,
                    yref="paper",
                    text=df_chr_len.loc[xi, "chr"],
                    showarrow=False,
                    xanchor="center",
                    font=dict(size=11, color=colors_hex[xi]),
                    textangle=330
                ),
                xref="x"
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
