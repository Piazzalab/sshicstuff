import re
from dash import html, dcc, callback
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State, ALL
import plotly.graph_objs as go

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
        dbc.Col([
            dcc.Dropdown(id='oligo-dropdown',
                         placeholder="Select capture oligos file",
                         multi=False),
        ], width=6, style={'margin-top': '10px', 'margin-bottom': '20px'}),

        dbc.Col([
            dcc.Dropdown(id='coord-dropdown',
                         placeholder="Select chr. coordinates file",
                         multi=False),
        ], width=6, style={'margin-top': '10px', 'margin-bottom': '20px'}),
    ]),

    dbc.Row([
        dbc.Col([
            dcc.Dropdown(
                id='samples-dropdown',
                placeholder="Select sample file",
                multi=False
            ),
        ], width=6, style={'margin-top': '10px', 'margin-bottom': '20px'}),

        dbc.Col([
            dcc.Dropdown(
                id='probes-dropdown',
                placeholder="Select probe(s) or group of probes",
                multi=True
            ),
        ], width=6, style={'margin-top': '10px', 'margin-bottom': '20px'}),
    ]),

    dbc.Row([
        dbc.Col([
            html.Div(id='slider-output-container',
                     style={'margin-top': '10px', 'font-size': '16px', 'margin-bottom': '10px'}),
            dcc.Slider(
                id='binning-slider',
                min=0,
                max=100,
                step=1,
                value=10,
                marks={i: str(i) for i in range(0, 101, 10)},
                included=False,
            ),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '10px', 'margin-left': '0px'}),

        dbc.Col([
            html.Button(id="plot-button", className="plot-button", children="Plot"),
        ], width=2, style={'margin-top': '20px', 'margin-bottom': '0px', 'margin-left': '0px'}),
    ]),

    html.Div(id='graphs', children=[], style={'margin-top': '20px', 'margin-bottom': '20px'}),
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
    Output('graphs', 'children'),
    Input('plot-button', 'n_clicks'),
    [State('binning-slider', 'value'),
     State('coord-dropdown', 'value'),
     State('samples-dropdown', 'value'),
     State('probes-dropdown', 'value')]
)
def update_graph(n_clicks, binning_value, coords_value, samples_value, probes_value):
    if n_clicks is None or n_clicks == 0:
        return None
    if not samples_value or not probes_value:
        return None

    fig = go.Figure()
    df_coords = pd.read_csv(coords_value, sep='\t')
    df_samples = pd.read_csv(samples_value, sep='\t')
    sample_name = samples_value.split('/')[-1].split('.')[0]

    df = df_samples[["chr", "start", "sizes", "genome_start"] + probes_value]
    if binning_value > 0:
        x_col = "genome_bins"
        binning_value *= 1000  # kbp convert to bp
        df = rebin_live(df, binning_value, df_coords)
    else:
        x_col = "genome_start"

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
        width=1500,
        height=600,
        title=f"{sample_name}",
        xaxis=dict(domain=[0.0, 0.9], title="Genomic position (bp)"),
        yaxis=dict(title="Contact frequency"),
        hovermode='closest',
        plot_bgcolor='white',
        paper_bgcolor='white'
    )

    df_chr_len = df_coords[["chr", "length"]]
    df_chr_len["chr_start"] = df_chr_len["length"].shift().fillna(0).astype("int64")
    df_chr_len["cumu_start"] = df_chr_len["chr_start"].cumsum()

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

    graph_layout = dcc.Graph(
        config={'displayModeBar': True, 'scrollZoom': True},
        style={'height': 'auto', 'width': '100%'},
        figure=fig
    )
    return graph_layout

