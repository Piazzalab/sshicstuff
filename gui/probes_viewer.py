import pandas as pd
import dash
import json
from os import listdir
from dash import html
from dash import dcc
from dash import callback
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State, ALL, MATCH
import plotly.graph_objs as go

from gui.common import *
from binning import rebin_live

if not os.path.exists(TEMPORARY_DIRECTORY):
    os.makedirs(TEMPORARY_DIRECTORY)


layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H2('Probes Viewer'),
        ], width=12, style={'margin-top': '20px', 'margin-bottom': '20px'})
    ]),

    dbc.Row([
        dbc.Col([
            html.H6('Upload files'),
            dcc.Upload(
                id="pv-upload-files",
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
                id="pv-clear-list",
                className="btn btn-danger",
                children="Clear list",
            )
        ], width=2, style={'margin-top': '50px', 'margin-bottom': '20px', 'margin-left': '20px'}),
    ]),

    dbc.Row([
        dbc.Col([
            dcc.Dropdown(id='pv-oligo-dropdown',
                         placeholder="Select capture oligos file",
                         multi=False),
        ], width=6, style={'margin-top': '10px', 'margin-bottom': '20px'}),

        dbc.Col([
            dcc.Dropdown(id='pv-coord-dropdown',
                         placeholder="Select chr. coordinates file",
                         multi=False),
        ], width=6, style={'margin-top': '10px', 'margin-bottom': '20px'}),
    ]),

    dbc.Row([
        dbc.Col([
            dcc.Input(id='pv-number-probes', type='number', value=2, step='1',
                      placeholder='How many cards :',
                      style={
                          'width': '100%',
                          'border': '1px solid #ccc',
                          'border-radius': '4px',
                          'padding': '5px',
                          'font-size': '16px',
                          'background-color': '#fff',
                          'color': '#333'
                      }),
        ], width=2, style={'margin-top': '20px'}),

        dbc.Col([
            html.Div(id='pv-slider-output-container',
                     style={'margin-top': '10px', 'font-size': '16px', 'margin-bottom': '10px'}),
            dcc.Slider(
                id='pv-binning-slider',
                min=0,
                max=100,
                step=1,
                value=10,
                marks={i: str(i) for i in range(0, 101, 10)},
                included=False,
            ),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '10px', 'margin-left': '0px'}),

        dbc.Col([
            dcc.Checklist(
                id="pv-sync-box",
                options=[{"label": "Synchronize axis", "value": "sync"}],
                value=[],
                inline=True,
                className='custom-checkbox-label',
                labelStyle={"margin": "5px"}
            )
        ], width=2, style={'margin-top': '20px', 'margin-bottom': '10px', 'margin-left': '0px'}),

        dbc.Col([
            html.Button(id="pv-plot-button", className="plot-button", children="Plot"),
        ], width=2, style={'margin-top': '20px', 'margin-bottom': '0px', 'margin-left': '0px'}),

        dcc.Store(id='pv-stored-graphs-axis-range', data={}),
        html.Div(id='pv-dynamic-probes-cards', children=[], style={'margin-top': '20px', 'margin-bottom': '20px'}),
        html.Div(id='pv-graphs', children=[], style={'margin-top': '20px', 'margin-bottom': '20px'}),
    ]),
])


@callback(
    [Output("pv-oligo-dropdown", "options"),
     Output("pv-coord-dropdown", "options"),
     Output("pv-clear-list", "n_clicks")],
    [Input("pv-upload-files", "filename"),
     Input("pv-upload-files", "contents"),
     Input("pv-clear-list", "n_clicks")],
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
        return files,  n_clicks
    else:
        options = []
        for f in files:
            if "profile" in f:
                continue
            options.append({'label': f, 'value': os.path.join(TEMPORARY_DIRECTORY, f)})
        return options, options, n_clicks


def update_table(file_path, delim):
    if file_path and delim:
        df = pd.read_csv(file_path, sep=delim)
        data = df.to_dict('records')
        columns = [{"name": i, "id": i} for i in df.columns]
        return data, columns
    return None, None


@callback(
    Output('pv-slider-output-container', 'children'),
    [Input('pv-binning-slider', 'value')])
def update_output(value):
    return f'Binning resolution : {value} kb'


def create_card(index, sample_options, probe_options, graph_options, sample_value, probe_value, graph_value):

    card = dbc.Col(
        dbc.Card([
            dbc.CardHeader(html.Div(
                id={'type': 'pv-probe-card-header', 'index': index},
                style={'font-size': '12px'})),

            dbc.CardBody([
                dbc.Row([
                    dbc.Col([
                        dcc.Dropdown(
                            options=sample_options,
                            value=sample_value,
                            placeholder="Select sample",
                            id={'type': 'sample-dropdown', 'index': index},
                            multi=False,
                        )
                    ], style={'margin-top': '10px', 'margin-bottom': '20px'}),
                ]),
                dbc.Row([
                    dbc.Col([
                        dcc.Dropdown(
                            options=probe_options,
                            value=probe_value,
                            placeholder="Select probe",
                            id={'type': 'probe-dropdown', 'index': index},
                            multi=False,
                        )
                    ], width=8, style={'margin-top': '10px', 'margin-bottom': '20px'}),

                    dbc.Col([
                        dcc.Dropdown(
                            id={'type': 'graph-dropdown', 'index': index},
                            options=graph_options,
                            value=graph_value,
                            placeholder="Select a graph",
                            multi=False
                        )
                    ], width=4, style={'margin-top': '10px'}),
                ])
            ])
        ])
    )

    return card


@callback(
    Output('pv-dynamic-probes-cards', 'children'),
    Input('pv-number-probes', 'value'),
    Input('pv-oligo-dropdown', 'value'),
    State('pv-dynamic-probes-cards', 'children'),
)
def update_probes_cards(n_cards, capture_oligos, cards_children):
    if n_cards is None or n_cards == 0:
        return []

    all_samples_items = sorted([f for f in listdir(TEMPORARY_DIRECTORY) if "profile" in f])
    samples_options = [{'label': s, 'value': s} for s in all_samples_items]
    graph_options = [{'label': f'graph {x}', 'value': f'graph {x}'} for x in range(n_cards)]

    probes_options = []
    if capture_oligos:
        df = pd.read_csv(capture_oligos)
        probes = df['name'].to_list()
        fragments = df['fragment'].to_list()
        probes_options = [{'label': p, 'value': f} for p, f in zip(probes, fragments)]

    existing_cards = []
    displaying_cards = []
    if cards_children:
        existing_cards = cards_children[0]['props']['children']
        for ii, item in enumerate(existing_cards):
            cardbody = item['props']['children']['props']['children'][1]['props']['children']
            displaying_cards.append(
                create_card(
                    index=ii,
                    sample_options=samples_options,
                    probe_options=probes_options,
                    graph_options=graph_options,
                    sample_value=cardbody[0]['props']['children'][0]['props']['children'][0]['props']['value'],
                    probe_value=cardbody[1]['props']['children'][0]['props']['children'][0]['props']['value'],
                    graph_value=cardbody[1]['props']['children'][1]['props']['children'][0]['props']['value']
                ))

    if len(existing_cards) > n_cards:
        displaying_cards = existing_cards[:n_cards]
    if len(existing_cards) < n_cards:
        cards_to_add = n_cards - len(existing_cards)
        for i in range(cards_to_add):
            displaying_cards.append(create_card(
                index=len(existing_cards) + i,
                sample_options=samples_options,
                probe_options=probes_options,
                graph_options=graph_options,
                sample_value=None,
                probe_value=None,
                graph_value=None
            ))

    rows = []
    for i in range(0, len(displaying_cards), 2):
        row = dbc.Row(displaying_cards[i:i + 2], style={'margin-top': '20px', 'margin-bottom': '20px'})
        rows.append(row)
    return rows


@callback(
    Output({'type': 'pv-probe-card-header', 'index': MATCH}, 'children'),
    Input({'type': 'probe-dropdown', 'index': MATCH}, 'value'),
    Input({'type': 'sample-dropdown', 'index': MATCH}, 'value'),
)
def update_card_header(probe_value, sample_value):
    if sample_value is None or probe_value is None:
        return None

    samp_id = sample_value.split('.')[0]
    return f"{samp_id} - {probe_value}"


def update_figure(
        graph_id: int, graph_dict: dict, traces_colors: list, binning: int,
        df_coords: pd.DataFrame, x_range=None, y_range=None):

    fig = go.Figure()
    trace_id = 0

    df_chr_len = df_coords[["chr", "length"]]
    df_chr_len["chr_start"] = df_chr_len["length"].shift().fillna(0).astype("int64")
    df_chr_len["cumu_start"] = df_chr_len["chr_start"].cumsum()

    for j in range(graph_dict['size']):
        samp = graph_dict['samples'][j]
        frag = graph_dict['fragments'][j]
        filepath = graph_dict['filepaths'][j]
        df = pd.read_csv(filepath, sep='\t')[["chr", "start", "sizes", "genome_start", frag]]

        if binning > 0:
            x_col = "genome_bins"
            binning *= 1000  # kbp convert to bp
            df = rebin_live(df, binning, df_coords)
        else:
            x_col = "genome_start"

        fig.add_trace(
            go.Scattergl(
                x=df[x_col],
                y=df[frag],
                name=f"{samp} - {frag}",
                mode='lines+markers',
                line=dict(width=1, color=traces_colors[trace_id]),
                marker=dict(size=4)
            )
        )

        fig.update_layout(
            width=1500,
            height=500,
            title=f"Graphe {graph_id}",
            xaxis=dict(domain=[0.0, 0.9], title="Genomic position (bp)"),
            yaxis=dict(title="Contact frequency"),
            hovermode='closest'
        )
        trace_id += 1

    if x_range:
        fig.update_xaxes(range=x_range)
    if y_range:
        fig.update_yaxes(range=y_range)

    for xi, x_pos in enumerate(df_chr_len.cumu_start.to_list()):
        name_pos = x_pos + 100
        fig.add_shape(type='line',
                      yref='paper',
                      xref='x',
                      x0=x_pos, x1=x_pos,
                      y0=0, y1=1,
                      line=dict(color='gray', width=1, dash='dot'))

        fig.add_annotation(
            go.layout.Annotation(
                x=name_pos,
                y=1.07,
                yref="paper",
                text=df_chr_len.chr[xi],
                showarrow=False,
                xanchor="center",
                font=dict(size=11, color=colors_hex[xi]),
                textangle=330
            ),
            xref="x"
        )
    return fig


@callback(
    Output('pv-graphs', 'children'),
    Input('pv-plot-button', 'n_clicks'),
    Input('pv-stored-graphs-axis-range', 'data'),
    State('pv-binning-slider', 'value'),
    State('pv-coord-dropdown', 'value'),
    State({'type': 'sample-dropdown', 'index': ALL}, 'value'),
    State({'type': 'probe-dropdown', 'index': ALL}, 'value'),
    State({'type': 'graph-dropdown', 'index': ALL}, 'value')
)
def update_graphs( n_clicks, axis_range, binning_value, coords_value, samples_value, probes_value, graphs_values):
    ctx = dash.callback_context
    triggerd_input = ctx.triggered[0]['prop_id'].split('.')[0]

    if n_clicks is None or n_clicks == 0:
        return None

    graphs_info = {}
    nb_graphs = 0

    x_range = None
    y_range = None
    if triggerd_input == 'pv-stored-graphs-axis-range':
        if axis_range:
            x_range = axis_range['x_range']
            y_range = axis_range['y_range']

    for i, graph in enumerate(graphs_values):
        if graph is None:
            continue
        graph_id = int(graph.split(' ')[-1])
        if graph_id not in graphs_info:
            nb_graphs += 1
            graphs_info[graph_id] = {
                'samples': [],
                'fragments': [],
                'filepaths': [],
                'size': 0,
            }
        graphs_info[graph_id]['samples'].append(samples_value[i])
        graphs_info[graph_id]['fragments'].append(str(probes_value[i]))
        graphs_info[graph_id]['filepaths'].append(join(TEMPORARY_DIRECTORY, samples_value[i]))
        graphs_info[graph_id]['size'] += 1

    chr_coords_path = str(join(TEMPORARY_DIRECTORY, coords_value))
    df_coords = pd.read_csv(chr_coords_path, sep='\t')

    figures = {}
    traces_count = 0
    for i in graphs_info:
        traces_to_add = graphs_info[i]['size']
        figures[i] = update_figure(
            graph_id=i,
            graph_dict=graphs_info[i],
            traces_colors=colors_rgba[traces_count:traces_count + traces_to_add],
            binning=binning_value,
            df_coords=df_coords,
            x_range=x_range,
            y_range=y_range
        )
        traces_count += traces_to_add

    graphs_layout = []
    for i in sorted(figures.keys()):
        graphs_layout.append(
            dcc.Graph(id={'type': 'graph', 'index': i},
                      config={'displayModeBar': True, 'scrollZoom': True},
                      style={'height': 'auto', 'width': '100%'},
                      figure=figures[i])
        )
    return graphs_layout


@callback(
    Output('pv-stored-graphs-axis-range', 'data'),
    Input({'type': 'graph', 'index': ALL}, 'relayoutData'),
    State('pv-sync-box', 'value')
)
def update_figure_range(relayout_data, sync_value):
    if not sync_value:
        return dash.no_update

    if len(relayout_data) == 1:
        return dash.no_update

    if not any(relayout_data):
        return [None, None]

    ctx = dash.callback_context
    input_id = json.loads(ctx.triggered[0]['prop_id'].split('.')[0])['index']

    updated_range = {'x_range': None, 'y_range': None}
    if 'xaxis.range[0]' in relayout_data[input_id]:
        updated_range['x_range'] = [relayout_data[input_id]['xaxis.range[0]'],
                                    relayout_data[input_id]['xaxis.range[1]']]

    if 'yaxis.range[0]' in relayout_data[input_id]:
        updated_range['y_range'] = [relayout_data[input_id]['yaxis.range[0]'],
                                    relayout_data[input_id]['yaxis.range[1]']]
    return updated_range
