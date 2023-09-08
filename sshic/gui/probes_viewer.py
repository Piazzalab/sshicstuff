import os
import pandas as pd
import numpy as np
import dash
import json
from os.path import join
from os import listdir
from dash import html
from dash import dcc
from dash import callback
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State, ALL, MATCH
import plotly.graph_objs as go

from common import generate_data_table, prepare_dataframe_for_output
import core.utils

colors = [
    'rgba(0, 0, 255, 0.2)',  # blue
    'rgba(255, 0, 0, 0.2)',  # red
    'rgba(249, 172, 37, 0.2)',  # yellow
    'rgba(245, 0, 87, 0.2)',  # pink
    'rgba(29, 233, 182, 0.2)',  # green
    'rgba(255, 234, 0, 0.2)',  # yellow 2
    'rgba(255, 11, 0, 0.2)',  # orange
    'rgba(141, 110, 99, 0.2)',  # brown
    'rgba(255, 64, 129, 0.2)',  # pink 2
    'rgba(120, 144, 156, 0.2)',  # blue grey
    'rgba(0, 131, 143, 0.2)',  # cyan
    'rgba(171, 71, 188, 0.2)',  # purple
    'rgba(255, 152, 0, 0.2)',  # amber
    'rgba(0, 150, 136, 0.2)',  # teal
    'rgba(0, 184, 212, 0.2)',  # cyan 2
    'rgba(0, 200, 83, 0.2)',  # green 2
    'rgba(229, 115, 115, 0.2)',  # red 2
    'rgba(255, 167, 38, 0.2)',  # orange 2
    'rgba(61, 90, 254, 0.2)',  # indigo
    'rgba(68, 138, 255, 0.2)',  # blue 2
    'rgba(121, 134, 203, 0.2)',  # deep purple
    'rgba(170, 102, 68, 0.2)',  # deep orange
    'rgba(255, 171, 145, 0.2)',  # pink 3
    'rgba(255, 209, 128, 0.2)'  # amber 2
]

chr_names = [f"chr{i}" for i in range(1, 17)] + ["2_micron", "mitochondrion", "chr_artificial"]
chr_pos = [230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751,
           666816, 1078177, 924431, 784333, 1091291, 948066, 6318, 85779, 7828]
chr_colors = ['#000000', '#0c090a', '#2c3e50', '#34495e', '#7f8c8d', '#8e44ad', '#2ecc71', '#2980b9',
              '#f1c40f', '#d35400', '#e74c3c', '#c0392b', '#1abc9c', '#16a085', '#bdc3c7', '#2c3e50',
              '#7f8c8d', '#f39c12', '#27ae60']

layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H2('Probes Viewer'),
        ], width=12, style={'margin-top': '20px', 'margin-bottom': '20px'})
    ]),

    dbc.Row([
        dbc.Col([
            html.Label("Select number of cards to display:")
        ], width=4, style={'margin-top': '10px', 'margin-bottom': '0px'}),

        dbc.Col([
            html.Label("Select capture oligos table:")
        ], width=4, style={'margin-top': '10px', 'margin-bottom': '0px'}),

        dbc.Col([
            html.Label("Additional groups of probes (if any) :"),
        ], width=4, style={'margin-top': '10px', 'margin-bottom': '0px'})
    ]),

    dbc.Row([
        dbc.Col([
            dcc.Input(id='pv-number-probes', type='number', value=2, step='1',
                      placeholder='How many probes to compare :',
                      style={
                          'width': '100%',
                          'border': '1px solid #ccc',
                          'border-radius': '4px',
                          'padding': '6px',
                          'font-size': '16px',
                          'background-color': '#fff',
                          'color': '#333'
                      }),

            dbc.Tooltip("Specify the number of probes you want to compare",
                        target="pv-number-probes", className="custom-tooltip", placement="right"),
        ], width=4, style={'margin-top': '0px', 'margin-bottom': '20px'}),

        dbc.Col([
            dcc.Dropdown(id='pv-oligo-selector', multi=False),
        ], width=4, style={'margin-top': '0px', 'margin-bottom': '20px'}),

        dbc.Col([
            dcc.Dropdown(id='pv-probe-groups-selector', options=[], value=None, multi=False),
        ], width=4, style={'margin-top': '0px', 'margin-bottom': '30px'}),
    ]),

    dbc.Row([
        dbc.Col([
            html.Div(id='pv-p2f-dataframe-title', style={'margin-top': '20px', 'margin-bottom': '20px'}),
            dcc.Loading(generate_data_table('pv-p2f-dataframe', [], [], 5))
        ], width=4, style={'margin-top': '0px', 'margin-bottom': '20px'}),

        dbc.Col([
            html.Div(id='pv-groups-dataframe-title', style={'margin-top': '20px', 'margin-bottom': '20px'}),
            dcc.Loading(generate_data_table('pv-groups-dataframe', [], [], 5))
        ], width=7, style={'margin-top': '0px', 'margin-bottom': '20px', 'margin-left': '30px'}),
    ]),

    dbc.Row([
        dbc.Col([
            html.Label("Select binning :", style={'margin-top': '0px', 'margin-bottom': '10px'}),
            dcc.Slider(
                id='pv-binning-slider',
                min=0,
                max=100,
                step=1,
                value=10,
                marks={i: str(i) for i in range(0, 101, 10)},
                included=False,
            ),
            html.Div(id='pv-slider-output-container',
                     style={'margin-top': '10px', 'font-size': '14px'}),
            html.Br(),
        ], width=4, style={'margin-top': '0px', 'margin-bottom': '10px', 'margin-left': '20px'}),

        dbc.Col([
            dcc.Checklist(
                id="pv-sync-box",
                options=[{"label": "Sync axis", "value": "sync"}],
                value=[],
                inline=True,
                className='custom-checkbox-label',
                labelStyle={"margin": "5px"}
            )
        ], width=2, style={'margin-top': '15px', 'margin-bottom': '10px', 'margin-left': '20px'}),

        dbc.Col([
            html.Button(id="pv-plot-buttom", className="plot-button", children="Submit"),
        ], width=2, style={'margin-top': '20px', 'margin-bottom': '0px', 'margin-left': '20px'}),
    ]),

    dcc.Store(id='pv-stored-graphs-axis-range', data={}),
    html.Div(id='pv-dynamic-probes-cards', children=[], style={'margin-top': '20px', 'margin-bottom': '20px'}),
    html.Div(id='pv-graphs', children=[], style={'margin-top': '20px', 'margin-bottom': '20px'}),
])


@callback(
    Output('pv-oligo-selector', 'options'),
    Output('pv-probe-groups-selector', 'options'),
    Input('data-basedir', 'data')
)
def update_oligo_selector(data_basedir):
    if data_basedir is None:
        return [], []
    inputs_dir = join(data_basedir, 'inputs')
    inputs_files = sorted([f for f in listdir(inputs_dir) if os.path.isfile(join(inputs_dir, f))],
                          key=lambda x: x.lower())

    options = [{'label': f, 'value': join(inputs_dir, f)} for f in inputs_files]
    return options, options


@callback(
     Output('pv-p2f-dataframe-title', 'children'),
     Output('pv-p2f-dataframe', 'data'),
     Output('pv-p2f-dataframe', 'columns'),
     Input('pv-oligo-selector', 'value')
)
def oligo_and_fragments(oligo_file):

    if oligo_file is None:
        return None, [], []

    df_oli = pd.read_csv(oligo_file, sep=core.utils.detect_delimiter(oligo_file))
    title = html.H6("Oligo probes VS. Fragments ID:")
    if 'fragment' not in df_oli.columns:
        return "Select a capture oligos file first", [], []

    data, columns = prepare_dataframe_for_output(df_oli, ["name", "fragment"])
    return title, data, columns


@callback(
    [Output('pv-groups-dataframe-title', 'children'),
     Output('pv-groups-dataframe', 'data'),
     Output('pv-groups-dataframe', 'columns')],
    [Input('pv-probe-groups-selector', 'value')]
)
def display_df_probe_groups(groups_file):
    if groups_file is None:
        return None, [], []

    df_groups = pd.read_csv(groups_file, sep='\t')
    data, columns = prepare_dataframe_for_output(df_groups)
    title = html.H6("Probe groups :")
    return title, data, columns


@callback(
    Output('pv-slider-output-container', 'children'),
    [Input('pv-binning-slider', 'value')])
def update_output(value):
    return f'You have selected a binning of {value} kb'


def create_card(
        index,
        sample_options,
        probe_options,
        pcr_options,
        weight_options,
        graph_options,
        sample_value,
        probe_value,
        pcr_value,
        weight_value,
        graph_value
):

    card = dbc.Col(
        dbc.Card([
            dbc.CardHeader(html.Div(id={'type': 'pv-probe-card-header', 'index': index})),
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
                    ]),

                    dbc.Col([
                        dcc.Dropdown(
                            options=probe_options,
                            value=probe_value,
                            placeholder="Select probe",
                            id={'type': 'probe-dropdown', 'index': index},
                            multi=False,
                        )
                    ]),
                ]),

                dbc.Row([
                    dbc.Col([
                        dcc.Checklist(
                            options=pcr_options,
                            value=pcr_value,
                            id={'type': 'pcr-checkboxes', 'index': index},
                            inline=True,
                            className='custom-checkbox-label',
                            labelStyle={"margin": "5px"}
                        )
                    ]),

                    dbc.Col([
                        dcc.Checklist(
                            options=weight_options,
                            value=weight_value,
                            id={'type': 'weight-checkboxes', 'index': index},
                            inline=True,
                            className='custom-checkbox-label',
                            labelStyle={"margin": "5px"}
                        )
                    ]),
                ]),

                dbc.Row([
                    dbc.Col([
                        dcc.Dropdown(
                            id={'type': 'graph-dropdown', 'index': index},
                            options=graph_options,
                            value=graph_value,
                            placeholder="Select a graph",
                            multi=False
                        )
                    ])
                ])
            ])
        ])
    )

    return card


@callback(
    Output('pv-dynamic-probes-cards', 'children'),
    Input('pv-number-probes', 'value'),
    State('data-basedir', 'data'),
    State('pv-dynamic-probes-cards', 'children')
)
def update_probes_cards(n_cards, data_basedir, cards_children):
    if n_cards is None or n_cards == 0:
        return []

    pp_outputs_dir = join(data_basedir, 'outputs')
    all_samples_items = sorted(listdir(pp_outputs_dir))
    samples_options = [{'label': s, 'value': s} for s in all_samples_items if os.path.isdir(join(pp_outputs_dir, s))]
    graph_options = [{'label': f'graph {x}', 'value': f'graph {x}'} for x in range(n_cards)]

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
                    probe_options=cardbody[0]['props']['children'][1]['props']['children'][0]['props']['options'],
                    pcr_options=cardbody[1]['props']['children'][0]['props']['children'][0]['props']['options'],
                    weight_options=cardbody[1]['props']['children'][1]['props']['children'][0]['props']['options'],
                    graph_options=graph_options,
                    sample_value=cardbody[0]['props']['children'][0]['props']['children'][0]['props']['value'],
                    probe_value=cardbody[0]['props']['children'][1]['props']['children'][0]['props']['value'],
                    pcr_value=cardbody[1]['props']['children'][0]['props']['children'][0]['props']['value'],
                    weight_value=cardbody[1]['props']['children'][1]['props']['children'][0]['props']['value'],
                    graph_value=cardbody[2]['props']['children'][0]['props']['children'][0]['props']['value']
                ))

    if len(existing_cards) > n_cards:
        displaying_cards = existing_cards[:n_cards]
    if len(existing_cards) < n_cards:
        cards_to_add = n_cards - len(existing_cards)
        for i in range(cards_to_add):
            displaying_cards.append(create_card(
                index=len(existing_cards) + i,
                sample_options=samples_options,
                probe_options=[],
                pcr_options=[],
                weight_options=[],
                graph_options=graph_options,
                sample_value=None,
                probe_value=None,
                pcr_value=None,
                weight_value=None,
                graph_value=None
            ))

    rows = []
    for i in range(0, len(displaying_cards), 3):
        row = dbc.Row(displaying_cards[i:i+3], style={'margin-top': '20px', 'margin-bottom': '20px'})
        rows.append(row)
    return rows


@callback(
    Output({'type': 'pcr-checkboxes', 'index': MATCH}, 'options'),
    Input({'type': 'sample-dropdown', 'index': MATCH}, 'value'),
    State('data-basedir', 'data')
)
def update_pcr_checkboxes_options(sample_value, data_basedir):
    pp_outputs_dir = join(data_basedir, 'outputs')
    if sample_value is None:
        return []
    all_items = listdir(join(pp_outputs_dir, sample_value))
    pcr_dirs = [item for item in all_items if os.path.isdir(join(pp_outputs_dir, sample_value, item))]
    return [{'label': d, 'value': d} for d in pcr_dirs if 'pcr' in d.lower()]


@callback(
    Output({'type': 'pcr-checkboxes', 'index': MATCH}, 'value'),
    Input({'type': 'pcr-checkboxes', 'index': MATCH}, 'value'),
)
def update_pcr_checkboxes(pcr_value):
    if not pcr_value:
        return []
    else:
        return [pcr_value[-1]]


@callback(
    Output({'type': 'weight-checkboxes', 'index': MATCH}, 'options'),
    Input({'type': 'pcr-checkboxes', 'index': MATCH}, 'value'),
    Input({'type': 'sample-dropdown', 'index': MATCH}, 'value'),
    State('data-basedir', 'data')
)
def update_weight_checkboxes_options(pcr_value, sample_value, data_basedir):
    ctx = dash.callback_context
    triggerd_input = ctx.triggered[0]['prop_id'].split('.')[0]
    if triggerd_input == '':
        return []

    pp_outputs_dir = join(data_basedir, 'outputs')
    if not sample_value or not pcr_value or pcr_value == []:
        return []

    pcr_dir = join(pp_outputs_dir, sample_value, pcr_value[-1])
    weight_dirs = [w for w in listdir(pcr_dir) if os.path.isdir(join(pcr_dir, w))]
    return [{'label': d, 'value': d} for d in weight_dirs]


@callback(
    Output({'type': 'weight-checkboxes', 'index': MATCH}, 'value'),
    Input({'type': 'weight-checkboxes', 'index': MATCH}, 'value')
)
def update_weight_checkboxes(weight_value):
    if not weight_value:
        return []
    else:
        return [weight_value[-1]]


@callback(
    Output({'type': 'probe-dropdown', 'index': MATCH}, 'options'),
    Input({'type': 'weight-checkboxes', 'index': MATCH}, 'value'),
    Input({'type': 'sample-dropdown', 'index': MATCH}, 'value'),
    Input({'type': 'pcr-checkboxes', 'index': MATCH}, 'value'),
    State('data-basedir', 'data')
)
def update_probe_dropdown_options(weight_value, sample_value, pcr_value, data_basedir):
    ctx = dash.callback_context
    triggerd_input = ctx.triggered[0]['prop_id'].split('.')[0]
    if triggerd_input == '':
        return []

    pp_outputs_dir = join(data_basedir, 'outputs')
    if sample_value is None:
        return []
    if pcr_value is None or pcr_value == [] or weight_value is None or weight_value == []:
        return []

    items_dir = join(pp_outputs_dir, sample_value, pcr_value[-1], weight_value[-1])
    df = pd.read_csv(join(items_dir, f"{sample_value}_unbinned_contacts.tsv"), sep='\t')
    probes = [c for c in df.columns if c not in ['chr', 'start', 'sizes', 'genome_start', 'end']]
    return [{'label': f, 'value': f} for f in probes]


@callback(
    Output({'type': 'pv-probe-card-header', 'index': MATCH}, 'children'),
    Input({'type': 'probe-dropdown', 'index': MATCH}, 'value'),
    Input({'type': 'sample-dropdown', 'index': MATCH}, 'value'),
    Input({'type': 'pcr-checkboxes', 'index': MATCH}, 'value'),
    Input({'type': 'weight-checkboxes', 'index': MATCH}, 'value')
)
def update_card_header(probe_value, sample_value, pcr_value, weight_value):
    if sample_value is None or probe_value is None:
        return None
    if pcr_value is None or pcr_value == [] or weight_value is None or weight_value == []:
        return None

    samp_id = sample_value.split('/')[-1]
    return f"{samp_id} - {probe_value} - {pcr_value[-1]} - {weight_value[-1]}"


def update_figure(
    graph_id: int,
    graph_dict: dict,
    traces_colors: list,
    binning: int,
    chr_boundaries: list,
    x_range=None,
    y_range=None
):
    fig = go.Figure()
    trace_id = 0
    for j in range(graph_dict['size']):
        samp = graph_dict['samples'][j]
        frag = graph_dict['fragments'][j]
        pcr = graph_dict['pcr'][j]
        weight = graph_dict['weight'][j]
        filepath = graph_dict['filepaths'][j]
        df = pd.read_csv(filepath, sep='\t')

        x_col = "genome_bins" if binning > 0 else "genome_start"
        fig.add_trace(
            go.Scattergl(
                x=df[x_col],
                y=df[frag],
                name=f"{samp} - {frag} - {pcr} - {weight}",
                mode='lines+markers',
                line=dict(width=1, color=traces_colors[trace_id]),
                marker=dict(size=4)
            )
        )

        fig.update_layout(
            width=1500,
            height=500,
            title=f"Graphe {graph_id}",
            xaxis=dict(domain=[0.0, 0.9], title="Genome bins"),
            yaxis=dict(title="Contact frequency"),
            hovermode='closest'
        )
        trace_id += 1

    if x_range:
        fig.update_xaxes(range=x_range)
    if y_range:
        fig.update_yaxes(range=y_range)

    for xi, x_pos in enumerate(chr_boundaries):
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
                text=chr_names[xi],
                showarrow=False,
                xanchor="center",
                font=dict(size=11, color=chr_colors[xi]),
                textangle=330
            ),
            xref="x"
        )
    return fig


@callback(
    Output('pv-graphs', 'children'),
    Input('pv-plot-buttom', 'n_clicks'),
    Input('pv-stored-graphs-axis-range', 'data'),
    State('pv-binning-slider', 'value'),
    State({'type': 'sample-dropdown', 'index': ALL}, 'value'),
    State({'type': 'pcr-checkboxes', 'index': ALL}, 'value'),
    State({'type': 'weight-checkboxes', 'index': ALL}, 'value'),
    State({'type': 'probe-dropdown', 'index': ALL}, 'value'),
    State({'type': 'graph-dropdown', 'index': ALL}, 'value'),
    State('data-basedir', 'data')
)
def update_graphs(
        n_clicks,
        axis_range,
        binning_value,
        samples_value,
        pcr_value,
        weight_value,
        probes_value,
        graphs_values,
        data_basedir
):
    ctx = dash.callback_context
    triggerd_input = ctx.triggered[0]['prop_id'].split('.')[0]

    if n_clicks is None or n_clicks == 0:
        return None

    pp_outputs_dir = join(data_basedir, 'outputs')
    graphs_info = {}
    nb_cards = len(samples_value)
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
                'pcr': [],
                'weight': [],
                'filepaths': [],
                'size': 0,
            }
        graphs_info[graph_id]['samples'].append(samples_value[i])
        graphs_info[graph_id]['fragments'].append(probes_value[i])
        graphs_info[graph_id]['pcr'].append(pcr_value[i][-1])
        graphs_info[graph_id]['weight'].append(weight_value[i][-1])

        filedir = join(pp_outputs_dir, samples_value[i], pcr_value[i][-1], weight_value[i][-1])
        if binning_value == 0:
            filepath = join(filedir, f"{samples_value[i]}_unbinned_frequencies.tsv")
            graphs_info[graph_id]['filepaths'].append(filepath)
        else:
            filepath = join(filedir, f"{samples_value[i]}_{binning_value}kb_binned_frequencies.tsv")
            graphs_info[graph_id]['filepaths'].append(filepath)
        graphs_info[graph_id]['size'] += 1

    # TODO: use a file that stores chr data

    chr_cum_pos = list(np.cumsum(chr_pos))
    chr_boundaries = [0] + chr_cum_pos[:-1]

    figures = {}
    traces_count = 0
    for i in graphs_info:
        traces_to_add = graphs_info[i]['size']
        figures[i] = update_figure(
            graph_id=i,
            graph_dict=graphs_info[i],
            traces_colors=colors[traces_count:traces_count+traces_to_add],
            binning=binning_value,
            chr_boundaries=chr_boundaries,
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
