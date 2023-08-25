import os
import re
import pandas as pd
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
                min=1,
                max=100,
                step=1,
                value=10,
                marks={1: "1"} | {i: str(i) for i in range(10, 101, 10)},
                included=False,
            ),
            html.Div(id='pv-slider-output-container',
                     style={'margin-top': '10px', 'font-size': '14px'}),
            html.Br(),
        ], width=4, style={'margin-top': '0px', 'margin-bottom': '10px', 'margin-left': '20px'}),
    ]),

    html.Div(id='pv-dynamic-probes-cards', children=[],
             style={'margin-top': '20px', 'margin-bottom': '20px'}),

    html.Div(id='pv-graphs', children=[],
             style={'margin-top': '20px', 'margin-bottom': '20px'})
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


@callback(
    Output('pv-dynamic-probes-cards', 'children'),
    Input('pv-number-probes', 'value'),
    State('data-basedir', 'data')
)
def update_probes_cards(n_cards, data_basedir):
    if n_cards is None or n_cards == 0:
        return []

    pp_outputs_dir = join(data_basedir, 'outputs')
    all_samples_items = sorted(listdir(pp_outputs_dir))
    samples_options = [{'label': s, 'value': s} for s in all_samples_items if os.path.isdir(join(pp_outputs_dir, s))]

    probes_cards = []
    for i in range(n_cards):
        sample_card = dbc.Col(
            dbc.Card([
                dbc.CardHeader(html.Div(id={'type': 'pv-probe-card-header', 'index': i})),
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            dcc.Dropdown(
                                options=samples_options,
                                value=None,
                                placeholder="Select sample",
                                id={'type': 'sample-dropdown', 'index': i},
                                multi=False,
                            )
                        ]),

                        dbc.Col([
                            dcc.Dropdown(
                                options=[],
                                value=None,
                                placeholder="Select probe",
                                id={'type': 'probe-dropdown', 'index': i},
                                multi=False,
                            )
                        ]),
                    ]),

                    dbc.Row([
                        dbc.Col([
                            dcc.Checklist(
                                options=[],
                                value=[],
                                id={'type': 'pcr-checkboxes', 'index': i},
                                inline=True,
                                className='custom-checkbox-label',
                                labelStyle={"margin": "5px"}
                            )
                        ]),

                        dbc.Col([
                            dcc.Checklist(
                                options=[],
                                value=[],
                                id={'type': 'weight-checkboxes', 'index': i},
                                inline=True,
                                className='custom-checkbox-label',
                                labelStyle={"margin": "5px"}
                            )
                        ]),
                    ]),

                    dbc.Row([
                        dbc.Col([
                            html.Div(id={'type': 'pv-display-graph-selector', 'index': i})
                        ])
                    ])
                ])
            ])
        )

        probes_cards.append(sample_card)

    rows = []
    for i in range(0, len(probes_cards), 3):
        row = dbc.Row(probes_cards[i:i+3], style={'margin-top': '20px', 'margin-bottom': '20px'})
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
    Output({'type': 'pv-display-graph-selector', 'index': MATCH}, 'children'),
    Input({'type': 'probe-dropdown', 'index': MATCH}, 'value'),
    State({'type': 'sample-dropdown', 'index': MATCH}, 'value'),
    State({'type': 'pcr-checkboxes', 'index': MATCH}, 'value'),
    State({'type': 'weight-checkboxes', 'index': MATCH}, 'value'),
    State({'type': 'sample-dropdown', 'index': ALL}, 'value')
)
def update_graph_selector(probe_value, sample_value, pcr_value, weight_value, samples_value):
    ctx = dash.callback_context
    triggerd_input = ctx.triggered[0]['prop_id'].split('.')[0]
    if triggerd_input != '':
        triggering_input_id = json.loads(triggerd_input)
        index = int(triggering_input_id['index'])
    else:
        return None

    nb_samples = len(samples_value)
    if sample_value is None or probe_value is None:
        return None
    if pcr_value is None or pcr_value == [] or weight_value is None or weight_value == []:
        return None

    graph_selector_dropdown = dbc.Row([
        dbc.Col([
            dcc.Dropdown(
                options=[{'label': f'graph {x}', 'value': f'graph {x}'} for x in range(nb_samples)],
                value=None,
                placeholder="Select graph",
                id={'type': 'graph-selector', 'index': index},
                multi=False
            )
        ])
    ])

    return graph_selector_dropdown


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
