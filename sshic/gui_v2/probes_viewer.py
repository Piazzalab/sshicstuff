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
from dash.dependencies import Input, Output, State, ALL
import plotly.graph_objs as go


layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            dcc.Input(id='pv-number-probes', type='number', value=2, step='1',
                      placeholder='How many probes to compare :',
                      style={
                          'width': '100%',
                          'border': '1px solid #ccc',
                          'border-radius': '4px',
                          'padding': '10px',
                          'font-size': '16px',
                          'background-color': '#fff',
                          'color': '#333'
                      }),

            dbc.Tooltip("Specify the number of probes you want to compare",
                        target="pv-number-probes", className="custom-tooltip", placement="right"),
        ], width=4, className="ml-auto mt-2"),
    ]),
    html.Div(id='pv-dynamic-probes-cards', children=[],
             style={'margin-top': '20px', 'margin-bottom': '20px'})
])


@callback(
    Output('pv-dynamic-probes-cards', 'children'),
    Input('pv-number-probes', 'value'),
    State('data-basedir', 'data')
)
def update_probes_cards(n_cards, data_basedir):
    if n_cards is None or n_cards == 0:
        return []

    pp_outputs_dir = join(data_basedir, 'outputs')
    all_samples_items = sorted(os.listdir(pp_outputs_dir))
    samples_options = [{'label': s, 'value': s} for s in all_samples_items]

    probes_cards = []
    for i in range(n_cards):
        sample_card = dbc.Col(
            dbc.Card([
                dbc.CardHeader(html.Div(id={'type': 'probe-card-header', 'index': i})),
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
                                options=[{'label': "TBD", 'value': "TBD"}],
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
                            html.Div(id={'type': 'probe-card-missing-file-output', 'index': i})
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
    Output({'type': 'pcr-checkboxes', 'index': ALL}, 'options'),
    Input({'type': 'sample-dropdown', 'index': ALL}, 'value'),
    State('data-basedir', 'data')
)
def update_pcr_checkboxes_options(samples_values, data_basedir):
    nb_cards = len(samples_values)
    pp_outputs_dir = join(data_basedir, 'outputs')
    pcr_options = []
    for i in range(nb_cards):
        if samples_values[i] is None:
            pcr_options.append([])
            continue
        all_items = os.listdir(join(pp_outputs_dir, samples_values[i]))
        pcr_dirs = [item for item in all_items if os.path.isdir(join(pp_outputs_dir, samples_values[i], item))]
        pcr_options.append([{'label': d, 'value': d} for d in pcr_dirs if 'pcr' in d.lower()])
    return pcr_options


@callback(
    Output({'type': 'pcr-checkboxes', 'index': ALL}, 'value'),
    Input({'type': 'pcr-checkboxes', 'index': ALL}, 'value'),
)
def update_pcr_checkboxes(pcr_values):
    for i in range(len(pcr_values)):
        if not pcr_values[i]:
            continue
        pcr_values[i] = [pcr_values[i][-1]]
    return pcr_values


@callback(
    Output({'type': 'weight-checkboxes', 'index': ALL}, 'options'),
    Input({'type': 'pcr-checkboxes', 'index': ALL}, 'value'),
    State({'type': 'sample-dropdown', 'index': ALL}, 'value'),
    State('data-basedir', 'data')
)
def update_weight_checkboxes_options(pcr_values, samples_values, data_basedir):
    nb_cards = len(samples_values)
    pp_outputs_dir = join(data_basedir, 'outputs')
    weight_options = []
    for i in range(nb_cards):
        if samples_values[i] is None:
            weight_options.append([])
            continue
        if pcr_values[i] is None or pcr_values[i] == []:
            weight_options.append([])
            continue
        pcr_dir = join(pp_outputs_dir, samples_values[i], pcr_values[i][-1])
        weight_dirs = [w for w in os.listdir(pcr_dir) if os.path.isdir(join(pcr_dir, w))]
        weight_options.append([{'label': d, 'value': d} for d in weight_dirs])
    return weight_options


@callback(
    Output({'type': 'weight-checkboxes', 'index': ALL}, 'value'),
    Input({'type': 'weight-checkboxes', 'index': ALL}, 'value')
)
def update_weight_checkboxes(weight_values):
    for i in range(len(weight_values)):
        if not weight_values[i]:
            continue
        weight_values[i] = [weight_values[i][-1]]
    return weight_values


@callback(
    Output({'type': 'probe-card-header', 'index': ALL}, 'children'),
    Input({'type': 'sample-dropdown', 'index': ALL}, 'value'),
    Input({'type': 'probe-dropdown', 'index': ALL}, 'value'),
    Input({'type': 'pcr-checkboxes', 'index': ALL}, 'value'),
    Input({'type': 'weight-checkboxes', 'index': ALL}, 'value')
)
def update_probe_card_header(sample_values, probe_values, pcr_values, weight_values):
    nb_samples = len(sample_values)
    probes_headers = ["" for _ in range(nb_samples)]
    for i in range(nb_samples):
        if sample_values[i] is None:
            continue
        if probe_values[i] is None:
            continue
        if pcr_values[i] is None or pcr_values[i] == []:
            continue
        if weight_values[i] is None or weight_values[i] == []:
            continue

        samp_id = sample_values[i].split('/')[-1]
        probes_headers[i] = f"{samp_id} - {probe_values[i]} - {pcr_values[i][-1]} - {weight_values[i][-1]}"
    return probes_headers


