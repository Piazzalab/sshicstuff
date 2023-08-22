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
            html.Button(id="pv-add-card-button", className="blue-button", children="Add probe"),
            dbc.Tooltip("Click to select another probe to display",
                        target="pv-add-card-button", className="custom-tooltip", placement="right"),
        ], width=2, className="ml-auto mt-2"),
    ]),
    html.Div(id='pv-dynamic-probes-cards', children=[], style={'margin-top': '20px', 'margin-bottom': '20px'})
])


@callback(
    Output('pv-dynamic-probes-cards', 'children'),
    Output('pv-add-card-button', 'n_clicks'),
    Input('pv-add-card-button', 'n_clicks')
)
def update_probes_cards(n_clicks):
    if n_clicks is None:
        return [], dash.no_update

    probes_cards = []
    for i in range(n_clicks):
        sample_card = dbc.Col(
            dbc.Card([
                dbc.CardHeader(html.Div(id={'type': 'probe-card-header', 'index': i})),
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            dcc.Dropdown(
                                options=[],
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
                                options=[
                                    {'label': 'pcrfree', 'value': 'pcrfree'},
                                    {'label': 'pcrdupkept', 'value': 'pcrdupkept'}
                                ],
                                value=[],
                                id={'type': 'pcr-checkboxes', 'index': i},
                                inline=True,
                                className='custom-checkbox-label',
                                labelStyle={"margin": "5px"}
                            )
                        ]),

                        dbc.Col([
                            dcc.Checklist(
                                options=[
                                    {'label': 'weighted', 'value': 'weighted'},
                                    {'label': 'not_weighted', 'value': 'not_weighted'}
                                ],
                                value=[],
                                id={'type': 'weight-checkboxes', 'index': i},
                                inline=True,
                                className='custom-checkbox-label',
                                labelStyle={"margin": "5px"}
                            )
                        ]),
                    ]),
                ])
            ])
        )

        probes_cards.append(sample_card)

    rows = []
    for i in range(0, len(probes_cards), 3):
        row = dbc.Row(probes_cards[i:i+3], style={'margin-top': '20px', 'margin-bottom': '20px'})
        rows.append(row)

    return rows, n_clicks


@callback(
    Output({'type': 'sample-dropdown', 'index': ALL}, 'options'),
    Input('pv-add-card-button', 'n_clicks'),
    State('data-basedir', 'data')
)
def update_sample_dropdown(n_clicks, data_basedir):
    if n_clicks is None or n_clicks == 0:
        return [dash.no_update] * dash.callback_context.inputs_list[0]["value"]

    samples_results_dir = join(data_basedir, 'outputs')
    samples_dirs_path = [join(samples_results_dir, d) for d in listdir(samples_results_dir)]
    samples_dirs_path = sorted([d for d in samples_dirs_path if os.path.isdir(d)])
    samples_id = [d.split("/")[-1] for d in samples_dirs_path]

    options = [{'label': s, 'value': p} for s, p in zip(samples_id, samples_dirs_path)]
    return [options] * dash.callback_context.inputs_list[0]["value"]


@callback(
    Output({'type': 'probe-dropdown', 'index': ALL}, 'options'),
    Input('pv-add-card-button', 'n_clicks'),
    State('data-basedir', 'data')
)
def update_probe_dropdown(n_clicks, data_basedir):
    if n_clicks is None or n_clicks == 0:
        return [dash.no_update] * dash.callback_context.inputs_list[0]["value"]

    options = [{'label': "TBD", 'value': "TBD"}]
    return [options] * dash.callback_context.inputs_list[0]["value"]


@callback(
    Output({'type': 'pcr-checkboxes', 'index': ALL}, 'value'),
    Input({'type': 'pcr-checkboxes', 'index': ALL}, 'value'),
)
def update_pcr_checkboxes(pcr_values):
    for i in range(len(pcr_values)):
        if not pcr_values[i]:
            continue
        if pcr_values[i][-1] == 'pcrfree':
            pcr_values[i] = ['pcrfree']
        elif pcr_values[i][-1] == 'pcrdupkept':
            pcr_values[i] = ['pcrdupkept']
        else:
            pcr_values[i] = []
    return pcr_values


@callback(
    Output({'type': 'weight-checkboxes', 'index': ALL}, 'value'),
    Input({'type': 'weight-checkboxes', 'index': ALL}, 'value')
)
def update_weight_checkboxes(weight_values):
    for i in range(len(weight_values)):
        if not weight_values[i]:
            continue
        if weight_values[i][-1] == 'weighted':
            weight_values[i] = ['weighted']
        elif weight_values[i][-1] == 'not_weighted':
            weight_values[i] = ['not_weighted']
        else:
            weight_values[i] = []
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
