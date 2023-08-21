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
        dbc.Col(
            dbc.Card([
                dbc.CardHeader("Samples"),
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            html.Label("Select samples : ", style={'margin-bottom': '10px'}),
                            dcc.Dropdown(id='pv-samples-selector', options=[], value=None, multi=True),
                        ], width=7, style={'margin-top': '0px', 'margin-bottom': '0px'}),

                        dbc.Col([
                            html.Label("Capture oligos table:", style={'margin-bottom': '10px'}),
                            dcc.Dropdown(id='pv-oligo-selector', multi=False),
                        ], width=5, style={'margin-top': '0px', 'margin-bottom': '0px'}),

                    ]),
                ])
            ])
        )
    ], style={'margin-top': '20px', 'margin-bottom': '50px'}),
    html.Div(id='pv-dynamic-probes-cards', children=[]),
])


@callback(
    Output('pv-samples-selector', 'options'),
    Input('selected-samples', 'data')
)
def update_samples_selector(selected_samples):
    if selected_samples is None:
        return []

    options = [{'label': s, 'value': s} for s in selected_samples]
    return options


@callback(
    Output('pv-dynamic-probes-cards', 'children'),
    Input('pv-samples-selector', 'value'),
    State('data-basedir', 'data')
)
def update_samples_cards(selected_samples, data_basedir):
    if selected_samples is None:
        return []

    samples_cards = []
    for i, sample in enumerate(selected_samples):
        sample_card = dbc.Row([
            dbc.Col(
                dbc.Card([
                    dbc.CardHeader(f"Sample {sample}"),
                    dbc.CardBody([
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
                            ], width=3),
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
                            ], width=3),
                        ]),
                        html.Br(),
                        dbc.Row([
                            dbc.Col([
                                dcc.Dropdown(
                                    options=[],
                                    value=None,
                                    placeholder="Select a probe",
                                    id={'type': 'probe-dropdown', 'index': i}
                                )
                            ], width=6),
                        ]),
                    ])
                ])
            )
        ], style={'margin-top': '10px', 'margin-bottom': '10px'})

        samples_cards.append(sample_card)
    return samples_cards


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
