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
    html.Div(id='pv-dynamic-probes-cards', children=[],
             style={'margin-top': '20px', 'margin-bottom': '20px'}),
],  fluid=True,  className='custom-container')


@callback(
    Output('pv-dynamic-probes-cards', 'children'),
    Input('selected-samples', 'data'),
    State('data-basedir', 'data')
)
def update_samples_cards(selected_samples, data_basedir):
    if selected_samples is None:
        return []

    samples_cards = []
    output_dir = join(data_basedir, "output")

    for i, sample in enumerate(selected_samples):
        sample_card = dbc.Col(
            dbc.Card([
                dbc.CardHeader(f"{sample}"),
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

                    dbc.Row([
                        dbc.Col([
                            dcc.Dropdown(
                                options=[],
                                value=None,
                                placeholder="Select probes",
                                id={'type': 'probe-dropdown', 'index': i},
                                multi=True,
                            )
                        ]),
                    ]),
                ])
            ])
        )

        samples_cards.append(sample_card)

    rows = []
    for i in range(0, len(samples_cards), 3):
        row = dbc.Row(samples_cards[i:i+3], style={'margin-top': '20px', 'margin-bottom': '20px'})
        rows.append(row)

    return rows


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
