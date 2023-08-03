import os
import re
import dash
from os.path import join, dirname, isdir
from dash import html, dcc
from dash import callback
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output


data_dir = join(dirname(dirname(os.getcwd())), 'data', 'samples')

layout = html.Div([
    dbc.Container([
        dbc.Row([
            dbc.Col([
                html.H1('Welcome to our single strand Hi-C (sshic) platform'),
            ], style={'margin-top': '50px', 'margin-bottom': '50px'}),
        ]),
        dbc.Row([
            dbc.Col([
                html.Label("Select a PCR mode:"),
                dcc.Dropdown(
                    id='pcr-selector',
                    options=[{'label': d, 'value': d} for d in os.listdir(data_dir) if isdir(join(data_dir, d))],
                    multi=False,
                ),
                html.Br(),
            ]),
        ]),
        dbc.Row([
            dbc.Col([
                html.Label("Select a Sample:"),
                dcc.Dropdown(
                    id='sample-selector',
                    options=[],
                    multi=False,
                ),
                html.Br(),
            ]),
        ]),
    ]),
])


@callback(
    Output('sample-selector', 'options'),
    [Input('pcr-selector', 'value')]
)
def update_sample_selector(pcr_value):
    if pcr_value:
        samples_dir = join(data_dir, pcr_value)
        samples = sorted([s for s in os.listdir(samples_dir) if isdir(join(samples_dir, s))],
                         key=lambda x: int(re.search(r'AD(\d+)', x).group(1)))
        return [{'label': s, 'value': s} for s in samples]
    return dash.no_update


@callback(
    Output('sample-path', 'data'),
    [Input('pcr-selector', 'value'),
     Input('sample-selector', 'value')]
)
def get_sample_path(pcr_value, sample_value):
    if pcr_value and sample_value:
        return join(data_dir, pcr_value, sample_value)
    return dash.no_update
