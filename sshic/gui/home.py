import os
import re
from os.path import join
from dash import html, dcc
from dash import callback
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output


data_dir = "../../data/samples"

layout = html.Div([
    dbc.Container([
        dbc.Row([
            dbc.Col([
                html.H1('Welcome to our single strand Hi-C (sshic) platform'),
                html.P('N.B : This platform is only a visualisation tool, it does not perform any analysis'),
            ], style={'margin-top': '50px', 'margin-bottom': '50px'}),
        ]),
        dbc.Row([
            dbc.Col([
                html.Label("Select a PCR mode:"),
                dcc.Dropdown(
                    id='pcr-selector',
                    options=[{'label': d, 'value': d} for d in os.listdir(data_dir) if os.path.isdir(join(data_dir, d))],
                    value=os.listdir(data_dir)[0],
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
    samples_dir = os.path.join(data_dir, pcr_value)
    samples = sorted([s for s in os.listdir(samples_dir) if os.path.isdir(join(samples_dir, s))],
                     key=lambda x: int(re.search(r'AD(\d+)', x).group(1)))
    return [{'label': s, 'value': s} for s in samples]
