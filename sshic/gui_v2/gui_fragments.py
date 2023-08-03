import re
import pandas as pd
import dash
from os.path import join
from os import listdir
from dash import html
from dash import dcc
from dash import callback
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
import plotly.graph_objs as go


layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.Div(id='sample-id-output-container'),
            html.Br(),
            html.Label("Select the digested fragments list :",
                       style={'margin-top': '20px', 'margin-bottom': '20px'}),
            dcc.Dropdown(
                id='fragments-list-selector',
                multi=False,
            ),
            html.Label("Select oligos related information table :",
                       style={'margin-top': '20px', 'margin-bottom': '20px'}),
            dcc.Dropdown(
                id='capture-oligo-selector',
                multi=False,
            ),
        ], width=4, style={'position': 'absolute', 'top': '100px', 'left': '25px'}),
    ]),
])


@callback(
    Output('sample-id-output-container', 'children'),
    [Input('sample-id', 'data')])
def update_output(value):
    return f"You are working on sample {value}"
