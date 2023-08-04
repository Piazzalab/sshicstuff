import os
import re
import pandas as pd
import dash
from os.path import join
from os import listdir
from dash import html
from dash import dcc
from dash import callback
from dash import dash_table
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
import plotly.graph_objs as go


layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.Div(id='sample-id-output-container',
                     style={'margin-top': '20px', 'margin-bottom': '20px'}),
            html.Br(),
            html.Label("Select the digested fragments list :"),
            dcc.Dropdown(
                id='fragments-list-selector',
                multi=False,
            ),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '25px'}),
    ]),
    dbc.Row([
        dbc.Col([
            html.Label("Select oligos related information table :"),
            dcc.Dropdown(
                id='capture-oligo-selector',
                multi=False,
            ),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '30px'}),
    ]),
])


@callback(
    Output('sample-id-output-container', 'children'),
    [Input('sample-id', 'data')])
def update_sample_id_output(value):
    return f"You are working on sample {value}"


@callback(
    [Output('fragments-list-selector', 'options'),
     Output('capture-oligo-selector', 'options')],
    [Input('data-inputs-path', 'data')])
def update_fragments_list_selector(inputs_dir_data):
    if inputs_dir_data:
        inputs_list = sorted([f for f in listdir(inputs_dir_data) if os.path.isfile(join(inputs_dir_data, f))])
        return [{'label': f, 'value': f} for f in inputs_list], [{'label': f, 'value': f} for f in inputs_list]
    return dash.no_update, dash.no_update
