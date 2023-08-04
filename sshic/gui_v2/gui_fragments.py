import os
import re
import pandas as pd
import dash
from os.path import join, isfile
from os import listdir
from dash import html, dcc, callback, dash_table
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
import plotly.graph_objs as go

layout = dbc.Container([
    dbc.Row(
        dbc.Col([
            html.Div(id='sample-id-output', style={'margin-top': '20px', 'margin-bottom': '20px'}),
            html.Br(),
            html.Label("Select the digested fragments list:"),
            dcc.Dropdown(id='fragments-selector', multi=False),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '25px'})),

    dbc.Row(
        dbc.Col([
            html.Label("Select oligos related information table:"),
            dcc.Dropdown(id='oligo-selector', multi=False),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '30px'})),

    dbc.Row([
        dbc.Col(dash_table.DataTable(id='fragments-table'), width=4),
        dbc.Col(dash_table.DataTable(id='oligo-table'), width=4),
        dbc.Col(dash_table.DataTable(id='fragments-contacts-table'), width=4),
    ]),
])


@callback(
    Output('sample-id-output', 'children'),
    Input('sample-id', 'data')
)
def display_sample_id(value):
    return f"You are working on sample {value}"


@callback(
    [Output('fragments-selector', 'options'),
     Output('oligo-selector', 'options')],
    Input('data-inputs-path', 'data')
)
def update_input_selectors(inputs_dir):
    if inputs_dir:
        input_files = sorted([file for file in listdir(inputs_dir) if isfile(join(inputs_dir, file))])
        options = [{'label': file, 'value': file} for file in input_files]
        return options, options
    return dash.no_update, dash.no_update
