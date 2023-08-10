import os
import re
import dash
import pandas as pd
from os.path import isfile, join
from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc
from dash import callback
from dash import html, dcc, dash_table
from dash.dependencies import Input, Output, State

import filter


layout = dbc.Container([
    dbc.Row([
        html.Div(id='pp-sample-id-output',  style={'margin-top': '20px', 'margin-bottom': '20px'}),
    ]),
    dbc.Row([
        dbc.Col([
            html.Label("Digested fragments list:"),
            dcc.Dropdown(id='pp-fragments-selector', multi=False),
        ], width=4, style={'margin-top': '0px', 'margin-bottom': '25px'}),

        dbc.Col([
            html.Label("Capture oligos table:"),
            dcc.Dropdown(id='pp-oligo-selector', multi=False),
        ], width=4, style={'margin-top': '0px', 'margin-bottom': '30px'}),

        dbc.Col([
            html.Label("Chromosome coordinates:"),
            dcc.Dropdown(id='pp-chr-coords', multi=False),
        ], width=4, style={'margin-top': '0px', 'margin-bottom': '30px'})
    ]),

    dbc.Row([
        dbc.Col([
            html.Button(
                id="pp-filter",
                className="blue-button",
                children="Filter",
            ),
            dbc.Tooltip(
                "Filter the contacts based on the oligos and fragments data, "
                "and save the filtered contacts to a TSV file.",
                target="pp-filter",
                className="custom-tooltip",
                placement="right",
            ),
        ], width=1, style={'margin-top': '0px', 'margin-bottom': '30px'}),

        dbc.Col([
            html.Div(id='pp-filter-output', style={'margin-top': '20px', 'margin-bottom': '20px'}),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '30px'})
    ])
])


@callback(
    Output('pp-sample-id-output', 'children'),
    Input('this-sample-id', 'data')
)
def display_sample_id(sample_id):
    return f"You are working on sample : {sample_id}"


@callback(
    Output('pp-fragments-selector', 'options'),
    Output('pp-oligo-selector', 'options'),
    Output('pp-chr-coords', 'options'),
    Input('data-basedir', 'data')
)
def update_dropdowns(data_basedir):
    if data_basedir is None:
        return [], [], []
    inputs_dir = join(data_basedir, "inputs")
    inputs_files = sorted([f for f in os.listdir(inputs_dir) if isfile(join(inputs_dir, f))])
    options = [{'label': f, 'value': join(inputs_dir, f)} for f in inputs_files]
    return options, options, options


@callback(
    [Output('pp-filter', 'n_clicks'),
     Output('pp-filter-output', 'children')],
    [Input('pp-filter', 'n_clicks')],
    [State('this-sample-out-dir-path', 'data'),
     State('this-sample-path', 'data'),
     State('pp-fragments-selector', 'value'),
     State('pp-oligo-selector', 'value')]
)
def filter_contacts(n_clicks, output_dir, sparse_matrix, fragments_file, oligos_file):
    if output_dir is None:
        return dash.no_update, "Please select a sample sparse matrix (in home tab)"

    pattern = re.compile(r'.+_filtered\.tsv')
    if n_clicks == 1:
        for file in os.listdir(output_dir):
            if pattern.match(file):
                return n_clicks, "Filtered contacts file already exists (click again to overwrite)"

    if n_clicks is None or n_clicks == 0:
        return 0, dash.no_update

    if fragments_file is None:
        return 0, "Select a digested fragments file"
    if oligos_file is None:
        return 0, "Select a capture oligos file"
    if sparse_matrix is None:
        return 0, "Please select a sample sparse matrix (in home tab)"

    filter.filter_contacts(oligos_file, fragments_file, sparse_matrix, output_dir)
    return 0, "Filtered contacts file created successfully"
