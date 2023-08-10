import os
import pandas as pd
from os.path import isfile, join
from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc
from dash import callback
from dash import html, dcc, dash_table
from dash.dependencies import Input, Output


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


