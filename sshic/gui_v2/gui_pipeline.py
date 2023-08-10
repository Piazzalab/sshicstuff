import os
import re
import dash
import pandas as pd
from os.path import isfile, join
import dash_bootstrap_components as dbc
from dash import callback
from dash import html, dcc, dash_table
from dash.dependencies import Input, Output, State

import utils
import filter
import coverage
import probe2fragment


def generate_data_table(id, data, columns):
    return dash_table.DataTable(
        id=id,
        data=data,
        columns=columns,
        style_table={'overflowX': 'auto'},
        page_size=8,
        style_header={
            'backgroundColor': '#eaecee',
            'color': ' #3498db ',
            'fontWeight': 'bold'},
        sort_action='native',
        sort_mode='multi',
    )


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
            html.Div(id='pp-p2f-dataframe-title',  style={'margin-top': '20px', 'margin-bottom': '20px'}),
            dcc.Loading(generate_data_table('pp-p2f-dataframe', [], []))
        ], width=4, style={'margin-top': '20px', 'margin-bottom': '25px'}),
    ]),

    dbc.Row([
        dbc.Col([
            html.Button(
                id="pp-p2f",
                className="blue-button",
                children="Probes to fragments",
            ),
            dbc.Tooltip(
                "Create a columns in the oligo table with the corresponding fragment",
                target="pp-p2f",
                className="custom-tooltip",
                placement="right",
            ),
        ], width=3, style={'margin-top': '0px', 'margin-bottom': '30px'}),

        dbc.Col([
            html.Div(id='pp-p2f-output', style={'margin-top': '20px', 'margin-bottom': '20px'}),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '30px'})
    ]),

    dbc.Row([
        dbc.Col([
            html.Button(
                id="pp-filter",
                className="blue-button",
                children="Filter",
            ),
            dbc.Tooltip(
                "This module filters the contacts by removing contacts "
                "that do not concern digested fragments containing oligos",
                target="pp-filter",
                className="custom-tooltip",
                placement="right",
            ),
        ], width=2, style={'margin-top': '0px', 'margin-bottom': '30px'}),

        dbc.Col([
            html.Div(id='pp-filter-output', style={'margin-top': '20px', 'margin-bottom': '20px'}),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '30px'})
    ]),

    dbc.Row([
        dbc.Col([
            html.Button(
                id="pp-coverage",
                className="blue-button",
                children="Coverage",
            ),
            dbc.Tooltip(
                "Calculate the coverage per oligo fragment and save the result as a bed-graph file",
                target="pp-coverage",
                className="custom-tooltip",
                placement="right",
            ),
        ], width=2, style={'margin-top': '0px', 'margin-bottom': '30px'}),

        dbc.Col([
            html.Div(id='pp-coverage-output', style={'margin-top': '20px', 'margin-bottom': '20px'}),
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


def prepare_dataframe_for_output(dataframe):
    selected_columns = ['name', 'fragment']
    df_output = dataframe[selected_columns]
    data = df_output.to_dict('records')
    columns = [{"name": col, "id": col} for col in selected_columns]
    return data, columns


@callback(
    [Output('pp-p2f', 'n_clicks'),
     Output('pp-p2f-output', 'children'),
     Output('pp-p2f-dataframe-title', 'children'),
     Output('pp-p2f-dataframe', 'data'),
     Output('pp-p2f-dataframe', 'columns')],
    [Input('pp-p2f', 'n_clicks')],
    [State('pp-fragments-selector', 'value'),
     State('pp-oligo-selector', 'value')]
)
def oligo_and_fragments(n_clicks, fragments_file, oligo_file):
    if n_clicks is None or n_clicks == 0:
        return 0, dash.no_update, dash.no_update, dash.no_update, dash.no_update

    if fragments_file is None or oligo_file is None:
        return 0, "Select both fragments and capture oligos files", dash.no_update, dash.no_update, dash.no_update

    df_oli = pd.read_csv(oligo_file, sep=utils.detect_delimiter(oligo_file))
    df_frag = pd.read_csv(fragments_file, sep=utils.detect_delimiter(fragments_file))

    title = html.H6("Oligo probes VS. Fragments ID:")
    if 'fragment' in df_oli.columns:
        data, columns = prepare_dataframe_for_output(df_oli)
        return 0, "Capture oligos table already contains a fragment column", title, data, columns
    else:
        probe2fragment.associate_probes_to_fragments(fragments_file, oligo_file)
        df_p2f = pd.read_csv(oligo_file, sep=utils.detect_delimiter(oligo_file))
        data, columns = prepare_dataframe_for_output(df_p2f)
        return 0, "Associated probes to fragments successfully", title, data, columns


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

    filter.filter_contacts(oligos_file, fragments_file, sparse_matrix, output_dir)
    return 0, "Filtered contacts file created successfully"


@callback(
    [Output('pp-coverage', 'n_clicks'),
     Output('pp-coverage-output', 'children')],
    [Input('pp-coverage', 'n_clicks')],
    [State('this-sample-out-dir-path', 'data'),
     State('this-sample-path', 'data'),
     State('pp-fragments-selector', 'value')]
)
def filter_contacts(n_clicks, output_dir, sparse_matrix, fragments_file):
    if output_dir is None:
        return dash.no_update, "Please select a sample sparse matrix (in home tab)"

    pattern = re.compile(r'.+_coverage_')
    if n_clicks == 1:
        for file in os.listdir(output_dir):
            if pattern.match(file):
                return n_clicks, "Coverage bed-graph file already exists (click again to overwrite)"

    if n_clicks is None or n_clicks == 0:
        return 0, dash.no_update

    if fragments_file is None:
        return 0, "Select a digested fragments file"

    coverage.coverage(sparse_matrix, fragments_file, output_dir)
    return 0, "Coverage file created successfully"

