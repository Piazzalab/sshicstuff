import dash
import os
from os.path import join, isfile
import pandas as pd
from dash import callback
from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output


def generate_data_table(id, data, columns):
    return dash_table.DataTable(
        id=id,
        data=data,
        columns=columns,
        style_table={'overflowX': 'auto'},
        page_size=16,
        style_header={
            'backgroundColor': '#eaecee',
            'color': ' #3498db ',
            'fontWeight': 'bold'},
        sort_action='native',
        sort_mode='multi',
    )


def generate_input_selector(id, value):
    return dcc.Dropdown(
        id=id,
        multi=False,
        value=value,
    )


def generate_delim_selector(id):
    return dcc.Dropdown(
        id=id,
        multi=False,
        options=[
            {'label': 'Tab', 'value': '\t'},
            {'label': 'Comma', 'value': ','},
            {'label': 'Space', 'value': ' '},
            {'label': 'Semicolon', 'value': ';'},
            {'label': 'Colon', 'value': ':'}
        ],
        value=None,
    )


def update_table(file_path, delim):
    if file_path and delim:
        df = pd.read_csv(file_path, sep=delim)
        data = df.to_dict('records')
        columns = [{"name": i, "id": i} for i in df.columns]
        return data, columns
    return None, None


layout = dbc.Container([
    dbc.Row([
        html.H3('Data Viewer', style={'margin-top': '20px', 'margin-bottom': '20px'})
    ]),
    dbc.Row([
        html.Div(id='data-viewer-display-samp-id', style={'margin-top': '20px', 'margin-bottom': '20px'})
    ]),
    dbc.Row([
        dbc.Col([
            html.Label("Select an input file you'd like to visualize:"),
            generate_input_selector('dataframe-selector-1', None),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '25px'}),
        dbc.Col([
            html.Label("Specify a delimiter: "),
            generate_delim_selector('dataframe-delim-selector-1'),
        ], width=2, style={'margin-top': '0px', 'margin-bottom': '25px'}),
        dbc.Col([
            dcc.Loading(generate_data_table('dataframe-1', [], []))
        ], width=8),
    ]),
    dbc.Row([
        dbc.Col([
            html.Label("Select an input file you'd like to visualize:"),
            generate_input_selector('dataframe-selector-2', None),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '25px'}),
        dbc.Col([
            html.Label("Specify a delimiter: "),
            generate_delim_selector('dataframe-delim-selector-2'),
        ], width=2, style={'margin-top': '0px', 'margin-bottom': '25px'}),
        dbc.Col([
            dcc.Loading(generate_data_table('dataframe-2', [], []))
        ], width=8),
    ]),
])


@callback(
    Output('data-viewer-display-samp-id', 'children'),
    Input('this-sample-id', 'data')
)
def display_sample_id(value):
    return f"You are currently working on sample {value}"


@callback(
    [Output('dataframe-selector-1', 'options'),
     Output('dataframe-selector-2', 'options')],
    Input('data-basedir', 'data')
)
def update_input_selectors(data_value):
    if data_value:
        inputs_dir = join(data_value, 'inputs')
        input_files = sorted([file for file in os.listdir(inputs_dir) if isfile(join(inputs_dir, file))])
        options = [{'label': file, 'value': join(inputs_dir, file)} for file in input_files]
        return options, options
    return dash.no_update, dash.no_update


@callback(
    [Output('dataframe-1', 'data'),
     Output('dataframe-1', 'columns')],
    [Input('dataframe-selector-1', 'value'),
     Input('dataframe-delim-selector-1', 'value')]
)
def update_dataframe_1(file_path, delim):
    return update_table(file_path, delim)


@callback(
    [Output('dataframe-2', 'data'),
     Output('dataframe-2', 'columns')],
    [Input('dataframe-selector-2', 'value'),
     Input('dataframe-delim-selector-2', 'value')]
)
def update_dataframe_2(file_path, delim):
    return update_table(file_path, delim)
