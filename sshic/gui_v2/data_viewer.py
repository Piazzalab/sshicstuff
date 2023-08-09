import dash
import os
from os.path import join, isfile
import pandas as pd
from dash import callback
from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State


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
            generate_input_selector('dataframe-input-selector', None),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '25px'}),
        dbc.Col([
            html.Label("Specify a delimiter: "),
            generate_delim_selector('dataframe-delim-input-selector'),
        ], width=2, style={'margin-top': '0px', 'margin-bottom': '25px'}),
    ]),
    dbc.Row([
        dbc.Col([
            dcc.Loading(generate_data_table('dataframe-input', [], []))
        ], width=8),
    ]),
    dbc.Row([
        dbc.Col([
            html.Label("Select an output file you'd like to visualize:"),
            generate_input_selector('dataframe-output-selector', None),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '25px'}),
        dbc.Col([
            html.Label("Specify a delimiter: "),
            generate_delim_selector('dataframe-delim-output-selector'),
        ], width=2, style={'margin-top': '0px', 'margin-bottom': '25px'}),
    ]),
    dbc.Row([
        dbc.Col(id='output-subdir-dropdown')
    ]),
    dbc.Row([
        dbc.Col([
            dcc.Loading(generate_data_table('dataframe-output', [], []))
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
    Output('dataframe-input-selector', 'options'),
    Input('data-basedir', 'data')
)
def update_input_selector(data_dir):
    if data_dir:
        inputs_dir = join(data_dir, 'inputs')
        input_files = sorted([file for file in os.listdir(inputs_dir) if isfile(join(inputs_dir, file))])
        options = [{'label': file, 'value': join(inputs_dir, file)} for file in input_files]
        return options
    return dash.no_update


@callback(
    Output('dataframe-output-selector', 'options'),
    Input('this-sample-out-dir-path', 'data')
)
def update_outputs_selector(sample_dir):
    if sample_dir:
        options = [{'label': file, 'value': os.path.join(sample_dir, file)} for file in os.listdir(sample_dir)]
        return options
    return dash.no_update


@callback(
    Output('output-subdir-dropdown', 'children'),
    Input('dataframe-output-selector', 'value')
)
def update_outputs_sub_folder_selector(value):
    if value and os.path.isdir(value):
        subdir_options = [{'label': file, 'value': os.path.join(value, file)} for file in os.listdir(value)]
        return dcc.Dropdown(id='output-subdir-selector', options=subdir_options)
    return dash.no_update



# @callback(
#     Output('subfolder-selector', 'options'),
#     Input('subfolder-selector', 'value'),
#     State('this-sample-out-dir-path', 'data')
# )
# def update_subfolder_selectors(selected_subfolder, sample_dir):
#     if selected_subfolder:
#         subfolder_path = os.path.join(sample_dir, selected_subfolder)
#         subfolder_options = [{'label': item, 'value': os.path.join(subfolder_path, item)} for item in os.listdir(subfolder_path)]
#         return subfolder_options
#     return dash.no_update


@callback(
    [Output('dataframe-input', 'data'),
     Output('dataframe-input', 'columns')],
    [Input('dataframe-input-selector', 'value'),
     Input('dataframe-delim-input-selector', 'value')]
)
def update_dataframe_1(file_path, delim):
    return update_table(file_path, delim)


@callback(
    [Output('dataframe-output', 'data'),
     Output('dataframe-output', 'columns')],
    [Input('dataframe-output-selector', 'value'),
     Input('dataframe-delim-output-selector', 'value')]
)
def update_dataframe_2(file_path, delim):
    return update_table(file_path, delim)