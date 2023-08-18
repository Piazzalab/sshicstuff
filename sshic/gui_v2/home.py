import os
import re
import dash
from os.path import join, dirname, isdir, isfile
from dash import html, dcc
from dash import callback
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output


layout = html.Div([
    dbc.Container([
        dbc.Row(
            dbc.Col(
                html.H1('Welcome to our single strand Hi-C (sshic) platform'),
                width={'size': 6, 'offset': 3},
                style={'text-align': 'center', 'margin-top': '50px', 'margin-bottom': '50px'}
            )
        ),
        dbc.Row([
            dbc.Col([
                html.Label('Please specify the location of your data (absolute path):',
                           style={'margin-top': '10px', 'margin-bottom': '10px'}),
                dcc.Input(
                    id='data-dir-input',
                    type='text',
                    placeholder='Input the folder path here',
                    value=join(dirname(dirname(os.getcwd())), "data"),
                    style={
                        'width': '100%',
                        'border': '1px solid #ccc',
                        'border-radius': '4px',
                        'padding': '6px',
                        'font-size': '14px',
                        'background-color': '#fff',
                        'color': '#333',
                    }
                ),
            ], width=6, style={'margin-top': '50px', 'margin-bottom': '50px'}),
        ]),
        dbc.Row([
            dbc.Col([
                html.Label("Select a PCR duplicates filter folder:"),
                dcc.Dropdown(
                    id='pcr-selector',
                    multi=False,
                ),
                html.Br(),
            ]),
        ]),
        dbc.Row([
            dbc.Col([
                html.Label("Select a Sample:"),
                dcc.Dropdown(
                    id='sample-file-selector',
                    options=[],
                    multi=False,
                ),
                html.Br(),
            ]),
        ]),
    ]),
])


def get_files_from_dir(directory, filter_string='', stamp="f"):
    if stamp == "f":
        return [
            f for f in os.listdir(directory)
            if isfile(join(directory, f)) and filter_string in f.lower()
        ]
    elif stamp == "d":
        return [
            f for f in os.listdir(directory)
            if isdir(join(directory, f)) and filter_string in f.lower()
        ]


@callback(
    Output('pcr-selector', 'options'),
    Input('data-dir-input', 'value')
)
def update_pcr_selector(data_value):
    if data_value:
        samples_dir = join(data_value, "samples")
        dir_list = get_files_from_dir(samples_dir, filter_string="pcr", stamp='d')
        return [{'label': s, 'value': s} for s in dir_list]
    return dash.no_update


@callback(
    Output('sample-file-selector', 'options'),
    [Input('data-dir-input', 'value'),
     Input('pcr-selector', 'value')]
)
def update_sample_selector(data_value, pcr_value):
    if data_value and pcr_value:
        samples_dir = join(data_value, "samples", pcr_value)
        samples = sorted(get_files_from_dir(samples_dir, stamp='f'),
                         key=lambda x: int(re.search(r'AD(\d+)', x).group(1)))
        return [{'label': s, 'value': s} for s in samples]
    return dash.no_update


@callback(
    Output('reference-selector', 'options'),
    [Input('data-dir-input', 'value'),
     Input('sample-file-selector', 'value')]
)
def update_reference_selector(data_value, sample_file_value):
    if data_value and sample_file_value:
        refs_dir = join(data_value, "inputs", "references")
        references = sorted(get_files_from_dir(refs_dir, stamp='f'))
        return [{'label': r, 'value': r} for r in references]
    return dash.no_update


@callback(
    Output('data-basedir', 'data'),
    Input('data-dir-input', 'value')
)
def get_data_basedir(data_value):
    if data_value:
        return data_value
    return dash.no_update


@callback(
    Output('this-sample-path', 'data'),
    [Input('data-dir-input', 'value'),
     Input('pcr-selector', 'value'),
     Input('sample-file-selector', 'value')]
)
def get_sample_path_value(data_value, pcr_value, sample_path_value):
    if data_value and pcr_value and sample_path_value:
        return join(data_value, "samples", pcr_value, sample_path_value)
    return dash.no_update


@callback(
    Output('this-sample-id', 'data'),
    Input('sample-file-selector', 'value')
)
def get_sample_id(sample_value):
    if sample_value:
        return re.search(r'AD\d+', sample_value).group()
    return dash.no_update


@callback(
    Output('this-sample-out-dir-path', 'data'),
    [Input('this-sample-id', 'data'),
     Input('this-sample-path', 'data')]
)
def create_samp_dir(sample_name, sample_path_value):
    if sample_name:
        samp_dir = join(dirname(sample_path_value), sample_name)
        samp_input_dir = join(samp_dir, "inputs")
        if not isdir(samp_dir):
            os.mkdir(samp_dir)
            os.mkdir(samp_input_dir)

        return samp_dir
    return dash.no_update
