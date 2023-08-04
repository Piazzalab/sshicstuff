import os
import re
import dash
import shutil
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
        dbc.Row(
            dbc.Col([
                html.Label('Please specify the location of your data (absolute path) :',
                           style={'margin-top': '10px', 'margin-bottom': '10px'}),
                dcc.Input(
                    id='data-dir-input',
                    type='text',
                    placeholder='Input the folder path here',
                    value=join(dirname(dirname(os.getcwd())), "data"),
                    style={'width': '100%'}
                ),
            ], width=6, style={'margin-top': '50px', 'margin-bottom': '50px'}),
        ),
        dbc.Row([
            dbc.Col([
                html.Label("Select a PCR duplicates filter folder :"),
                dcc.Dropdown(
                    id='pcr-selector',
                    multi=False,
                ),
                html.Br(),
            ]),
        ]),
        dbc.Row([
            dbc.Col([
                html.Label("Select a Sample :"),
                dcc.Dropdown(
                    id='sample-file-selector',
                    options=[],
                    multi=False,
                ),
                html.Br(),
            ]),
        ]),
        dbc.Row([
            dbc.Col([
                html.Label("Please indicate a WT reference if you wish to weight your contacts :"),
                dcc.Dropdown(
                    id='reference-selector',
                    options=[],
                    value=None,
                    multi=False,
                ),
                html.Br(),
            ]),
        ]),
    ]),
])


@callback(
    Output('pcr-selector', 'options'),
    Input('data-samples-path', 'data')
)
def update_pcr_filter_selector(samples_dir_data):
    if samples_dir_data:
        pcr_filters_dir_list = [pcr for pcr in os.listdir(samples_dir_data)
                                if isdir(join(samples_dir_data, pcr)) and 'pcr' in pcr]
        return [{'label': s, 'value': s} for s in pcr_filters_dir_list]
    return dash.no_update


@callback(
    Output('sample-file-selector', 'options'),
    [Input('data-samples-path', 'data'),
     Input('pcr-selector', 'value')]
)
def update_sample_selector(samples_dir_data, pcr_value):
    if samples_dir_data and pcr_value:
        samples_dir = join(samples_dir_data, pcr_value)
        samples = sorted([
            s for s in os.listdir(samples_dir) if isfile(join(samples_dir, s))
        ],
            key=lambda x: int(re.search(r'AD(\d+)', x).group(1))
        )
        return [{'label': s, 'value': s} for s in samples]
    return dash.no_update


@callback(
    Output('reference-selector', 'options'),
    [Input('data-inputs-path', 'data'),
     Input('sample-file-selector', 'value')]
)
def update_reference_selector(inputs_dir_data, sample_file_value):
    if inputs_dir_data and sample_file_value:
        refs_dir = join(inputs_dir_data, "references")
        references = sorted([
            r for r in os.listdir(refs_dir) if isfile(join(refs_dir, r))
        ])
        return [{'label': r, 'value': r} for r in references]
    return dash.no_update


@callback(
    Output('data-samples-path', 'data'),
    Output('data-inputs-path', 'data'),
    Input('data-dir-input', 'value')
)
def get_data_samples_path(data_dir_value):
    if data_dir_value:
        samples_dir = join(data_dir_value, 'samples')
        inputs_dir = join(data_dir_value, 'inputs')
        return samples_dir, inputs_dir
    return dash.no_update, dash.no_update


@callback(
    Output('this-sample-path', 'data'),
    [Input('data-samples-path', 'data'),
     Input('pcr-selector', 'value'),
     Input('sample-file-selector', 'value')]
)
def get_samples_path_value(samples_dir_data, pcr_value, sample_path_value):
    if samples_dir_data and pcr_value and sample_path_value:
        return join(samples_dir_data, pcr_value, sample_path_value)
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
    Output('this-sample-ref-path', 'data'),
    [Input('data-inputs-path', 'data'),
     Input('reference-selector', 'value')]
)
def get_reference(inputs_dir_data, reference_value):
    if inputs_dir_data and reference_value:
        return join(inputs_dir_data, "references", reference_value)
    return None


@callback(
    Output('this-sample-out-dir-path', 'data'),
    [Input('this-sample-id', 'data'),
     Input('this-sample-path', 'data'),
     Input('this-sample-ref-path', 'data')]
)
def create_samp_dir(sample_name, sample_path_value, reference_path):
    if sample_name:
        samp_dir = join(dirname(sample_path_value), sample_name)
        samp_in_dir = join(samp_dir, "inputs")
        if not isdir(samp_dir):
            os.mkdir(samp_dir)
            os.mkdir(samp_in_dir)

        if reference_path:
            shutil.copy(reference_path, samp_in_dir)

        return samp_dir
    return dash.no_update

