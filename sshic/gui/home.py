import os
import dash
import numpy as np
from os.path import join, dirname
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
                html.Label("Select your sample(s):"),
                dcc.Input(
                    id='sample-search-input',
                    type='text',
                    placeholder='Search...',
                    style={
                        'width': '20%',
                        'border': '1px solid #ccc',
                        'border-radius': '4px',
                        'margin-left': '20px',
                        'padding': '4px',
                        'font-size': '14px',
                        'background-color': '#fff',
                        'color': '#333',
                    }
                ),
                dcc.Checklist(
                    id='select-all-checkbox',
                    options=[{"label": "Select All", "value": "All"}],
                    value=[],
                    inline=True,
                    className='custom-checkbox-label',
                    labelStyle={"margin": "5px"}
                ),

                html.Div(
                    children=[
                        dcc.Checklist(
                            id='samples-checklist',
                            options=[],
                            value=[],
                            inline=True,
                            className='custom-checkbox-label',
                            labelStyle={"margin": "5px"}
                        )
                    ],
                ),
            ], width=10,),
        ]),
    ]),
])


@callback(
    Output('samples-checklist', 'options'),
    Input('data-dir-input', 'value'),
    Input('sample-search-input', 'value')
)
def update_samples_checklist(data_value, search_value):
    if data_value:
        samples_dir = join(data_value, "samples")
        samples = sorted(np.unique([f.split("_")[0] for f in os.listdir(samples_dir)]))
        if search_value:
            samples = [s for s in samples if search_value in s]
        return [{'label': s, 'value': s} for s in samples]
    return dash.no_update


@callback(
    Output('samples-checklist', 'value'),
    Input('samples-checklist', 'options'),
    Input('select-all-checkbox', 'value'),
    prevent_initial_call=True
)
def update_checklist_selection(checkbox_options, all_selected):
    if all_selected:
        return [option['value'] for option in checkbox_options]
    return []


@callback(
    Output('data-basedir', 'data'),
    Input('data-dir-input', 'value')
)
def get_data_basedir(data_value):
    if data_value:
        return data_value
    return dash.no_update


@callback(
    Output('selected-samples', 'data'),
    Input('samples-checklist', 'value')
)
def get_selected_samples(samples_value):
    if samples_value:
        return samples_value
    return None

