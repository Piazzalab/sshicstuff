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


def update_table(file_path, delim):
    if file_path and delim:
        df = pd.read_csv(file_path, sep=delim)
        data = df.to_dict('records')
        columns = [{"name": i, "id": i} for i in df.columns]
        return data, columns
    return None, None


layout = dbc.Container([
    html.H3('Data Viewer', style={'margin-top': '20px', 'margin-bottom': '20px'}),
    dbc.Row([
        dbc.Col([
            dcc.Upload(
                id="upload-data",
                children=html.Div(
                    ["Drag and drop or click to select a file to upload."]
                ),
                style={
                    "width": "100%",
                    "height": "60px",
                    "lineHeight": "60px",
                    "borderWidth": "1px",
                    "borderStyle": "dashed",
                    "borderRadius": "20px",
                    "textAlign": "center",
                    "margin": "10px",
                },
                multiple=True,
            ),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '25px'}),
        dbc.Col([
            html.Label('Select a delimiter :'),
            dcc.Dropdown(
                id='delim-selector',
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
        ], width=2, style={'margin-top': '0px', 'margin-bottom': '25px'}),
    ]),
    dbc.Row([
        dbc.Col([
            dcc.Loading(generate_data_table('dataframe-input', [], []))
        ], width=8),
    ]),
])


@callback(
    Output('upload-data', 'multiple'),
    Input('data-basedir', 'data')
)
def update_upload(multiple):
    return not multiple
