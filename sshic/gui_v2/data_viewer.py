import dash
import pandas as pd
from dash import callback
from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State

from sshic.core.filter import filter_contacts
from callbacks import display_sample_id, update_input_selectors

layout = dbc.Container([
    dbc.Row([
        html.Div(id='sample-id-output', style={'margin-top': '20px', 'margin-bottom': '20px'})
    ]),
    dbc.Row([
        dbc.Col([
            html.Label("Select the digested fragments list table:"),
            dcc.Dropdown(id='fragments-list-selector', multi=False),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '25px'}),
        dbc.Col([
            html.Label("Specify a delimiter: "),
            dcc.Dropdown(
                id='frag-list-delim-selector',
                multi=False,
                options=[{'label': 'Tab', 'value': '\t'},
                         {'label': 'Comma', 'value': ','},
                         {'label': 'Space', 'value': ' '},
                         {'label': 'Semicolon', 'value': ';'},
                         {'label': 'Colon', 'value': ':'}],
                value=None,
            ),
        ], width=2, style={'margin-top': '0px', 'margin-bottom': '25px'})
    ]),

    dbc.Row([
        dbc.Col([
            html.Label("Select oligos related information table:"),
            dcc.Dropdown(id='oligo-selector', multi=False),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '30px'}),
        dbc.Col([
            html.Label("Specify a delimiter: "),
            dcc.Dropdown(
                id='oligo-delim-selector',
                multi=False,
                options=[{'label': 'Tab', 'value': '\t'},
                         {'label': 'Comma', 'value': ','},
                         {'label': 'Space', 'value': ' '},
                         {'label': 'Semicolon', 'value': ';'},
                         {'label': 'Colon', 'value': ':'}],
                value=None,
            ),
        ], width=2, style={'margin-top': '0px', 'margin-bottom': '25px'})
    ]),

    dbc.Row([
        dbc.Col([
            dcc.Loading(
                dash_table.DataTable(
                    id='fragments-list-table',
                    style_table={'overflowX': 'auto'},
                    page_size=16,
                    style_header={
                        'backgroundColor': '#eaecee',
                        'color': ' #3498db ',
                        'fontWeight': 'bold'},
                    sort_action='native',
                    sort_mode='multi',
                )
            )
        ], width=4),

        dbc.Col([
            dcc.Loading(
                dash_table.DataTable(
                    id='oligo-table',
                    style_table={'overflowX': 'auto'},
                    page_size=16,
                    style_header={
                        'backgroundColor': '#eaecee',
                        'color': ' #3498db ',
                        'fontWeight': 'bold'},
                    sort_action='native',
                    sort_mode='multi',
                )
            )
        ], width=4),
    ], style={'margin-top': '20px', 'margin-bottom': '20px'}),
], style={'max-width': '90%', 'margin': '0 auto'})


@callback(
    Output('fragments-list-table', 'data'),
    Output('fragments-list-table', 'columns'),
    [Input('fragments-list-selector', 'value'),
     Input('frag-list-delim-selector', 'value')]
)
def update_fragments_list_table(file_path, delim):
    if file_path and delim:
        df = pd.read_csv(file_path, sep=delim)
        data = df.to_dict('records')
        columns = [{"name": i, "id": i} for i in df.columns]
        return data, columns
    return None, None


@callback(
    Output('oligo-table', 'data'),
    Output('oligo-table', 'columns'),
    [Input('oligo-selector', 'value'),
     Input('oligo-delim-selector', 'value')]
)
def update_oligo_table(file_path, delim):
    if file_path and delim:
        df = pd.read_csv(file_path, sep=delim)
        data = df.to_dict('records')
        columns = [{"name": i, "id": i} for i in df.columns]
        return data, columns
    return None, None
