import dash
import pandas as pd
from dash import callback
from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output

from sshic.core.filter import filter_contacts
from callbacks import display_sample_id, update_input_selectors

layout = dbc.Container([
    dbc.Row(
        dbc.Col([
            html.Div(id='sample-id-output', style={'margin-top': '20px', 'margin-bottom': '20px'}),
            html.Br(),
            html.Label("Select the digested fragments list:"),
            dcc.Dropdown(id='fragments-selector', multi=False),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '25px'})),

    dbc.Row(
        dbc.Col([
            html.Label("Select oligos related information table:"),
            dcc.Dropdown(id='oligo-selector', multi=False),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '30px'})),

    dbc.Row([
        dbc.Col(dash_table.DataTable(id='fragments-table'), width=4),
        dbc.Col(dash_table.DataTable(id='oligo-table'), width=4),
        dbc.Col(dash_table.DataTable(id='fragments-contacts-table'), width=4),
    ]),
])


@callback(
    Output('fragments-table', 'data'),
    Output('fragments-table', 'columns'),
    Input('fragments-selector', 'value'))
def update_fragments_table(file_path):
    if file_path is not None:
        df = pd.read_csv(file_path, sep='\t')
        data = df.to_dict('records')
        columns = [{"name": i, "id": i} for i in df.columns]
        return data, columns
    return dash.no_update, dash.no_update


@callback(
    Output('oligo-table', 'data'),
    Output('oligo-table', 'columns'),
    Input('oligo-selector', 'value'))
def update_oligo_table(file_path):
    if file_path is not None:
        df = pd.read_csv(file_path)
        data = df.to_dict('records')
        columns = [{"name": i, "id": i} for i in df.columns]
        return data, columns
    return dash.no_update, dash.no_update


# @callback(
#     Output('fragments-contacts-table', 'data'),
#     Output('fragments-contacts-table', 'columns'),
#     Input('fragments-table', 'data'),
#     Input('oligo-table', 'data'))
# def update_fragments_contacts_table(fragments_data, oligo_data):
#     if fragments_data and oligo_data:
#         fragments_df = pd.DataFrame(fragments_data)
#         oligo_df = pd.DataFrame(oligo_data)
#
#         # use your function to filter contacts
#         contacts_df = filter_contacts(fragments_df, oligo_df)
#
#         data = contacts_df.to_dict('records')
#         columns = [{"name": i, "id": i} for i in contacts_df.columns]
#         return data, columns
#     return [], []
