from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc

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
