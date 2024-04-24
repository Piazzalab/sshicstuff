from dash import html
import dash_bootstrap_components as dbc


layout = html.Div([
    dbc.Container([
        dbc.Row(
            dbc.Col(
                html.H1('Welcome to our single stranded DNA specific Hi-C (ssHiC) platform'),
                width={'size': 10, 'offset': 1},
                style={'text-align': 'center', 'margin-top': '50px', 'margin-bottom': '50px'}
            )
        )
    ]),
])
