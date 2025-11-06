import os

import dash
import dash_bootstrap_components as dbc
from dash import dcc, html
from flask import Flask

import sshicstuff.gui.layout_4c_profile as lb
import sshicstuff.gui.layout_design as lo

# DO NOT REMOVE EVEN IF NOT USED.
import sshicstuff.gui.callbacks_4c_profile
import sshicstuff.gui.callbacks_design

server = Flask(__name__)
server.config['MAX_CONTENT_LENGTH'] = 1024 * 1024 * 1024  # 1 Go

prefix = os.environ.get("SHINYPROXY_PUBLIC_PATH", "/")  # fallback local

app = dash.Dash(
    __name__,
    server=server,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    requests_pathname_prefix=None if prefix in ("", "/") else prefix,
    routes_pathname_prefix=None if prefix in ("", "/") else prefix
)

app.config.suppress_callback_exceptions = True

# Layout avec Tabs
app.layout = html.Div([
    dbc.Row(
        dbc.Col(html.H1("ssDNA specific Hi-C graphical suite"), width=12, style={'textAlign': 'center', 'margin': '20px'})
    ),
    dcc.Tabs(id="tabs", value='oligo-tab', children=[
        dcc.Tab(label='Oligo Designer', value='oligo-tab'),
        dcc.Tab(label='ssHiC Browser', value='browser-tab'),
    ], style={'bottom': '50px', 'margin': '20px'}),
    html.Div(id='tabs-content')
])

@app.callback(
    dash.dependencies.Output('tabs-content', 'children'),
    dash.dependencies.Input('tabs', 'value')
)
def render_tab_content(tab):
    return html.Div([
        html.Div(lo.layout, style={'display': 'block' if tab == 'oligo-tab' else 'none'}),
        html.Div(lb.layout, style={'display': 'block' if tab == 'browser-tab' else 'none'}),
    ])
