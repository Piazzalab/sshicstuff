import os

import dash
import dash_bootstrap_components as dbc
from dash import dcc, html
from os.path import dirname, join
from flask import Flask
from pathlib import Path

import sshicstuff.gui.layout_4c_profile as lb
import sshicstuff.gui.layout_design as lo

# DO NOT REMOVE EVEN IF NOT USED.
import sshicstuff.gui.callbacks_4c_profile
import sshicstuff.gui.callbacks_design

server = Flask(__name__)
server.config['MAX_CONTENT_LENGTH'] = 1024 * 1024 * 1024  # 1 Go

# NEW: secret key so Flask session works (for per-user cache)
server.secret_key = os.environ.get("FLASK_SECRET_KEY", "dev-sshicstuff-secret")

prefix = os.environ.get("SHINYPROXY_PUBLIC_PATH", "/")  # fallback local

def get_app_version():
    """
    Try to read version from pyproject.toml *if it exists*
    """
    # Get the dir where THIS file is located → pas le cwd !
    base_dir = Path(__file__).resolve().parent.parent.parent.parent  # hop back to /gui/..
    pyproject = base_dir / "pyproject.toml"

    if pyproject.exists():
        with pyproject.open("r") as f:
            for line in f:
                if line.strip().startswith("version"):
                    return line.split("=")[1].strip().strip('"\'')
    else:
        return "unknown"

app = dash.Dash(
    __name__,
    server=server,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    requests_pathname_prefix=None if prefix in ("", "/") else prefix,
    routes_pathname_prefix=None if prefix in ("", "/") else prefix
)

app.config.suppress_callback_exceptions = True

version = get_app_version()

# Layout avec Tabs
app.layout = html.Div([
    dbc.Row(
        dbc.Col(html.H1(f"ssDNA specific Hi-C graphical suite V{version}"), width=12, style={'textAlign': 'center', 'margin': '20px'})
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
