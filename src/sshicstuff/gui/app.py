"""
Dash application entry point for the sshicstuff GUI.

The app hosts two tabs:
- Oligo Designer  → runs oligo4sshic and generates capture tables
- ssHiC Browser   → interactive 4C-like profile viewer
"""

from __future__ import annotations

import os
from pathlib import Path

import dash
import dash_bootstrap_components as dbc
from dash import dcc, html
from flask import Flask

import sshicstuff.gui.layout_4c_profile as lb
import sshicstuff.gui.layout_design as lo

# Side-effect imports: register all Dash callbacks.
# DO NOT REMOVE even if the symbols are not referenced directly.
import sshicstuff.gui.callbacks_4c_profile  # noqa: F401
import sshicstuff.gui.callbacks_design      # noqa: F401


# ---------------------------------------------------------------------------
# Flask / Dash setup
# ---------------------------------------------------------------------------

server = Flask(__name__)
server.config["MAX_CONTENT_LENGTH"] = 1024 * 1024 * 1024  # 1 GiB upload limit
server.secret_key = os.environ.get("FLASK_SECRET_KEY", "dev-sshicstuff-secret")

_prefix = os.environ.get("SHINYPROXY_PUBLIC_PATH", "/")
_use_prefix = None if _prefix in ("", "/") else _prefix

app = dash.Dash(
    __name__,
    server=server,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    requests_pathname_prefix=_use_prefix,
    routes_pathname_prefix=_use_prefix,
)
app.config.suppress_callback_exceptions = True


# ---------------------------------------------------------------------------
# Version helper
# ---------------------------------------------------------------------------

def _get_version() -> str:
    """Read the package version from pyproject.toml, if present."""
    base = Path(__file__).resolve().parent.parent.parent.parent
    pyproject = base / "pyproject.toml"
    if pyproject.exists():
        for line in pyproject.read_text().splitlines():
            if line.strip().startswith("version"):
                return line.split("=")[1].strip().strip("\"'")
    return "unknown"


_version = _get_version()


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------

app.layout = html.Div([
    dbc.Row(
        dbc.Col(
            html.H1(f"ssDNA specific Hi-C graphical suite v{_version}"),
            width=12,
            style={"textAlign": "center", "margin": "20px"},
        )
    ),
    dcc.Tabs(
        id="tabs",
        value="oligo-tab",
        children=[
            dcc.Tab(label="Oligo Designer", value="oligo-tab"),
            dcc.Tab(label="ssHiC Browser",  value="browser-tab"),
        ],
        style={"bottom": "50px", "margin": "20px"},
    ),
    html.Div(id="tabs-content"),
])


@app.callback(
    dash.dependencies.Output("tabs-content", "children"),
    dash.dependencies.Input("tabs", "value"),
)
def render_tab(tab: str):
    return html.Div([
        html.Div(lo.layout, style={"display": "block" if tab == "oligo-tab"    else "none"}),
        html.Div(lb.layout, style={"display": "block" if tab == "browser-tab"  else "none"}),
    ])