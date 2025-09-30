import dash
import dash_bootstrap_components as dbc
from dash import dcc, html
from flask import Flask

# Import your page layouts here
import sshicstuff.gui.browser_callbacks as bc
import sshicstuff.gui.oligomaker_callbacks as oc
import sshicstuff.gui.layout_browser as lb
import sshicstuff.gui.layout_oligomaker as lo



server = Flask(__name__)
app = dash.Dash(__name__, server=server, external_stylesheets=[dbc.themes.BOOTSTRAP])
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


if __name__ == '__main__':
    # Bind to all interfaces so Docker can publish it
    app.run_server(host="0.0.0.0", port=8050, debug=True)
