import dash
import dash_bootstrap_components as dbc
from dash import dcc, html
from flask import Flask

# Import your page layouts here
import layout_browser
import layout_oligomaker
import browser_callbacks
import oligomaker_callbacks

server = Flask(__name__)
app = dash.Dash(__name__, server=server, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.config.suppress_callback_exceptions = True

# Layout avec Tabs
app.layout = html.Div([
    dbc.Row(
        dbc.Col(html.H1("ssDNA specific Hi-C graphical suite"), width=12, style={'textAlign': 'center', 'margin': '20px'})
    ),
    dcc.Tabs(id="tabs", value='browser-tab', children=[
        dcc.Tab(label='ssHiC Browser', value='browser-tab'),
        dcc.Tab(label='Oligo Designer', value='oligo-tab'),
    ], style={'bottom': '50px', 'margin': '20px'}),
    html.Div(id='tabs-content')
])

@app.callback(
    dash.dependencies.Output('tabs-content', 'children'),
    dash.dependencies.Input('tabs', 'value')
)
def render_tab_content(tab):
    return html.Div([
        html.Div(layout_browser.layout, style={'display': 'block' if tab == 'browser-tab' else 'none'}),
        html.Div(layout_oligomaker.layout, style={'display': 'block' if tab == 'oligo-tab' else 'none'}),
    ])


if __name__ == '__main__':
    app.run_server(debug=True)
