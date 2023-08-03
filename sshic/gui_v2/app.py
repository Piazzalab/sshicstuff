from flask import Flask
import dash
import dash_bootstrap_components as dbc
from dash import html, dcc
from dash.dependencies import Input, Output

# Import your page layouts here
from home import layout as home_layout
from gui_binning import layout as binning_layout
from gui_fragments import layout as fragments_layout

# Create a Dash application instance:
server = Flask(__name__, template_folder='templates', )
app = dash.Dash(__name__, server=server, external_stylesheets=[dbc.themes.BOOTSTRAP, 'assets/style.css'])
app.config.suppress_callback_exceptions = True

# Set up the app layout
app.layout = html.Div([
    dcc.Tabs(id="tabs", value='home', children=[
        dcc.Tab(label='Home', value='home'),
        dcc.Tab(label='Fragments', value='fragments'),
        dcc.Tab(label='Binning', value='binning'),
    ]),
    html.Div(id='page-content'),
    dcc.Store(id='data-samples-path'),
    dcc.Store(id='data-inputs-path'),
    dcc.Store(id='sample-path'),
    dcc.Store(id='sample-id'),
    dcc.Store(id='reference-path'),
])


# Callback to update the page content based on the selected tab
@app.callback(
    Output('page-content', 'children'),
    [Input('tabs', 'value')])
def display_page(value):
    if value == 'home':
        return home_layout
    elif value == 'fragments':
        return fragments_layout
    elif value == 'binning':
        return binning_layout


if __name__ == '__main__':
    app.run_server(debug=True)
