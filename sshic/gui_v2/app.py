from flask import Flask
import dash
import dash_bootstrap_components as dbc
from dash import html, dcc
from dash.dependencies import Input, Output

# Import your page layouts here
import home
import data_viewer
import gui_fragments
import gui_binning

# Create a Dash application instance:
server = Flask(__name__, template_folder='templates', )
app = dash.Dash(__name__, server=server, external_stylesheets=[dbc.themes.BOOTSTRAP, 'assets/style.css'])
app.config.suppress_callback_exceptions = True

# Set up the app layout
app.layout = html.Div([
    dcc.Tabs(id="tabs",
             value='home',
             parent_className='custom-tabs',
             className='custom-tabs-container',
             children=[
                 dcc.Tab(label='Home', value='home',
                         className='custom-tab', selected_className='custom-tab--selected'),
                 dcc.Tab(label='Data Viewer', value='data-viewer',
                         className='custom-tab', selected_className='custom-tab--selected'),
                 dcc.Tab(label='Organize Contacts', value='fragments',
                         className='custom-tab', selected_className='custom-tab--selected'),
                 dcc.Tab(label='Bin Contacts', value='binning',
                         className='custom-tab', selected_className='custom-tab--selected'),
             ]),
    html.Div(id='page-content'),
    dcc.Store(id='data-samples-path'),
    dcc.Store(id='data-inputs-path'),
    dcc.Store(id='this-sample-path'),
    dcc.Store(id='this-sample-id'),
    dcc.Store(id='this-sample-out-dir-path'),
    dcc.Store(id='this-sample-ref-path'),
])


# Callback to update the page content based on the selected tab
@app.callback(
    Output('page-content', 'children'),
    [Input('tabs', 'value')])
def display_page(value):
    if value == 'home':
        return home.layout
    elif value == 'data-viewer':
        return data_viewer.layout
    elif value == 'fragments':
        return gui_fragments.layout
    elif value == 'binning':
        return gui_binning.layout


if __name__ == '__main__':
    app.run_server(debug=True)
