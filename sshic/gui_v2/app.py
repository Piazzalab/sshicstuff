from flask import Flask
import dash
import dash_bootstrap_components as dbc
from dash import html, dcc
from dash.dependencies import Input, Output, State

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
    dcc.Store(id='data-basedir'),
    dcc.Store(id='this-sample-path'),
    dcc.Store(id='this-sample-id'),
    dcc.Store(id='this-sample-out-dir-path'),
    dcc.Store(id='this-sample-ref-path'),
    dcc.Store(id='home-store'),
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


@app.callback(
    Output('home-store', 'data'),
     Input('pcr-selector', 'value'),
     Input('sample-file-selector', 'value'),
     Input('reference-selector', 'value')
)
def update_home_store(pcr_value, sample_value, reference_value):
    return {
        'pcr-selector': pcr_value,
        'sample-file-selector': sample_value,
        'reference-selector': reference_value,
    }


@app.callback(
     Output('pcr-selector', 'value'),
     Output('sample-file-selector', 'value'),
     Output('reference-selector', 'value'),
    Input('home-store', 'data'),
    State('tabs', 'value')
)
def restore_home_state(home_data, current_tab):
    if current_tab == 'home' and home_data:
        return (
            # home_data.get('data-dir-input', None),
            home_data.get('pcr-selector', None),
            home_data.get('sample-file-selector', None),
            home_data.get('reference-selector', None),
        )
    else:
        return None, None, None


if __name__ == '__main__':
    app.run_server(debug=True)
