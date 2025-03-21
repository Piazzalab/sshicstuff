from flask import Flask
import dash
import dash_bootstrap_components as dbc
from dash import html

# Import your page layouts here
import sshicstuff.gui.browser as browser
import sshicstuff.gui.layout as layout

# Create a Dash application instance:
server = Flask(__name__, template_folder='templates')
app = dash.Dash(__name__, server=server, external_stylesheets=[dbc.themes.BOOTSTRAP, 'assets/style.css'])
app.config.suppress_callback_exceptions = True

# Set up the app layout
app.layout = html.Div([
    dbc.Row(
        dbc.Col(
            html.H1('Single stranded DNA specific Hi-C (ssHiC) browser'),
            width={'size': 10, 'offset': 1},
            style={'text-align': 'center', 'margin-top': '50px', 'margin-bottom': '50px'}
        )
    ),
    layout.layout,
])


if __name__ == '__main__':
    app.run_server(debug=True)
