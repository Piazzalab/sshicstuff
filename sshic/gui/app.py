# app.py
import dash
import dash_bootstrap_components as dbc
from dash import html, dcc
from dash.dependencies import Input, Output

from home import layout as home_layout
from binning import layout as binning_layout

# Create a Dash application instance:
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.config.suppress_callback_exceptions = True

app.index_string = '''
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
        <style>
            body {{
                overflow-x: scroll;
            }}
        </style>
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
'''


app.layout = html.Div([
    dcc.Tabs(id="tabs", value='home', children=[
        dcc.Tab(label='Home', value='home'),
        dcc.Tab(label='Binning', value='binning'),
    ]),
    html.Div(id='page-content')
])


@app.callback(
    Output('page-content', 'children'),
    [Input('tabs', 'value')])
def display_page(value):
    if value == 'home':
        return home_layout
    elif value == 'binning':
        return binning_layout


if __name__ == '__main__':
    app.run_server(debug=True)
