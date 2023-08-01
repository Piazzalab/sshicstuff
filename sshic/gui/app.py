# app.py
import dash
import dash_bootstrap_components as dbc
from dash import html, dcc
from dash.dependencies import Input, Output

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
    dcc.Tabs(id="tabs", value='tab-1', children=[
        dcc.Tab(label='Binning', value='tab-1'),
    ]),
    html.Div(id='tabs-content')
])


@app.callback(
    Output('tabs-content', 'children'),
    Input('tabs', 'value'))
def render_content(tab):
    if tab == 'tab-1':
        return binning_layout


if __name__ == '__main__':
    app.run_server(debug=True)
