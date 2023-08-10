import dash
import os
import base64
from os.path import join, dirname
import pandas as pd
from dash import callback
from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output


TEMPORARY_DIRECTORY = join(dirname(dirname(os.getcwd())), "data", "__cache__")

if not os.path.exists(TEMPORARY_DIRECTORY):
    os.makedirs(TEMPORARY_DIRECTORY)


def generate_data_table(id, data, columns):
    return dash_table.DataTable(
        id=id,
        data=data,
        columns=columns,
        style_table={'overflowX': 'auto'},
        page_size=16,
        style_header={
            'backgroundColor': '#eaecee',
            'color': ' #3498db ',
            'fontWeight': 'bold'},
        sort_action='native',
        sort_mode='multi',
    )


layout = dbc.Container([
    html.H2('Data Viewer', style={'margin-top': '20px', 'margin-bottom': '20px'}),
    dbc.Row([
        dbc.Col([
            dcc.Upload(
                id="upload-data",
                children=html.Div(
                    ["Drag and drop or click to select a file to upload."]
                ),
                style={
                    "width": "100%",
                    "height": "80px",
                    "lineHeight": "80px",
                    "borderWidth": "2px",
                    "borderStyle": "dashed",
                    "borderRadius": "20px",
                    "textAlign": "center",
                    "margin": "10px",
                },
                multiple=True,
            ),
        ], width=8, style={'margin-top': '0px', 'margin-bottom': '25px'}),
    ]),
    dbc.Row([
        dbc.Col([
            html.H4("File List"),
            dcc.Dropdown(
                id='file-list-selector',
                options=[],
                multi=False,
            ),
        ]),
        dbc.Col([
            html.H4('Select a delimiter :'),
            dcc.Dropdown(
                id='delim-selector',
                multi=False,
                options=[
                    {'label': 'Tab', 'value': '\t'},
                    {'label': 'Comma', 'value': ','},
                    {'label': 'Space', 'value': ' '},
                    {'label': 'Semicolon', 'value': ';'},
                    {'label': 'Colon', 'value': ':'}
                ],
                value=None,
            )
        ], width=3, style={'margin-top': '0px', 'margin-bottom': '0px', 'margin-left': '40px'}),

        dbc.Col([
            html.Button(
                id="clear-list",
                className="btn btn-danger",
                children="Clear List",
            )
        ], style={'margin-top': '35px', 'margin-bottom': '0px', 'margin-left': '40px'}),
    ]),
    dbc.Row([
        dbc.Col([
            html.H4("Data Preview"),
            dcc.Loading(generate_data_table('dataframe', [], []))
        ], width=12, style={'margin-top': '20px', 'margin-bottom': '25px'}),
    ]),
])


def save_file(name, content):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(os.path.join(TEMPORARY_DIRECTORY, name), "wb") as fp:
        fp.write(base64.decodebytes(data))


def uploaded_files():
    """List the files in the upload directory."""
    files = []
    for filename in os.listdir(TEMPORARY_DIRECTORY):
        path = os.path.join(TEMPORARY_DIRECTORY, filename)
        if os.path.isfile(path):
            files.append(filename)
    return files


@callback(
    Output("file-list-selector", "options"),
    Output("clear-list", "n_clicks"),
    [Input("upload-data", "filename"),
     Input("upload-data", "contents"),
     Input("clear-list", "n_clicks")],
)
def update_file_list(uploaded_filenames, uploaded_file_contents, n_clicks):
    if uploaded_filenames is not None and uploaded_file_contents is not None:
        for name, data in zip(uploaded_filenames, uploaded_file_contents):
            save_file(name, data)

    files = uploaded_files()
    if n_clicks is not None:
        if n_clicks > 0:
            for filename in files:
                os.remove(os.path.join(TEMPORARY_DIRECTORY, filename))
            files = []

    n_clicks = 0
    if len(files) == 0:
        return files, n_clicks
    else:
        options = []
        for filename in files:
            options.append({'label': filename, 'value': os.path.join(TEMPORARY_DIRECTORY, filename)})
        return options, n_clicks


def update_table(file_path, delim):
    if file_path and delim:
        df = pd.read_csv(file_path, sep=delim)
        data = df.to_dict('records')
        columns = [{"name": i, "id": i} for i in df.columns]
        return data, columns
    return None, None


@callback(
    [Output('dataframe', 'data'),
     Output('dataframe', 'columns')],
    [Input('file-list-selector', 'value'),
     Input('delim-selector', 'value')]
)
def update_dataframe(file_path, delim):
    if file_path is not None and delim is not None:
        return update_table(file_path, delim)
    return dash.no_update, dash.no_update
