import dash
import os
import base64
from urllib.parse import quote as urlquote
from os.path import join, dirname
import pandas as pd
from dash import callback
from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output


TEMPORARY_DIRECTORY = join(dirname(dirname(os.getcwd())), "data", "__cache__")

if not os.path.exists(TEMPORARY_DIRECTORY):
    os.makedirs(TEMPORARY_DIRECTORY)


layout = dbc.Container([
    html.H3('Data Viewer', style={'margin-top': '20px', 'margin-bottom': '20px'}),
    dbc.Row([
        dbc.Col([
            dcc.Upload(
                id="upload-data",
                children=html.Div(
                    ["Drag and drop or click to select a file to upload."]
                ),
                style={
                    "width": "100%",
                    "height": "60px",
                    "lineHeight": "60px",
                    "borderWidth": "1px",
                    "borderStyle": "dashed",
                    "borderRadius": "20px",
                    "textAlign": "center",
                    "margin": "10px",
                },
                multiple=True,
            ),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '25px'}),
        dbc.Col([
            html.Label('Select a delimiter :'),
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
        ], width=2, style={'margin-top': '0px', 'margin-bottom': '25px', 'margin-left': '40px'}),
    ]),
    dbc.Row([
        dbc.Col([
            html.H2("File List"),
            html.Ul(id="file-list"),
        ]),
        dbc.Col([
            html.Button(
                id="clear-list",
                className="btn btn-danger",
                children="Clear List",
                style={'margin-top': '10px', 'margin-bottom': '10px'},
            )
        ]),
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


def file_download_link(filename):
    """Create a Plotly Dash 'A' element that downloads a file from the app."""
    location = "/download/{}".format(urlquote(filename))
    return html.A(filename, href=location)


@callback(
    Output("file-list", "children"),
    Output("clear-list", "n_clicks"),
    [Input("upload-data", "filename"),
     Input("upload-data", "contents"),
     Input("clear-list", "n_clicks")],
)
def update_output(uploaded_filenames, uploaded_file_contents, n_clicks):
    """Save uploaded files, clear the file list, and regenerate the file list."""

    if uploaded_filenames is not None and uploaded_file_contents is not None:
        for name, data in zip(uploaded_filenames, uploaded_file_contents):
            save_file(name, data)

    if n_clicks is not None:
        if n_clicks > 0:
            files = uploaded_files()
            for filename in files:
                os.remove(os.path.join(TEMPORARY_DIRECTORY, filename))
            files = []

    n_clicks = 0
    files = uploaded_files()
    if len(files) == 0:
        return [html.Li("No files yet!")], n_clicks
    else:
        return [html.Li(file_download_link(filename)) for filename in files], n_clicks
