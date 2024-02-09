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


layout = dbc.Container([
    html.H2('Upload files', style={'margin-top': '20px', 'margin-bottom': '20px'}),
    dbc.Row([
        dbc.Col([
            html.Button(
                id="up-clear-list",
                className="btn btn-danger",
                children="Clear all",
            )
        ], width=2, style={'margin-top': '0px', 'margin-bottom': '20px', 'margin-left': '0px'}),
    ]),

    dbc.Row([
        dbc.Col([
            html.H6('Upload files (.tsv, .csv, .txt, .bed, etc ...)'),
            dcc.Upload(
                id="up-upload-files",
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
        ], width=5, style={'margin-top': '0px', 'margin-bottom': '25px'}),

        dbc.Col([
            html.H6('Upload samples (compressed .zip or tar.gz)'),
            dcc.Upload(
                id="up-upload-samples",
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
        ], width=7, style={'margin-top': '0px', 'margin-bottom': '25px'}),

    ]),
])




# @callback(
#     Output("dv-file-list-selector", "options"),
#     Output("dv-clear-list", "n_clicks"),
#     [Input("dv-upload-data", "filename"),
#      Input("dv-upload-data", "contents"),
#      Input("dv-clear-list", "n_clicks")],
# )
# def update_file_list(uploaded_filenames, uploaded_file_contents, n_clicks):
#     if uploaded_filenames is not None and uploaded_file_contents is not None:
#         for name, data in zip(uploaded_filenames, uploaded_file_contents):
#             save_file(name, data)
#
#     files = uploaded_files()
#     if n_clicks is not None:
#         if n_clicks > 0:
#             for filename in files:
#                 os.remove(os.path.join(TEMPORARY_DIRECTORY, filename))
#             files = []
#
#     n_clicks = 0
#     if len(files) == 0:
#         return files, n_clicks
#     else:
#         options = []
#         for filename in files:
#             options.append({'label': filename, 'value': os.path.join(TEMPORARY_DIRECTORY, filename)})
#         return options, n_clicks
#
#
# def update_table(file_path, delim):
#     if file_path and delim:
#         df = pd.read_csv(file_path, sep=delim)
#         data = df.to_dict('records')
#         columns = [{"name": i, "id": i} for i in df.columns]
#         return data, columns
#     return None, None
#
#
# @callback(
#     [Output('dv-dataframe', 'data'),
#      Output('dv-dataframe', 'columns')],
#     [Input('dv-file-list-selector', 'value'),
#      Input('dv-delim-selector', 'value')]
# )
# def update_dataframe(file_path, delim):
#     if file_path is not None and delim is not None:
#         return update_table(file_path, delim)
#     return [], []
