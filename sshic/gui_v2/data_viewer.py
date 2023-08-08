import dash
import os
from os.path import join, isfile, dirname
import pandas as pd
from dash import callback
from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State

from file_tree import FileTree

filetree = FileTree(join(dirname(dirname(os.getcwd())), "data"))
layout = dbc.Container([
    dbc.Row([
        html.H3('Data Viewer', style={'margin-top': '20px', 'margin-bottom': '20px'})
    ]),
    dbc.Row([
        dbc.Col([
            html.Label("Select an input file you'd like to visualize:"),
            filetree.render(),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '25px'}),
    ]),
])
