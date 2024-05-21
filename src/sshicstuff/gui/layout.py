# dash
from dash import html, dcc
import dash_bootstrap_components as dbc
import dash_daq as daq

# common.py
from sshicstuff.gui.common import empty_figure

layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H6('Upload files'),
            dcc.Upload(
                id="upload-files",
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
                },
                multiple=True,
            ),
        ], width=5, style={'margin-top': '0px', 'margin-bottom': '25px'}),

        dbc.Col([
            html.Button(
                id="clear-list",
                className="btn btn-danger",
                children="Clear list",
            )
        ], width=2, style={'margin-top': '50px', 'margin-bottom': '20px', 'margin-left': '20px'}),
    ]),

    dbc.Row([
        # Input files dropdown COLUMN
        dbc.Col([
            html.H6('Input files dropdown'),
            dbc.Row([
                dcc.Dropdown(
                    id='oligo-dropdown',
                    placeholder="Select capture oligos file",
                    multi=False),
            ], style={'margin-top': '10px', 'margin-bottom': '20px'}),

            dbc.Row([
                dcc.Dropdown(
                    id='coord-dropdown',
                    placeholder="Select chr. coordinates file",
                    multi=False),
            ], style={'margin-top': '10px', 'margin-bottom': '20px'}),

            dbc.Row([
                dcc.Dropdown(
                    id='samples-dropdown',
                    placeholder="Select sample file",
                    multi=False
                ),
            ], style={'margin-top': '10px', 'margin-bottom': '20px'}),

            dbc.Row([
                dcc.Dropdown(
                    id='probes-dropdown',
                    placeholder="Select probe(s) or group of probes",
                    multi=True
                ),
            ], style={'margin-top': '10px', 'margin-bottom': '20px'}),
        ], width=6),

        # Region & binning settings COLUMN
        dbc.Col([
            dbc.Row([
                html.Div(
                    id='slider-output-container', style={'font-size': '14px', 'margin-bottom': '10px'}),
                dcc.Slider(
                    id='binning-slider',
                    min=0,
                    max=100,
                    step=1,
                    value=10,
                    marks={i: str(i) for i in range(0, 101, 10)},
                    included=False,
                ),
            ]),

            dbc.Row([
                dbc.Label("Region", style={'font-size': '14px', 'margin-top': '20px'}),
            ]),

            dbc.Row([
                dbc.Col([
                    dcc.Dropdown(
                        id='region-dropdown',
                        value=None,
                        placeholder="Select chromosome",
                        multi=False),

                ], width=6),

                dbc.Col([
                    dcc.Input(
                        id='start-pos',
                        placeholder="Start",
                        value=None,
                        type='text',
                        className="custom-input"
                    ),
                ], width=3),

                dbc.Col([
                    dcc.Input(
                        id='end-pos',
                        placeholder="End",
                        value=None,
                        className="custom-input"
                    ),
                ], width=3),

            ]),
        ], width=6),
    ]),

    dbc.Row([
        dbc.Col([
            html.Button(
                id="plot-button", className="plot-button", children="Plot",
                style={'margin-top': '25px'}),
        ], width=2),

        dbc.Col([
            dbc.Label("Y min", style={'font-size': '14px', 'margin-top': '10px'}),
            dcc.Input(
                id='y-min', type='number', value=None,
                placeholder='Y min', className="custom-input"),
        ], width=1),

        dbc.Col([
            dbc.Label("Y max", style={'font-size': '14px', 'margin-top': '10px'}),
            dcc.Input(
                id='y-max', type='number', value=None,
                placeholder='Y max', className="custom-input"),
        ], width=1),

        dbc.Col([
            dbc.Label("Height", style={'font-size': '14px', 'margin-top': '10px'}),
            dcc.Input(
                id='height', type='number', value=600, step=20,
                placeholder='Height', className="custom-input"),
        ], width=1),

        dbc.Col([
            dbc.Label("Width", style={'font-size': '14px', 'margin-top': '10px'}),
            dcc.Input(
                id='width', type='number', value=1600, step=20,
                placeholder='Width', className="custom-input"),
        ], width=1),

        dbc.Col([
            daq.BooleanSwitch(
                id='re-scale-switch',
                on=False,
                label='Re-scale',

            ),
            html.Div(id='re-scale-output',
                     style={'margin-top': '10px', 'margin-left': '20px', 'font-size': '12px'}),

        ], width=1, style={'margin-top': '10px', 'margin-bottom': '10px', 'margin-left': '20px'}),
    ]),

    dbc.Row([
        # html.Div(id='graphs', children=[], style={'margin-top': '20px', 'margin-bottom': '20px'}),
        dcc.Graph(
            id='graph',
            config={
                'displayModeBar': True,
                'scrollZoom': True,
                'doubleClick': 'reset',
                'autosizable': False
            },
            style={'height': '100%', 'width': '100%'},
            figure=empty_figure
        )
    ], style={'margin-left': '-50px'}),
])
