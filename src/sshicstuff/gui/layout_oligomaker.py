# dash
import dash_bootstrap_components as dbc
from dash import html, dcc

layout = html.Div(id="oligo-page", children=[
    dbc.Container([
        dcc.Location(id="url", refresh=False),
        dbc.Row([
            html.Div(id="alert-version"),
            html.Div(id="alert-upload"),
            html.Div(id="alert-clean-cache"),
            html.Div(id="alert-submit"),
        ]),

        dbc.Row([
            dbc.Col([
                html.H6('Upload files'),
                dcc.Upload(
                    id="upload-files-oligo",
                    children=html.Div(
                        ["Drag and drop your genome fasta files."]
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
            ], width=4, style={'margin-top': '0px', 'margin-bottom': '20px'}),

            dbc.Col([
                html.Button(
                    id="clear-list-oligo",
                    className="btn btn-danger",
                    children="Clear list",
                )
            ], width=2, style={'margin-top': '50px', 'margin-bottom': '20px'}),

            dbc.Col([
                html.Button("Submit", id="submit-button", className="btn btn-primary custom-submit-btn")
            ], width=2, style={'text-align': 'center', 'margin-top': '10px', 'margin-bottom': '10px'}),
        ]),

        dbc.Row([
            # First column
            dbc.Col([
                dbc.Row([
                    dbc.Col([
                        html.H6('1. Select genome fasta file'),
                        dcc.Dropdown(
                            id='genome-fasta-dropdown',
                            placeholder="Select genome",
                            multi=False),
                    ], width=10),
                ], style={'margin-top': '20px'}),

                # Dynamic Chromosome Region Selection
                dbc.Row([
                    dbc.Col([
                        dbc.Row([
                            html.H6("2. Select chromosome region(s)"),
                            html.Div(id="chromosome-region-container"),  # Dynamic container
                        ]),

                        dbc.Row([
                            dbc.Col([
                                html.Button("+ Add ", id="add-chromosome-region", className="btn btn-success", n_clicks=0),
                            ], width=2),
                        ]),
                    ], width=10),
                ], style={'margin-top': '20px'}),

                dbc.Row([
                    dbc.Col([
                        html.H6('3. Primary site'),
                    ], width=4),

                    dbc.Col([
                        dcc.Input(
                            id='site',
                            placeholder="site",
                            value="GATC",
                            className="custom-input"
                        ),
                    ], width=5),
                ], style={'margin-top': '20px'}),


                # Dynamic Secondary Site Selection
                dbc.Row([
                    dbc.Col([
                        dbc.Row([
                            html.H6("4. Secondary site(s)"),
                            html.Div(id="secondary-site-container"),  # Conteneur dynamique
                        ]),

                        dbc.Row([
                            dbc.Col([
                                html.Button("+ Add", id="add-secondary-site", className="btn btn-success", n_clicks=0),
                            ], width=2),
                        ]),
                    ], width=10),
                ], style={'margin-top': '20px'}),

            ], width=6),

            # Second column
            dbc.Col([
                dbc.Row([
                    dbc.Col([
                        html.Div([
                            html.H6("5. Oligo size: ", style={'display': 'inline', 'text-align': 'center'}),
                            html.Span(id="size-value", children="80", style={'font-weight': 'bold'})
                        ]),
                        dcc.Slider(
                            id='size',
                            min=40,
                            max=100,
                            step=1,
                            value=80,
                            marks={i: str(i) for i in range(40, 101, 5)},
                            included=False,
                        ),
                    ], width=12),
                ], style={'margin-top': '10px'}),

                dbc.Row([
                    dbc.Col([
                        html.Div([
                            html.H6("6. Site start: ", style={'display': 'inline', 'text-align': 'center'}),
                            html.Span(id="site-start-value", children="70", style={'font-weight': 'bold'})
                        ]),
                        dcc.Slider(
                            id='site-start',
                            min=40,
                            max=100,
                            step=1,
                            value=70,
                            marks={i: str(i) for i in range(40, 101, 5)},
                            included=False,
                        ),
                    ], width=12),
                ], style={'margin-top': '10px'}),

                dbc.Row([
                    dbc.Col([html.H6('7. No SNP zone', style={'display': 'inline', 'text-align': 'center'})], width=4),
                    dbc.Col([html.H6('8. Max oligo complementary size', style={'display': 'inline', 'text-align': 'center'})], width=4),
                    dbc.Col([html.H6('9. Number of SNPs', style={'display': 'inline', 'text-align': 'center'})], width=4),
                ], style={'margin-top': '20px'}),

                dbc.Row([
                    dbc.Col([
                        dcc.Input(
                            id='np-snp-zone',
                            placeholder="No SNP zone",
                            value=5,
                            type='number',
                            className="custom-input"
                        ),
                    ], width=4),
                    dbc.Col([
                        dcc.Input(
                            id='complementary-size',
                            placeholder="Complementary size",
                            value=7,
                            type='number',
                            className="custom-input"
                        ),
                    ], width=4),
                    dbc.Col([
                        dcc.Input(
                            id='n-snps',
                            placeholder="Number of SNPs",
                            value=5,
                            type='number',
                            className="custom-input"
                        ),
                    ], width=4),
                ], style={'margin-top': '10px'}),

                dbc.Row([
                    dbc.Col([
                        html.Div([
                            html.H6("10. Trials: ", style={'display': 'inline'}),
                            html.Span(id="trials-value", children="20", style={'font-weight': 'bold'})
                        ]),
                        dcc.Slider(
                            id='trials',
                            min=10,
                            max=100,
                            step=5,
                            value=20,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            included=False,
                        ),
                    ], width=12),
                ], style={'margin-top': '20px'}),
            ]),
        ]),

        dbc.Row([
            dbc.Col([
                html.Hr(),
                html.H4("Output Preview"),
                html.Div(id="output-table-container"),
                html.Button("Download TSV", id="download-button", className="btn btn-primary custom-download-btn"),
                dcc.Download(id="download-dataframe-tsv")
            ], width=12)
        ])
    ])
], style={'margin-top': '20px', 'margin-bottom': '20px'})
