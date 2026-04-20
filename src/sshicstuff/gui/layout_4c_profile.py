"""
Layout for the ssHiC Browser (4C-like profile viewer) tab.
"""

from dash import dcc, html
import dash_bootstrap_components as dbc
import dash_daq as daq

from sshicstuff.plot import EMPTY_FIGURE


def _card(title: str, children, id: str | None = None, extra_style: dict | None = None):
    kwargs = {"className": "section-card", "style": extra_style or {}}
    if id:
        kwargs["id"] = id
    return html.Div([
        html.Div(title, className="section-card-title"),
        *children,
    ], **kwargs)


layout = html.Div(id="browser-page", children=[
    dbc.Container(fluid=True, style={"maxWidth": "1600px"}, children=[

        # ── Alerts ───────────────────────────────────────────────────────
        html.Div(id="alert-upload-4c", className="mt-2"),
        html.Div(id="alert-clean-cache-4c"),

        # ── Row 1 : Upload ───────────────────────────────────────────────
        _card("Upload files", [
            dbc.Row([
                dbc.Col([
                    dcc.Upload(
                        id="upload-files-4c",
                        children=html.Div("Drag & drop or click — oligo, coordinates, sample files"),
                        className="upload-zone",
                        multiple=True,
                    ),
                ], width=8),
                dbc.Col([
                    html.Button(
                        "Clear cache", id="clear-list-4c",
                        className="btn btn-danger mt-1",
                    ),
                ], width=2, className="d-flex align-items-center"),
            ]),
        ]),

        # ── Row 2 : Files | View settings | Region ───────────────────────
        dbc.Row([

            # Files
            dbc.Col(width=4, children=_card("Files", [
                dcc.Dropdown(id="oligo-dropdown",
                             placeholder="Capture oligos file (optional)",
                             className="mb-2"),
                html.Div(id="oligo-loaded-badge"),
                dcc.Dropdown(id="coord-dropdown",
                             placeholder="Chromosome coordinates file",
                             className="mb-2 mt-2"),
                html.Div(id="coord-loaded-badge"),
                dcc.Dropdown(id="samples-dropdown",
                             placeholder="Sample profile file",
                             className="mb-2 mt-2"),
                html.Div(id="sample-loaded-badge"),
                dcc.Dropdown(id="probes-dropdown",
                             placeholder="Probe(s) / group(s) to display",
                             multi=True,
                             className="mt-2"),
            ])),

            # View settings
            dbc.Col(width=4, children=_card("View settings", [
                html.Div(id="binning-slider-output-container", className="slider-label"),
                dcc.Slider(
                    id="binning-slider", min=0, max=100, step=1, value=10,
                    marks={0: "0"} | {i: str(i) for i in range(10, 101, 10)},
                    included=False,
                ),
                html.Div(id="window-slider-output-container",
                         className="slider-label mt-3"),
                dcc.Slider(
                    id="window-slider", min=1, max=20, step=1, value=1,
                    marks={1: "1"} | {i: str(i) for i in range(4, 21, 4)},
                    included=False,
                ),
                html.Hr(className="my-3"),
                dbc.Label("Normalisation", className="fw-semibold mb-1",
                          style={"fontSize": "13px"}),
                dcc.RadioItems(
                    id="normalization-radio",
                    options=[
                        {"label": " Raw counts",             "value": "raw"},
                        {"label": " Fraction of viewpoint",  "value": "fraction_viewpoint"},
                        {"label": " Fraction global",        "value": "fraction_global"},
                    ],
                    value="raw",
                    labelStyle={"display": "block", "fontSize": "13px",
                                "marginBottom": "3px"},
                ),
            ])),

            # Region
            dbc.Col(width=4, children=_card("Region", [
                dbc.Label("Chromosome", style={"fontSize": "13px"}),
                dcc.Dropdown(
                    id="region-dropdown", value=None,
                    placeholder="Whole genome (default)",
                    className="mb-3",
                ),
                dbc.Row([
                    dbc.Col([
                        dbc.Label("View start (bp)", style={"fontSize": "13px"}),
                        dcc.Input(id="start-pos", placeholder="auto",
                                  value=None, type="number",
                                  className="custom-input"),
                    ], width=6),
                    dbc.Col([
                        dbc.Label("View end (bp)", style={"fontSize": "13px"}),
                        dcc.Input(id="end-pos", placeholder="auto",
                                  value=None, type="number",
                                  className="custom-input"),
                    ], width=6),
                ]),
                html.Div(id="chr-length-info",
                         className="text-muted mt-2",
                         style={"fontSize": "12px"}),
            ])),
        ], className="g-3"),

        # ── Row 3 : Display options + Plot ───────────────────────────────
        _card("Display", [
            dbc.Row([
                dbc.Col([
                    daq.BooleanSwitch(id="log-scale-switch", on=False,
                                     label="Log scale", labelPosition="right"),
                ], width=2, className="d-flex align-items-center"),

                # Advanced settings collapse
                dbc.Col([
                    html.Button("⚙ Advanced settings",
                               id="advanced-toggle",
                               className="collapse-toggle",
                               n_clicks=0,
                               style={"background": "none", "border": "none",
                                      "padding": "0", "cursor": "pointer"}),
                    dbc.Collapse(id="advanced-collapse", is_open=False, children=[
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Y min", style={"fontSize": "13px"}),
                                dcc.Input(id="y-min", type="number", value=None,
                                          placeholder="auto", className="custom-input"),
                            ], width=3),
                            dbc.Col([
                                dbc.Label("Y max", style={"fontSize": "13px"}),
                                dcc.Input(id="y-max", type="number", value=None,
                                          placeholder="auto", className="custom-input"),
                            ], width=3),
                            dbc.Col([
                                dbc.Label("Width (px)", style={"fontSize": "13px"}),
                                dcc.Input(id="width", type="number", value=1400,
                                          step=50, className="custom-input"),
                            ], width=3),
                            dbc.Col([
                                dbc.Label("Height (px)", style={"fontSize": "13px"}),
                                dcc.Input(id="height", type="number", value=600,
                                          step=50, className="custom-input"),
                            ], width=3),
                        ], className="mt-2"),
                    ]),
                ], width=6),

                dbc.Col([
                    html.Button("Plot", id="plot-button", className="plot-button"),
                ], width=2, className="d-flex align-items-center justify-content-end"),

                dbc.Col([
                    html.Button("Export PDF", id="btn-figure-pdf",
                                className="export-button me-2"),
                    html.Button("Export SVG", id="btn-figure-svg",
                                className="export-button"),
                    dcc.Download(id="download-figure-pdf"),
                    dcc.Download(id="download-figure-svg"),
                ], width=2, className="d-flex align-items-center gap-2"),
            ], align="center"),
        ]),

        # ── Row 4 : Graph ─────────────────────────────────────────────────
        dcc.Loading(
            id="loading-graph",
            type="dot",
            color="#2245b7",
            children=dcc.Graph(
                id="graph",
                config={
                    "displayModeBar": True,
                    "scrollZoom": False,
                    "doubleClick": "reset",
                    "toImageButtonOptions": {"format": "svg", "filename": "profile"},
                },
                style={"width": "100%"},
                figure=EMPTY_FIGURE,
            ),
        ),

    ], className="pb-4"),
], style={"marginTop": "16px"})