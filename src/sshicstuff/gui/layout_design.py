"""
Layout for the Oligo Designer tab.
"""

from dash import dcc, html
import dash_bootstrap_components as dbc


def _card(title: str, children, id: str | None = None):
    kwargs = {"className": "section-card"}
    if id:
        kwargs["id"] = id
    return html.Div([
        html.Div(title, className="section-card-title"),
        *children,
    ], **kwargs)


layout = html.Div(id="oligo-page", children=[
    dbc.Container(fluid=True, style={"maxWidth": "1400px"}, children=[
        dcc.Location(id="url", refresh=False),

        # ── Alerts ──────────────────────────────────────────────────────
        html.Div(id="alert-version-o4s", className="mt-2"),
        html.Div(id="alert-upload-o4s"),
        html.Div(id="alert-clean-cache-o4s"),
        html.Div(id="alert-submit-o4s"),

        # ── Upload card ──────────────────────────────────────────────────
        _card("1 — Upload genome FASTA", [
            dbc.Row([
                dbc.Col([
                    dcc.Upload(
                        id="upload-files-o4s",
                        children=html.Div("Drag & drop or click — genome FASTA file(s)"),
                        className="upload-zone",
                        multiple=True,
                    ),
                ], width=8),
                dbc.Col([
                    html.Button("Clear cache", id="clear-list-o4s",
                                className="btn btn-danger mt-1"),
                ], width=2, className="d-flex align-items-center"),
            ], className="mb-2"),
            dcc.Dropdown(id="genome-fasta-dropdown",
                         placeholder="Select genome FASTA file",
                         className="mt-2"),
        ]),

        dbc.Row([
            # ── Chromosome regions ───────────────────────────────────────
            dbc.Col(width=6, children=_card("2 — Chromosome regions (optional)", [
                html.P(
                    "Leave empty to design oligos on the whole genome.",
                    className="text-muted mb-2",
                    style={"fontSize": "12px"},
                ),
                dcc.Store(id="region-count-store", data=1),
                html.Div(id="chromosome-region-container"),
                html.Button("+ Add region", id="add-chromosome-region",
                            className="btn btn-success btn-sm mt-2",
                            n_clicks=0),
            ])),

            # ── Enzyme sites ─────────────────────────────────────────────
            dbc.Col(width=6, children=_card("3 — Enzyme sites", [
                dbc.Row([
                    dbc.Col([
                        dbc.Label("Primary site", style={"fontSize": "13px"}),
                        dcc.Input(id="site", value="GATC",
                                  placeholder="e.g. GATC",
                                  className="custom-input"),
                    ], width=5),
                ], className="mb-3"),

                dbc.Label("Secondary site(s)", style={"fontSize": "13px"}),
                dcc.Store(id="secondary-site-count-store", data=1),
                html.Div(id="secondary-site-container"),
                html.Button("+ Add site", id="add-secondary-site",
                            className="btn btn-success btn-sm mt-2",
                            n_clicks=0),
            ])),
        ], className="g-3"),

        # ── Parameters ───────────────────────────────────────────────────
        _card("4 — Parameters", [
            dbc.Row([
                # Left: sliders
                dbc.Col(width=8, children=[
                    dbc.Row([
                        dbc.Col([
                            html.Div(id="size-label", className="slider-label"),
                            dcc.Slider(id="size", min=40, max=120, step=1, value=80,
                                       marks={i: str(i) for i in range(40, 121, 10)},
                                       included=False),
                        ], width=12),
                    ]),
                    dbc.Row([
                        dbc.Col([
                            html.Div(id="site-start-label", className="slider-label mt-2"),
                            dcc.Slider(id="site-start", min=10, max=120, step=1, value=70,
                                       marks={i: str(i) for i in range(10, 121, 10)},
                                       included=False),
                        ], width=12),
                    ]),
                    dbc.Row([
                        dbc.Col([
                            html.Div(id="trials-label", className="slider-label mt-2"),
                            dcc.Slider(id="trials", min=5, max=100, step=5, value=20,
                                       marks={i: str(i) for i in range(5, 101, 10)},
                                       included=False),
                        ], width=12),
                    ]),
                ]),

                # Right: numeric inputs
                dbc.Col(width=4, children=[
                    dbc.Label("No-SNP zone (bp)", style={"fontSize": "13px"}),
                    dcc.Input(id="np-snp-zone", value=5, type="number",
                              className="custom-input mb-3"),

                    dbc.Label("Max complementary size", style={"fontSize": "13px"}),
                    dcc.Input(id="complementary-size", value=7, type="number",
                              className="custom-input mb-3"),

                    dbc.Label("Number of SNPs", style={"fontSize": "13px"}),
                    dcc.Input(id="n-snps", value=5, type="number",
                              className="custom-input"),
                ]),
            ]),
        ]),

        # ── Submit row ───────────────────────────────────────────────────
        dbc.Row([
            dbc.Col([
                html.Button("Run oligo4sshic", id="submit-button",
                            className="custom-submit-btn"),
            ], width="auto"),
        ], className="mb-3 mt-1"),

        # ── Output preview ────────────────────────────────────────────────
        dcc.Loading(
            id="loading-design-output",
            type="dot",
            color="#2245b7",
            children=html.Div(
                id="output-table-container",
                className="section-card",
                style={"display": "none"},
                children=[
                    html.Div("5 — Output preview", className="section-card-title"),
                    html.Div(id="output-table-inner"),
                    html.Button("Download TSV", id="download-button",
                                className="custom-download-btn mt-3"),
                    dcc.Download(id="download-dataframe-tsv"),
                ],
            ),
        ),

    ], className="pb-5"),
], style={"marginTop": "16px"})