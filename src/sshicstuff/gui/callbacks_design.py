"""
Dash callbacks for the Oligo Designer tab.

Wraps the `oligo4sshic` binary via the sshicstuff design command.
The binary must be installed and available on PATH.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from uuid import uuid4

import pandas as pd
import dash_bootstrap_components as dbc

from dash import callback, dash_table, dcc, html, no_update
from dash.dependencies import Input, Output, State, ALL
from flask import session

from sshicstuff.gui.cache import (
    __CACHE_DIR__ as BASE_CACHE_DIR,
    APP_INSTANCE_ID,
    save_file_cache,
    uploaded_files_cache,
)
import sshicstuff.core.design as design

SESSION_KEY = "sshicstuff_session_id"


# ---------------------------------------------------------------------------
# Session cache
# ---------------------------------------------------------------------------

def _get_user_cache_dir() -> str:
    sid = session.get(SESSION_KEY)
    if sid is None:
        sid = str(uuid4())
        session[SESSION_KEY] = sid
    cache_dir = os.path.join(str(BASE_CACHE_DIR), APP_INSTANCE_ID, sid)
    os.makedirs(cache_dir, exist_ok=True)
    return cache_dir


# ---------------------------------------------------------------------------
# Slider labels
# ---------------------------------------------------------------------------

@callback(
    Output("size-label", "children"),
    Input("size", "value"),
)
def update_size_label(v):
    return [
        "Oligo size: ",
        html.Span(f"{v} bp", style={"fontWeight": "700", "color": "#2245b7"}),
    ]


@callback(
    Output("site-start-label", "children"),
    Input("site-start", "value"),
)
def update_site_start_label(v):
    return [
        "Site start: ",
        html.Span(f"{v} bp", style={"fontWeight": "700", "color": "#2245b7"}),
    ]


@callback(
    Output("trials-label", "children"),
    Input("trials", "value"),
)
def update_trials_label(v):
    return [
        "Trials: ",
        html.Span(str(v), style={"fontWeight": "700", "color": "#2245b7"}),
    ]


# ---------------------------------------------------------------------------
# File management
# ---------------------------------------------------------------------------

@callback(
    [
        Output("genome-fasta-dropdown", "options"),
        Output("clear-list-o4s",        "n_clicks"),
        Output("alert-clean-cache-o4s", "children"),
        Output("alert-upload-o4s",      "children"),
    ],
    [
        Input("upload-files-o4s",  "filename"),
        Input("upload-files-o4s",  "contents"),
        Input("clear-list-o4s",    "n_clicks"),
    ],
)
def update_file_list(filenames, contents, n_clicks):
    cache_dir = _get_user_cache_dir()
    upload_alert = None

    if filenames and contents:
        for name, data in zip(filenames, contents):
            save_file_cache(name, data, cache_dir)
        upload_alert = dbc.Alert(
            f"Uploaded {len(filenames)} file(s): {', '.join(filenames)}",
            color="success", dismissable=True, className="py-2",
        )

    files = uploaded_files_cache(cache_dir)
    clear_alert = None

    if n_clicks and n_clicks > 0:
        for f in files:
            os.remove(os.path.join(cache_dir, f))
        files = []
        clear_alert = dbc.Alert("Cache cleared.", color="info",
                                dismissable=True, className="py-2")

    fasta_opts = [
        {"label": f, "value": os.path.join(cache_dir, f)}
        for f in files
        if f.endswith((".fa", ".fasta", ".fna", ".gz"))
    ]
    return fasta_opts, 0, clear_alert, upload_alert


# ---------------------------------------------------------------------------
# Dynamic chromosome regions
# ---------------------------------------------------------------------------

def _region_row(index: int) -> html.Div:
    return html.Div([
        dbc.Row([
            dbc.Col([
                dcc.Input(
                    id={"type": "region-input", "index": index},
                    placeholder=f"Region {index + 1}  e.g. chr1:10000-500000",
                    type="text",
                    className="custom-input",
                    style={"width": "100%"},
                ),
            ], width=10),
            dbc.Col([
                html.Button(
                    "✕", id={"type": "delete-region", "index": index},
                    className="custom-delete-button",
                    n_clicks=0,
                ),
            ], width=2, className="d-flex align-items-center"),
        ], className="mb-2 align-items-center"),
    ], id={"type": "region-row-container", "index": index})


@callback(
    [
        Output("chromosome-region-container", "children"),
        Output("region-count-store", "data"),
    ],
    [
        Input("add-chromosome-region", "n_clicks"),
        Input({"type": "delete-region", "index": ALL}, "n_clicks"),
    ],
    [
        State("region-count-store", "data"),
        State("chromosome-region-container", "children"),
    ],
    prevent_initial_call=True,
)
def manage_regions(add_clicks, delete_clicks, count, current_children):
    from dash import ctx
    triggered = ctx.triggered_id

    children = current_children or []

    if triggered == "add-chromosome-region":
        children = children + [_region_row(count)]
        return children, count + 1

    # A delete button was clicked — remove the matching row
    if isinstance(triggered, dict) and triggered.get("type") == "delete-region":
        idx = triggered["index"]
        children = [
            c for c in children
            if not (isinstance(c, dict)
                    and c.get("props", {}).get("id", {}).get("index") == idx)
        ]
        return children, count

    return no_update, no_update


# Populate the container on page load with one empty row
@callback(
    Output("chromosome-region-container", "children", allow_duplicate=True),
    Input("url", "pathname"),
    prevent_initial_call=True,
)
def init_regions(_):
    return [_region_row(0)]


# ---------------------------------------------------------------------------
# Dynamic secondary sites
# ---------------------------------------------------------------------------

def _secondary_site_row(index: int) -> html.Div:
    return html.Div([
        dbc.Row([
            dbc.Col([
                dcc.Input(
                    id={"type": "secondary-site-input", "index": index},
                    placeholder=f"Site {index + 1}  e.g. CAATTG",
                    type="text",
                    className="custom-input",
                    style={"width": "100%"},
                ),
            ], width=10),
            dbc.Col([
                html.Button(
                    "✕", id={"type": "delete-secondary-site", "index": index},
                    className="custom-delete-button",
                    n_clicks=0,
                ),
            ], width=2, className="d-flex align-items-center"),
        ], className="mb-2 align-items-center"),
    ], id={"type": "secondary-site-row-container", "index": index})


@callback(
    [
        Output("secondary-site-container",   "children"),
        Output("secondary-site-count-store", "data"),
    ],
    [
        Input("add-secondary-site", "n_clicks"),
        Input({"type": "delete-secondary-site", "index": ALL}, "n_clicks"),
    ],
    [
        State("secondary-site-count-store", "data"),
        State("secondary-site-container",   "children"),
    ],
    prevent_initial_call=True,
)
def manage_secondary_sites(add_clicks, delete_clicks, count, current_children):
    from dash import ctx
    triggered = ctx.triggered_id

    children = current_children or []

    if triggered == "add-secondary-site":
        children = children + [_secondary_site_row(count)]
        return children, count + 1

    if isinstance(triggered, dict) and triggered.get("type") == "delete-secondary-site":
        idx = triggered["index"]
        children = [
            c for c in children
            if not (isinstance(c, dict)
                    and c.get("props", {}).get("id", {}).get("index") == idx)
        ]
        return children, count

    return no_update, no_update


@callback(
    Output("secondary-site-container", "children", allow_duplicate=True),
    Input("url", "pathname"),
    prevent_initial_call=True,
)
def init_secondary_sites(_):
    defaults = ["CAATTG", "AATATT", "GANTC"]
    return [_secondary_site_row(i) for i in range(len(defaults))]


# ---------------------------------------------------------------------------
# Submit — run oligo4sshic
# ---------------------------------------------------------------------------

@callback(
    [
        Output("output-table-container", "style"),
        Output("output-table-inner",     "children"),
        Output("alert-submit-o4s",       "children"),
    ],
    Input("submit-button", "n_clicks"),
    [
        State("genome-fasta-dropdown", "value"),
        State("site",                  "value"),
        State({"type": "secondary-site-input", "index": ALL}, "value"),
        State("size",                  "value"),
        State("site-start",            "value"),
        State("np-snp-zone",           "value"),
        State("complementary-size",    "value"),
        State("n-snps",                "value"),
        State("trials",                "value"),
        State({"type": "region-input", "index": ALL}, "value"),
    ],
    prevent_initial_call=True,
)
def run_design(
    n_clicks,
    genome_path,
    primary_site,
    secondary_site_inputs,
    size,
    site_start,
    no_snp_zone,
    comp_size,
    n_snps,
    trials,
    region_inputs,
):
    if not n_clicks:
        return no_update, no_update, no_update

    if not genome_path:
        alert = dbc.Alert("Please select a genome FASTA file.", color="warning")
        return {"display": "none"}, no_update, alert

    import shutil
    if shutil.which("oligo4sshic") is None:
        alert = dbc.Alert(
            "oligo4sshic binary not found in PATH. Please install it first.",
            color="danger",
        )
        return {"display": "none"}, no_update, alert

    # Build secondary-sites string
    secondary_sites = ",".join(
        s.strip() for s in (secondary_site_inputs or []) if s and s.strip()
    ) or "CAATTG,AATATT,GANTC"

    # Parse regions into BED-format interval files if provided
    fwd_intervals_path = rev_intervals_path = None
    regions = [r.strip() for r in (region_inputs or []) if r and r.strip()]
    if regions:
        fwd_intervals_path, rev_intervals_path = _write_interval_files(regions)

    outdir = tempfile.mkdtemp(prefix="o4s_")

    cmd = [
        "oligo4sshic",
        "--fasta",            genome_path,
        "--site",             primary_site or "GATC",
        "--secondary-sites",  secondary_sites,
        "--size",             str(size or 80),
        "--site-start",       str(site_start or 70),
        "--no-snp-zone",      str(no_snp_zone or 5),
        "--complementary-size", str(comp_size or 7),
        "--snp-number",       str(n_snps or 5),
        "--tries",            str(trials or 20),
        "--output-raw",       os.path.join(outdir, "raw.fa"),
        "--output-snp",       os.path.join(outdir, "snp.fa"),
    ]
    if fwd_intervals_path:
        cmd += ["--forward-intervals", fwd_intervals_path]
    if rev_intervals_path:
        cmd += ["--reverse-intervals", rev_intervals_path]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        alert = dbc.Alert(
            [html.B("oligo4sshic failed: "), html.Pre(exc.stderr[-2000:])],
            color="danger", dismissable=True,
        )
        return {"display": "none"}, no_update, alert

    # Parse output
    raw_fa  = os.path.join(outdir, "raw.fa")
    snp_fa  = os.path.join(outdir, "snp.fa")
    if not os.path.exists(raw_fa) or not os.path.exists(snp_fa):
        alert = dbc.Alert("oligo4sshic ran but produced no output.", color="warning")
        return {"display": "none"}, no_update, alert

    try:
        df_anneal = design.format_annealing_output(raw_fa, snp_fa)
    except Exception as exc:
        alert = dbc.Alert(f"Output parsing error: {exc}", color="danger")
        return {"display": "none"}, no_update, alert

    # Store result for download
    _store_result(df_anneal, outdir)

    table = _make_table(df_anneal)
    alert = dbc.Alert(
        f"Done — {len(df_anneal)} oligo(s) designed.", color="success", dismissable=True
    )
    return {"display": "block"}, table, alert


def _write_interval_files(regions: list[str]) -> tuple[str, str]:
    """Convert 'chr:start-end' strings to two temporary BED files (fwd / rev)."""
    rows = []
    for r in regions:
        r = r.replace(",", "")
        if ":" in r:
            chrom, rest = r.split(":", 1)
            start_s, end_s = rest.split("-")
            rows.append((chrom.strip(), int(start_s.strip()), int(end_s.strip())))
        else:
            rows.append((r, None, None))

    fwd = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
    rev = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
    for chrom, start, end in rows:
        line = f"{chrom}\t{start or 0}\t{end or ''}\n"
        fwd.write(line)
        rev.write(line)
    fwd.close()
    rev.close()
    return fwd.name, rev.name


# Module-level result store (single-user assumption for local GUI)
_last_result: pd.DataFrame | None = None


def _store_result(df: pd.DataFrame, outdir: str) -> None:
    global _last_result
    _last_result = df


def _make_table(df: pd.DataFrame) -> dash_table.DataTable:
    cols = [{"name": c, "id": c} for c in df.columns]
    return dash_table.DataTable(
        data=df.head(200).to_dict("records"),
        columns=cols,
        page_size=20,
        style_table={"overflowX": "auto"},
        style_header={
            "backgroundColor": "#eef1fb",
            "fontWeight": "700",
            "fontSize": "12px",
            "color": "#2245b7",
            "borderBottom": "2px solid #e0e4ef",
        },
        style_cell={
            "fontSize": "12px",
            "padding": "5px 10px",
            "border": "none",
            "borderBottom": "1px solid #f0f2f8",
            "maxWidth": "200px",
            "overflow": "hidden",
            "textOverflow": "ellipsis",
        },
        style_data_conditional=[
            {"if": {"row_index": "odd"}, "backgroundColor": "#fafbff"},
        ],
        tooltip_data=[
            {col["id"]: {"value": str(row.get(col["id"], "")), "type": "markdown"}
             for col in cols}
            for row in df.head(200).to_dict("records")
        ],
        tooltip_duration=None,
    )


# ---------------------------------------------------------------------------
# Download TSV
# ---------------------------------------------------------------------------

@callback(
    Output("download-dataframe-tsv", "data"),
    Input("download-button", "n_clicks"),
    prevent_initial_call=True,
)
def download_tsv(n_clicks):
    if not n_clicks or _last_result is None:
        return no_update
    return dcc.send_data_frame(_last_result.to_csv, "oligo_design.tsv",
                               sep="\t", index=False)