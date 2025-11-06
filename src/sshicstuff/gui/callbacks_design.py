import base64
import os
import subprocess
from dash import dash_table

import dash_bootstrap_components as dbc
from dash import callback, ctx, dcc, html
from dash.dependencies import Input, Output, State, ALL

from sshicstuff.core.methods import format_annealing_oligo_output, __CACHE_DIR__


def list_cached_files():
    return [{'label': f, 'value': f} for f in os.listdir(__CACHE_DIR__) if os.path.isfile(os.path.join(__CACHE_DIR__, f))]


@callback(
    Output("alert-version-o4s", "children"),
    Input("url", "pathname")
)
def show_oligo4sshic_version_alert(_):

    try:
        v = subprocess.check_output(["oligo4sshic", "--version"], text=True)
    except FileNotFoundError:
        return dbc.Alert("oligo4sshic not found in PATH. Please install it first.", color="warning", dismissable=True)
    
    v = v.split(" ")[-1]

    return dbc.Alert(f"oligo4sshic v{v} found", color="info", dismissable=True)

@callback(
    Output("alert-upload-o4s", "children"),
    Input("upload-files-o4s", "filename"),
    prevent_initial_call=True
)
def show_upload_alert(filenames):
    if filenames:
        return dbc.Alert(f"Uploaded {len(filenames)} files : {', '.join(filenames)}", color="success", dismissable=True)
    return None

@callback(
    Output("alert-clean-cache-o4s", "children"),
    Input("clear-list-o4s", "n_clicks"),
    prevent_initial_call=True
)
def show_clean_cache_alert(n_clicks):
    if n_clicks:
        return dbc.Alert("Cache cleared", color="info", dismissable=True)
    return None


@callback(
    Output("genome-fasta-dropdown", "options"),
    [Input("upload-files-o4s", "filename"),
     Input("upload-files-o4s", "contents"),
     Input("clear-list-o4s", "n_clicks")],
)
def update_dropdown(filenames, file_contents, clear_clicks):
    ctx_trigger = ctx.triggered_id

    if ctx_trigger == "clear-list-o4s":
        if os.path.exists(__CACHE_DIR__):
            for file in os.listdir(__CACHE_DIR__):
                file_path = os.path.join(__CACHE_DIR__, file)
                if os.path.isfile(file_path):
                    os.remove(file_path) 

    elif filenames and file_contents: 
        for filename, content in zip(filenames, file_contents):
            file_path = os.path.join(__CACHE_DIR__, filename)

            if os.path.exists(file_path):
                os.remove(file_path)

            data = content.split(",")[1]
            with open(file_path, "wb") as f:
                f.write(base64.b64decode(data))

    return list_cached_files()


@callback(
    Output("chromosome-region-container", "children"),
    [Input("add-chromosome-region", "n_clicks"),
     Input({'type': 'remove-region', 'index': ALL}, "n_clicks")],
    State("chromosome-region-container", "children"),
    prevent_initial_call=True
)
def update_chromosome_regions(add_clicks, remove_clicks, children):
    """
    Gère dynamiquement l'ajout et la suppression de régions chromosomiques.
    """
    ctx_trigger = ctx.triggered_id

    if children is None:
        children = []

    # Vérifier si c'est le bouton "+" qui a été cliqué
    if ctx_trigger == "add-chromosome-region":
        new_index = len(children)  # Utiliser la taille actuelle pour indexer la nouvelle ligne

 
        new_region = html.Div([
            # Première ligne : Chr, Strand et le bouton "X"
            dbc.Row([
                html.H6(f"Region {new_index + 1}"),
            ], style={'margin-top': '10px'}),

            dbc.Row([
                dbc.Col([
                    dcc.Input(
                        id={'type': 'chr', 'index': new_index},
                        placeholder="Chr",
                        type='text',
                        className="custom-input"
                    ),
                ], width=6),

                dbc.Col([
                    dcc.Dropdown(
                        id={'type': 'strand', 'index': new_index},
                        options=[{"label": "Forward", "value": "forward"},
                                 {"label": "Reverse", "value": "reverse"}],
                        placeholder="Strand",
                        multi=False
                    ),
                ], width=4),

                dbc.Col([
                    html.Button("✗ Del", id={'type': 'remove-region', 'index': new_index}, className="btn btn-danger"),
                ], width=2),
            ], style={'margin-top': '10px'}),

            # Deuxième ligne : Start et End
            dbc.Row([
                dbc.Col([
                    dcc.Input(
                        id={'type': 'start', 'index': new_index},
                        placeholder="Start",
                        type='number',
                        className="custom-input"
                    ),
                ], width=6),

                dbc.Col([
                    dcc.Input(
                        id={'type': 'end', 'index': new_index},
                        placeholder="End",
                        type='number',
                        className="custom-input"
                    ),
                ], width=6),
            ], style={'margin-top': '10px', 'margin-bottom': '10px'}),

        ], id=f"chr-region-{new_index}")

        return children + [new_region]


    # Vérifier si c'est un bouton "-" qui a été cliqué
    elif isinstance(ctx_trigger, dict) and ctx_trigger.get("type") == "remove-region":
        index_to_remove = ctx_trigger["index"]
        return [row for row in children if row["props"]["id"] != f"chr-region-{index_to_remove}"]

    return children




@callback(
    Output("secondary-site-container", "children"),
    [Input("add-secondary-site", "n_clicks"),
     Input({'type': 'remove-secondary', 'index': ALL}, "n_clicks")],
    State("secondary-site-container", "children"),
    prevent_initial_call=True
)
def update_secondary_sites(add_clicks, remove_clicks, children):
    """
    Gère dynamiquement l'ajout et la suppression des sites secondaires.
    """
    ctx_trigger = ctx.triggered_id

    if children is None:
        children = []

    if ctx_trigger == "add-secondary-site":
        new_index = len(children)

        new_site = html.Div([
            dbc.Row([
                dbc.Col([
                    dcc.Input(
                        id={'type': 'secondary-site', 'index': new_index},
                        placeholder="Enter secondary site",
                        type='text',
                        className="custom-input"
                    ),
                ], width=8),

                dbc.Col([
                    html.Button("✗ Del", id={'type': 'remove-secondary', 'index': new_index}, className="btn btn-danger"),
                ], width=2),
            ], style={'margin-top': '10px', 'margin-bottom': '10px'}),
        ], id=f"secondary-site-{new_index}")

        return children + [new_site]

    elif isinstance(ctx_trigger, dict) and ctx_trigger.get("type") == "remove-secondary":
        index_to_remove = ctx_trigger["index"]
        return [row for row in children if row["props"]["id"] != f"secondary-site-{index_to_remove}"]

    return children  



@callback(
    Output("size-value", "children"),
    Input("size", "value")
)
def update_size_value(value):
    return str(value)

@callback(
    Output("site-start-value", "children"),
    Input("site-start", "value")
)
def update_site_start_value(value):
    return str(value)

@callback(
    Output("trials-value", "children"),
    Input("trials", "value")
)
def update_trials_value(value):
    return str(value)



@callback(
    Output("alert-submit-o4s", "children"),
    Output("output-table-container", "children"),
    Input("submit-button", "n_clicks"),
    State("genome-fasta-dropdown", "value"),
    State("chromosome-region-container", "children"),
    State("site", "value"),
    State("secondary-site-container", "children"),
    State("size", "value"),
    State("site-start", "value"),
    State("np-snp-zone", "value"),
    State("complementary-size", "value"),
    State("n-snps", "value"),
    State("trials", "value"),
    prevent_initial_call=True
)
def run_oligo4sshic(
    n_clicks, 
    fasta,
    regions, 
    site, 
    secondary_sites, 
    size, 
    site_start, 
    no_snp_zone, 
    complementary_size, 
    snp_number, 
    trials
):
    
    missing_args = []
    
    if not fasta:
        missing_args.append("--fasta")
    if not site:
        missing_args.append("--site")
    if not size:
        missing_args.append("--size")
    if not site_start:
        missing_args.append("--site-start")
    if not no_snp_zone:
        missing_args.append("--no-snp-zone")
    if not complementary_size:
        missing_args.append("--complementary-size")
    if not snp_number:
        missing_args.append("--snp-number")
    if not trials:
        missing_args.append("--tries")
    if not regions:
        missing_args.append("--forward-intervals OR/AND --reverse-intervals")
    
    if missing_args:
        missing_alert = dbc.Alert(f"Missing required arguments: {', '.join(missing_args)}", color="danger", dismissable=True)
        return missing_alert, None, None

    forward_intervals = []
    reverse_intervals = []
    
    for region in regions:
        chromosome = region['props']['children'][1]['props']['children'][0]['props']['children'][0]['props']['value']
        strand = region['props']['children'][1]['props']['children'][1]['props']['children'][0]['props']['value']
        start = region['props']['children'][2]['props']['children'][0]['props']['children'][0]['props']['value']
        end = region['props']['children'][2]['props']['children'][1]['props']['children'][0]['props']['value']
        interval = f"{chromosome}:{start}-{end}"

        if strand == "forward":
            forward_intervals.append(interval)
        else:
            reverse_intervals.append(interval)

    if len(forward_intervals) > 1:
        forward_intervals = ",".join(forward_intervals)
    else:
        if forward_intervals:
            forward_intervals = forward_intervals[0]

    if len(reverse_intervals) > 1:
        reverse_intervals = ",".join(reverse_intervals)
    else:
        if reverse_intervals:
            reverse_intervals = reverse_intervals[0]

    if secondary_sites:
        secondary_sites_args = []
        for s in secondary_sites:
            upper_s = s['props']['children'][0]['props']['children'][0]['props']['children'][0]['props']['value'].upper()
            if upper_s not in secondary_sites_args:
                secondary_sites_args.append(upper_s)
        secondary_sites = ",".join(secondary_sites_args)
    else:
        secondary_sites = None

    fasta_path = os.path.join(__CACHE_DIR__, fasta)
    output_snp_path = os.path.join(__CACHE_DIR__, "output_tmp_snp.tsv")
    output_raw_path = os.path.join(__CACHE_DIR__, "output_tmp_raw.tsv")
    output_df_path = os.path.join(__CACHE_DIR__, "oligo4sshic_output.tsv")

    arguments = {}
    arguments["--fasta"] = fasta_path
    arguments["--forward-intervals"] = forward_intervals
    arguments["--reverse-intervals"] = reverse_intervals
    arguments["--output-snp"] = output_snp_path
    arguments["--output-raw"] = output_raw_path
    arguments["--site"] = site
    arguments["--secondary-sites"] = secondary_sites
    arguments["--size"] = size
    arguments["--site-start"] = site_start
    arguments["--no-snp-zone"] = no_snp_zone
    arguments["--complementary-size"] = complementary_size
    arguments["--snp-number"] = snp_number
    arguments["--tries"] = trials

    cmd = ["oligo4sshic"]
    for arg, value in arguments.items():
        if value:
            cmd.extend([arg, str(value)])  

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        success_alert = dbc.Alert(f"Execution successful! Output:\n{result.stdout}", color="success", dismissable=True)
    except subprocess.CalledProcessError as e:
        success_alert = dbc.Alert(f"Execution failed! Error:\n{e.stderr}", color="danger", dismissable=True)

    dcc.send_file(output_snp_path)
    dcc.send_file(output_raw_path)
    if not os.path.isfile(output_snp_path) or not os.path.isfile(output_raw_path):
        raise FileNotFoundError(f"Output file {output_snp_path} not found")

    df_oligo = format_annealing_oligo_output(
        output_raw_path,
        output_snp_path,
    )

    df_oligo.to_csv(output_df_path, sep="\t", index=False)

    preview_table = dash_table.DataTable(
        data=df_oligo.head(10).to_dict("records"),
        columns=[{"name": i, "id": i} for i in df_oligo.columns],
        page_size=10,
        style_table={'overflowX': 'auto'},
        style_cell={'textAlign': 'left'},
    )

    return success_alert, preview_table


@callback(
    Output("download-dataframe-tsv", "data"),
    Input("download-button", "n_clicks"),
    prevent_initial_call=True
)
def download_final_df(n_clicks):
    output_file= os.path.join(__CACHE_DIR__, "oligo4sshic_output.tsv")
    if os.path.isfile(output_file):
        return dcc.send_file(output_file)
    return None
