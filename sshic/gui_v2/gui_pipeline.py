import os
import re
import dash
import pandas as pd
from os.path import basename, join, isfile
from shutil import copyfile
import dash_bootstrap_components as dbc
from dash import callback
from dash import html, dcc, dash_table
from dash.dependencies import Input, Output, State

import utils
import filter
import coverage
import probe2fragment


def generate_data_table(id, data, columns, rows):
    return dash_table.DataTable(
        id=id,
        data=data,
        columns=columns,
        style_cell={'textAlign': 'left'},
        style_table={'overflowX': 'auto'},
        page_size=rows,
        style_header={
            'backgroundColor': '#eaecee',
            'color': ' #3498db ',
            'fontWeight': 'bold'},
        sort_action='native',
        sort_mode='multi',
    )


layout = dbc.Container([
    dbc.Row([
        html.Div(id='pp-sample-id-output',  style={'margin-top': '20px', 'margin-bottom': '20px'}),
    ]),

    dbc.Row([
        dbc.Col([
            html.Label("Digested fragments list:"),
            dcc.Dropdown(id='pp-fragments-selector', multi=False),
        ], width=4, style={'margin-top': '0px', 'margin-bottom': '25px'}),

        dbc.Col([
            html.Label("Capture oligos table:"),
            dcc.Dropdown(id='pp-oligo-selector', multi=False),
        ], width=4, style={'margin-top': '0px', 'margin-bottom': '30px'}),

        dbc.Col([
            html.Label("Chromosome coordinates:"),
            dcc.Dropdown(id='pp-chr-coords', multi=False),
        ], width=4, style={'margin-top': '0px', 'margin-bottom': '30px'})
    ]),

    dbc.Row([
        dbc.Col([
            html.Label("WT reference (if any) :"),
            dcc.Dropdown(id='pp-reference-selector', options=[], value=None, multi=False),
        ], width=4, style={'margin-top': '0px', 'margin-bottom': '30px'}),

        dbc.Col([
            html.Label("Additional groups of probes (if any) :"),
            dcc.Dropdown(id='pp-probe-groups', options=[], value=None, multi=False),
        ], width=4, style={'margin-top': '0px', 'margin-bottom': '30px'}),
    ]),

    dbc.Row([
        dbc.Col([
            html.Div(id='pp-p2f-dataframe-title',  style={'margin-top': '20px', 'margin-bottom': '20px'}),
            dcc.Loading(generate_data_table('pp-p2f-dataframe', [], [], 10))
        ], width=4, style={'margin-top': '0px', 'margin-bottom': '30px'}),

        dbc.Col([
            html.Div(id='pp-groups-dataframe-title', style={'margin-top': '20px', 'margin-bottom': '20px'}),
            dcc.Loading(generate_data_table('pp-groups-dataframe', [], [], 10))
        ], width=7, style={'margin-top': '0px', 'margin-bottom': '30px', 'margin-left': '30px'}),
    ]),

    dbc.Row([
        dbc.Col([
            html.Button(id="pp-copy-inputs", className="green-button", children="Copy files"),
            dbc.Tooltip(
                "Once you have selected all the inputs file you selected, you can copy them "
                "into the sample output directory to keep trace.",
                target="pp-copy-inputs", className="custom-tooltip", placement="right"),
        ], width=3, style={'margin-top': '0px', 'margin-bottom': '10px'}),
    ]),

    dbc.Row([
        dbc.Col([
            html.Button(id="pp-p2f", className="blue-button", children="Probes to fragments"),
            dbc.Tooltip(
                "Create a column in the oligo table with the corresponding fragment",
                target="pp-p2f", className="custom-tooltip", placement="right"),
        ], width=3, style={'margin-top': '0px', 'margin-bottom': '10px'}),
    ]),

    dbc.Row([
        dbc.Col([
            html.Div(id='pp-p2f-output', style={'margin-top': '10px', 'margin-bottom': '10px'}),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '10px'})
    ]),

    dbc.Row([
        dbc.Col([
            html.Button(id="pp-filter", className="blue-button", children="Filter"),
            dbc.Tooltip(
                "This module filters the contacts by removing contacts "
                "that do not concern digested fragments containing oligos",
                target="pp-filter", className="custom-tooltip", placement="right"),
        ], width=2, style={'margin-top': '0px', 'margin-bottom': '10px'}),
    ]),
    dbc.Row([
        dbc.Col([
            html.Div(id='pp-filter-output', style={'margin-top': '10px', 'margin-bottom': '10px'}),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '10px'})
    ]),

    dbc.Row([
        dbc.Col([
            html.Button(
                id="pp-coverage", className="blue-button", children="Coverage"),
            dbc.Tooltip(
                "Calculate the coverage per oligo fragment and save the result as a bed-graph file",
                target="pp-coverage", className="custom-tooltip", placement="right"),
        ], width=2, style={'margin-top': '0px', 'margin-bottom': '10px'}),
    ]),
    dbc.Row([
        dbc.Col([
            html.Div(id='pp-coverage-output', style={'margin-top': '10px', 'margin-bottom': '10px'}),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '10px'})
    ]),

    dbc.Row([
        dbc.Col([
            html.Button(id="pp-orga-contacts", className="blue-button", children="Organize contacts"),
            dbc.Tooltip(
                "Organize the contacts made by each probe with the genome and save "
                "the results as two .tsv files one for contacts and one for frequencies.",
                target="pp-orga-contacts", className="custom-tooltip", placement="right"),
        ], width=2, style={'margin-top': '0px', 'margin-bottom': '10px'}),
    ]),
    dbc.Row([
        dbc.Col([
            html.Div(id='pp-orga-contacts-output', style={'margin-top': '10px', 'margin-bottom': '10px'}),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '10px'})
    ])
])


@callback(
    Output('pp-sample-id-output', 'children'),
    Input('this-sample-id', 'data')
)
def display_sample_id(sample_id):
    return f"You are working on sample : {sample_id}"


@callback(
    Output('pp-copy-inputs', 'n_clicks'),
    [Input('pp-copy-inputs', 'n_clicks')],
    [State('this-sample-path', 'data'),
     State('pp-fragments-selector', 'value'),
     State('pp-oligo-selector', 'value'),
     State('pp-chr-coords', 'value'),
     State('pp-reference-selector', 'value'),
     State('this-sample-out-dir-path', 'data')]
)
def copy_input_files(
        n_clicks,
        sample_matrix,
        fragments_file, oligo_file,
        chr_coords_file,
        reference_file,
        sample_out_dir
):
    if n_clicks is None or n_clicks == 0:
        return None

    if n_clicks > 0:
        if sample_out_dir is None:
            return None

        files_to_copy = []
        inputs_dir = join(sample_out_dir, "inputs")
        if not os.path.exists(inputs_dir):
            os.makedirs(inputs_dir)
        if fragments_file is not None:
            files_to_copy.append(fragments_file)
        if oligo_file is not None:
            files_to_copy.append(oligo_file)
        if chr_coords_file is not None:
            files_to_copy.append(chr_coords_file)
        if sample_matrix is not None:
            files_to_copy.append(sample_matrix)

        for file in files_to_copy:
            copyfile(file, join(inputs_dir, basename(file)))

        if reference_file is not None:
            reference_file_name = basename(reference_file)
            reference_dir = join(inputs_dir, "references")
            if not os.path.exists(reference_dir):
                os.makedirs(reference_dir)
            copyfile(reference_file, join(reference_dir, reference_file_name))
        return 0


@callback(
    [Output('pp-fragments-selector', 'options'),
     Output('pp-oligo-selector', 'options'),
     Output('pp-chr-coords', 'options'),
     Output('pp-probe-groups', 'options'),
     Output('pp-reference-selector', 'options')],
    Input('data-basedir', 'data')
)
def update_dropdowns(data_basedir):
    if data_basedir is None:
        return [], [], [], [], []
    inputs_dir = join(data_basedir, "inputs")
    inputs_files = sorted([f for f in os.listdir(inputs_dir) if isfile(join(inputs_dir, f))],
                          key=lambda x: x.lower())
    options = [{'label': f, 'value': join(inputs_dir, f)} for f in inputs_files]

    reference_dir = join(inputs_dir, "references")
    references = sorted([f for f in os.listdir(reference_dir) if isfile(join(reference_dir, f))])
    ref_options = [{'label': f, 'value': join(reference_dir, f)} for f in references]
    return options, options, options, options, ref_options


def prepare_dataframe_for_output(dataframe, selected_columns=None):
    if dataframe is None:
        return None, None
    if selected_columns:
        df_output = dataframe[selected_columns]
        data = df_output.to_dict('records')
        columns = [{"name": col, "id": col} for col in selected_columns]
    else:
        data = dataframe.to_dict('records')
        columns = [{"name": col, "id": col} for col in dataframe.columns]
    return data, columns


@callback(
    [Output('pp-p2f', 'n_clicks'),
     Output('pp-p2f-output', 'children'),
     Output('pp-p2f-dataframe-title', 'children'),
     Output('pp-p2f-dataframe', 'data'),
     Output('pp-p2f-dataframe', 'columns')],
    [Input('pp-p2f', 'n_clicks')],
    [State('pp-fragments-selector', 'value'),
     State('pp-oligo-selector', 'value')]
)
def oligo_and_fragments(n_clicks, fragments_file, oligo_file):
    if n_clicks is None or n_clicks == 0:
        return 0, dash.no_update, dash.no_update, dash.no_update, dash.no_update

    if fragments_file is None or oligo_file is None:
        return 0, "Select both fragments and capture oligos files", dash.no_update, dash.no_update, dash.no_update

    df_oli = pd.read_csv(oligo_file, sep=utils.detect_delimiter(oligo_file))
    df_frag = pd.read_csv(fragments_file, sep=utils.detect_delimiter(fragments_file))

    title = html.H6("Oligo probes VS. Fragments ID:")
    if 'fragment' in df_oli.columns:
        data, columns = prepare_dataframe_for_output(df_oli, ["name", "fragment"])
        return 0, "Capture oligos table already contains a fragment column", title, data, columns
    else:
        probe2fragment.associate_probes_to_fragments(fragments_file, oligo_file)
        df_p2f = pd.read_csv(oligo_file, sep=utils.detect_delimiter(oligo_file))
        data, columns = prepare_dataframe_for_output(df_p2f, ["name", "fragment"])
        return 0, "Associated probes to fragments successfully", title, data, columns


@callback(
    [Output('pp-groups-dataframe-title', 'children'),
     Output('pp-groups-dataframe', 'data'),
     Output('pp-groups-dataframe', 'columns')],
    [Input('pp-probe-groups', 'value')]
)
def probe_groups(groups_file):
    if groups_file is None:
        return dash.no_update, dash.no_update, dash.no_update

    df_groups = pd.read_csv(groups_file, sep='\t')
    data, columns = prepare_dataframe_for_output(df_groups)
    title = html.H6("Probe groups:")
    return title, data, columns


@callback(
    [Output('pp-filter', 'n_clicks'),
     Output('pp-filter-output', 'children')],
    [Input('pp-filter', 'n_clicks')],
    [State('this-sample-out-dir-path', 'data'),
     State('this-sample-path', 'data'),
     State('pp-fragments-selector', 'value'),
     State('pp-oligo-selector', 'value')]
)
def filter_contacts(n_clicks, output_dir, sparse_matrix, fragments_file, oligos_file):
    if n_clicks is None or n_clicks == 0:
        return 0, dash.no_update

    pattern = re.compile(r'.+_filtered\.tsv')
    if n_clicks == 1:
        if output_dir is None:
            return dash.no_update, "You have to select a sample first"
        for file in os.listdir(output_dir):
            if pattern.match(file):
                return n_clicks, "Filtered contacts file already exists (click again to overwrite)"

    if fragments_file is None:
        return 0, "Select a digested fragments file"
    if oligos_file is None:
        return 0, "Select a capture oligos file"

    filter.filter_contacts(oligos_file, fragments_file, sparse_matrix, output_dir)
    return 0, "Filtered contacts file created successfully"


@callback(
    [Output('pp-coverage', 'n_clicks'),
     Output('pp-coverage-output', 'children')],
    [Input('pp-coverage', 'n_clicks')],
    [State('this-sample-out-dir-path', 'data'),
     State('this-sample-path', 'data'),
     State('pp-fragments-selector', 'value')]
)
def compute_cover(n_clicks, output_dir, sparse_matrix, fragments_file):
    if n_clicks is None or n_clicks == 0:
        return 0, dash.no_update

    pattern = re.compile(r'.+_coverage_')
    if n_clicks == 1:
        if output_dir is None:
            return dash.no_update, "You have to select a sample first"
        for file in os.listdir(output_dir):
            if pattern.match(file):
                return n_clicks, "Coverage bed-graph file already exists (click again to overwrite)"

    if fragments_file is None:
        return 0, "Select a digested fragments file"

    coverage.coverage(sparse_matrix, fragments_file, output_dir)
    return 0, "Coverage file created successfully"

