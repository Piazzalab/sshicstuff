import os
import re
import dash
import json

import pandas as pd
from os.path import basename, join, isfile, isdir
from shutil import copyfile
import dash_bootstrap_components as dbc
from dash import callback
from dash import html, dcc, dash_table
from dash.dependencies import Input, Output, State, ALL

import core.utils
import core.filter
import core.coverage
import core.fragments
import core.probe2fragment
import core.binning
import core.statistics
import core.weight
import core.aggregated
import utils

from common import generate_data_table, prepare_dataframe_for_output

layout = dbc.Container([
    dbc.Row([
        dbc.Col(
            dbc.Card([
                dbc.CardHeader("Samples"),
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            html.Label("Select a sample : "),
                            dcc.Dropdown(id='pp-sample-selector', options=[], value=None, multi=False),
                            html.Div(id='pp-current-sample-id-output',
                                     style={'margin-top': '20px', 'margin-bottom': '20px'})
                        ], width=4, style={'margin-top': '0px', 'margin-bottom': '0px'}),

                        dbc.Col([
                            html.Div(id='pp-current-sample-files',
                                     style={'margin-top': '0px', 'margin-bottom': '20px'}),
                            html.Div(id='pp-current-sample-pcr-output',
                                     style={'margin-top': '20px', 'margin-bottom': '20px'})
                        ], width=8),
                    ]),
                    dcc.Store(id='pp-current-sample-id'),
                    dcc.Store(id='pp-current-sample-file-path'),
                    dcc.Store(id='pp-current-sample-in-dir-path'),
                    dcc.Store(id='pp-current-sample-out-dir-path')
                ])
            ])
        )
    ], style={'margin-top': '20px', 'margin-bottom': '50px'}),

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
            html.Button(id="pp-copy-inputs-button", className="green-button", children="Copy files"),
            dbc.Tooltip(
                "Once you have selected all the inputs file you selected, you can copy them "
                "into the sample output directory to keep trace.",
                target="pp-copy-inputs-button", className="custom-tooltip", placement="right"),
        ], width=3, style={'margin-top': '0px', 'margin-bottom': '10px'}),
    ]),

    dbc.Row([
        dbc.Col(
            dbc.Card([
                dbc.CardHeader("Prepare Data"),
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            html.Button(id="pp-p2f-button", className="blue-button", children="Probes to frags"),
                            dbc.Tooltip(
                                "Create a column in the oligo table with the corresponding fragment",
                                target="pp-p2f-button", className="custom-tooltip", placement="right"),
                        ], width=3, style={'margin-top': '0px', 'margin-bottom': '0px'}),

                        dbc.Col([
                            html.Button(id="pp-filter-button", className="blue-button", children="Filter"),
                            dbc.Tooltip(
                                "This module filters the contacts by removing contacts "
                                "that do not concern digested fragments containing oligos",
                                target="pp-filter-button", className="custom-tooltip", placement="right"),
                        ], width=3, style={'margin-top': '0px', 'margin-bottom': '0px'}),

                        dbc.Col([
                            html.Button(
                                id="pp-coverage-button", className="blue-button", children="Coverage"),
                            dbc.Tooltip(
                                "Calculate the coverage per oligo fragment and "
                                "save the result as a bed-graph file",
                                target="pp-coverage-button", className="custom-tooltip", placement="right"),
                        ], width=3, style={'margin-top': '0px', 'margin-bottom': '0px'}),
                    ]),

                    dbc.Row([
                        dbc.Col([
                            html.Div(id='pp-p2f-output', style={'margin-top': '10px', 'margin-bottom': '10px'}),
                        ], width=3, style={'margin-top': '0px', 'margin-bottom': '10px'}),

                        dbc.Col([
                            html.Div(id='pp-filter-output', style={'margin-top': '10px', 'margin-bottom': '10px'}),
                        ], width=3, style={'margin-top': '0px', 'margin-bottom': '10px'}),

                        dbc.Col([
                            html.Div(id='pp-coverage-output', style={'margin-top': '10px', 'margin-bottom': '10px'}),
                        ], width=3, style={'margin-top': '0px', 'margin-bottom': '10px'})

                    ])
                ])
            ])
        )
    ], style={'margin-top': '0px', 'margin-bottom': '50px'}),

    dbc.Row([
        dbc.Col(
            dbc.Card([
                dbc.CardHeader("Contacts Resolution"),
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            html.Button(id="pp-orga-contacts-button",
                                        className="blue-button", children="Organize contacts"),
                            dbc.Tooltip(
                                "Organize the contacts made by each probe with the genome and save "
                                "the results as two .tsv files one for contacts and one for frequencies.",
                                target="pp-orga-contacts-button", className="custom-tooltip", placement="right"),
                        ], width=3, style={'margin-top': '0px', 'margin-bottom': '10px'}),

                        dbc.Col([
                            html.Button(id="pp-binning-button", className="blue-button", children="Binning"),
                            dbc.Tooltip(
                                "Change the resolution of contacts tables (1kb, 5kb, 10kb etc ...)",
                                target="pp-binning-button", className="custom-tooltip", placement="right"),
                        ], width=2, style={'margin-top': '0px', 'margin-bottom': '10px'}),

                        dbc.Col([
                            dcc.Input(id='pp-binning-input-box', type='number', value='', step='1',
                                      placeholder='Specify binning (in kb)',
                                      style={
                                          'width': '100%',
                                          'border': '1px solid #ccc',
                                          'border-radius': '4px',
                                          'padding': '10px',
                                          'font-size': '16px',
                                          'background-color': '#fff',
                                          'color': '#333'
                                      }),
                            dcc.Store(id='pp-stored-bins', data=[])
                        ], width=3, style={'margin-top': '0px', 'margin-bottom': '0px', 'margin-left': '0px'}),

                        dbc.Col([
                            html.Div(id='pp-bins-list', style={'display': 'flex', 'flexWrap': 'wrap'})
                        ], width=3, style={'margin-top': '0px', 'margin-bottom': '0px', 'margin-left': '0px'}),

                    ]),

                    dbc.Row([
                        dbc.Col([
                            html.Div(id='pp-orga-contacts-output',
                                     style={'margin-top': '10px', 'margin-bottom': '10px'}),
                        ], width=3, style={'margin-top': '0px', 'margin-bottom': '10px'}),

                        dbc.Col([
                            html.Div(id='pp-binning-output', style={'margin-top': '10px', 'margin-bottom': '10px'}),
                        ], width=3, style={'margin-top': '0px', 'margin-bottom': '10px'})
                    ]),
                ])
            ])
        )
    ], style={'margin-top': '0px', 'margin-bottom': '50px'}),

    dbc.Row([
        dbc.Col(
            dbc.Card([
                dbc.CardHeader("Capture Efficiency"),
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            html.Button(id="pp-stats-button", className="blue-button", children="Statistics"),
                            dbc.Tooltip("Generate statistics and normalization for contacts made by each probe",
                                        target="pp-stats-button", className="custom-tooltip", placement="right"),
                        ], width=2, style={'margin-top': '0px', 'margin-bottom': '10px'}),

                        dbc.Col([
                            dcc.Input(id='pp-stats-cis-range-input-box', type='number', value="", step='1',
                                      placeholder='Specify a cis range (in bp)',
                                      style={
                                          'width': '80%',
                                          'border': '1px solid #ccc',
                                          'border-radius': '4px',
                                          'padding': '10px',
                                          'font-size': '16px',
                                          'background-color': '#fff',
                                          'color': '#333'
                                      }),
                            dbc.Tooltip("Range of bp around the probes (both left and right) "
                                        "to consider as cis contacts",
                                        target="pp-stats-cis-range-input-box",
                                        className="custom-tooltip", placement="right"),
                        ], width=4, style={'margin-top': '0px', 'margin-bottom': '0px'}),

                        dbc.Col([
                            html.Button(id="pp-weight-button", className="blue-button", children="Weight"),
                            dbc.Tooltip("Weight a sample by normalizing its contacts"
                                        " by a coefficient of capture"
                                        "efficiency compared to a reference WT",
                                        target="pp-weight-button", className="custom-tooltip", placement="right"),
                        ], width=3, style={'margin-top': '0px', 'margin-bottom': '10px'}),
                    ]),

                    dbc.Row([
                        dbc.Col([
                            html.Div(id='pp-stats-output', style={'margin-top': '10px', 'margin-bottom': '10px'}),
                        ], width=6, style={'margin-top': '0px', 'margin-bottom': '10px'}),

                        dbc.Col([
                            html.Div(id='pp-weight-output', style={'margin-top': '10px', 'margin-bottom': '10px'}),
                        ], width=3, style={'margin-top': '0px', 'margin-bottom': '10px'})
                    ]),
                ])
            ])
        )
    ], style={'margin-top': '0px', 'margin-bottom': '50px'}),



    dbc.Row([
        dbc.Col(
            dbc.Card([
                dbc.CardHeader("Aggregated"),
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            html.Button(id="pp-aggregate-button", className="blue-button", children="Aggregate"),
                            dbc.Tooltip("Aggregate contacts made by probes around centromeres and/or telomeres",
                                        target="pp-aggregate-button", className="custom-tooltip", placement="right"),
                        ], width=2, style={'margin-top': '0px', 'margin-bottom': '10px'}),
                    ]),

                    dbc.Row([
                        dbc.Col([
                            dcc.Dropdown(id='pp-aggr-weight-selector', options=[], value=None, multi=False),
                        ], width=2, style={'margin-top': '0px', 'margin-bottom': '30px'}),

                        dbc.Col([
                            dcc.Dropdown(id='pp-aggr-on-selector', options=["centromeres", "telomeres"], value=None,
                                         multi=False),
                        ], width=2, style={'margin-top': '0px', 'margin-bottom': '30px'}),

                        dbc.Col([
                            dcc.Input(id='pp-aggr-window', type='number', value="", step='1',
                                      placeholder="Specify window region (in bp)",
                                      style={
                                          'width': '100%',
                                          'height': '36px',
                                          'border': '1px solid #ccc',
                                          'border-radius': '4px',
                                          'padding': '10px',
                                          'font-size': '14px',
                                          'background-color': '#fff',
                                          'color': '#333'
                                      }),
                            dbc.Tooltip("Window (in bp) that defines the centromere or "
                                        "telomere region (on 5' and 3')",
                                        target='pp-aggr-window', className="custom-tooltip", placement="right"),
                        ], width=3, style={'margin-top': '0px', 'margin-bottom': '0px', 'margin-left': '0px'}),

                        dbc.Col([
                            dcc.Dropdown(id='pp-aggr-chr-exclusion-selector', options=[], value=None, multi=True),
                            dbc.Tooltip("Enter the chromosome you want to exclude from aggregation.",
                                        target='pp-chr-exclusion-selector', className="custom-tooltip",
                                        placement="right"),
                        ], width=2, style={'margin-top': '0px', 'margin-bottom': '0px'}),

                    ]),

                    dbc.Row([
                        dbc.Col([
                            dcc.Checklist(id='pp-aggr-self-chr-checkbox',
                                          options=[{'label': ' Exclude probe located chr', 'value': 'checked'}],
                                          value=[]),
                        ], width=3, style={'margin-top': '0px', 'margin-bottom': '30px'}),

                        dbc.Col([
                            dcc.Checklist(id='pp-aggr-inter-norm-checkbox',
                                          options=[{'label': ' Inter-norm', 'value': 'checked'}], value=[]),
                        ], width=2, style={'margin-top': '0px', 'margin-bottom': '30px', 'margin-left': '-70px'}),

                        dbc.Col([
                            dcc.Checklist(id='pp-aggr-plot-checkbox',
                                          options=[{'label': ' Plot', 'value': 'checked'}], value=[]),
                        ], width=2, style={'margin-top': '0px', 'margin-bottom': '0px', 'margin-left': '-70px'}),
                    ]),

                    dbc.Row([
                        dbc.Col([
                            html.Div(id='pp-aggregate-output', style={'margin-top': '10px', 'margin-bottom': '10px'}),
                        ], width=6, style={'margin-top': '0px', 'margin-bottom': '10px'})
                    ]),
                ])
            ], style={'height': '240px'})
        )
    ])
])


@callback(
    Output('pp-sample-selector', 'options'),
    Input('selected-samples', 'data'),
)
def update_sample_selector(selected_samples):
    if selected_samples is None:
        return []
    options = [{'label': s, 'value': s} for s in selected_samples]
    return options


@callback(
    Output('pp-current-sample-id-output', 'children'),
    Output('pp-current-sample-id', 'data'),
    Input('pp-sample-selector', 'value')
)
def update_current_sample_id_output(sample_id):
    if sample_id is None:
        return None, None
    return f"Current sample ID: {sample_id}", sample_id


@callback(
    Output('pp-current-sample-files', 'children'),
    Input('pp-sample-selector', 'value'),
    State('data-basedir', 'data')
)
def display_samples_files(sample_id, data_basedir):
    if sample_id is None:
        return None

    samples_dir = join(data_basedir, "samples")
    current_samp_files = [
        f for f in os.listdir(samples_dir) if isfile(join(samples_dir, f)) and sample_id.lower() in f.lower()
    ]
    return html.Div([
        html.Label("Select a file : "),
        dcc.Dropdown(
            id='pp-current-sample-file-selector',
            options=current_samp_files, value=None, multi=False),
    ])


@callback(
    Output('pp-current-sample-file-path', 'data'),
    Output('pp-current-sample-in-dir-path', 'data'),
    Output('pp-current-sample-out-dir-path', 'data'),
    Output('pp-current-sample-pcr-output', 'children'),
    Input('pp-current-sample-file-selector', 'value'),
    State('data-basedir', 'data'),
    State('pp-current-sample-id', 'data')
)
def update_current_sample_paths(sample_file, data_basedir, sample_id):
    if sample_file is None or data_basedir is None or sample_id is None:
        return None, None, None, None

    sample_file_path = join(data_basedir, "samples", sample_file)

    sample_dir = join(data_basedir, "outputs", sample_id)
    sample_input_dir = join(sample_dir, "inputs")
    pcr_output = None
    if "pcrfree" in sample_file.lower():
        sample_output_dir = join(sample_dir, "pcrfree")
        pcr_output = "pcrfree"
    elif "pcrdupkept" in sample_file.lower():
        sample_output_dir = join(sample_dir, "pcrdupkept")
        pcr_output = "pcrdupkept"
    else:
        sample_output_dir = sample_dir

    if not os.path.exists(sample_output_dir):
        os.makedirs(sample_output_dir)
    if not os.path.exists(sample_input_dir):
        os.makedirs(sample_input_dir)
    return sample_file_path, sample_input_dir, sample_output_dir, pcr_output


@callback(
    Output('pp-copy-inputs-button', 'n_clicks'),
    [Input('pp-copy-inputs-button', 'n_clicks')],
    [State('pp-current-sample-file-path', 'data'),
     State('pp-fragments-selector', 'value'),
     State('pp-oligo-selector', 'value'),
     State('pp-chr-coords', 'value'),
     State('pp-reference-selector', 'value'),
     State('pp-current-sample-out-dir-path', 'data'),
     State('pp-current-sample-in-dir-path', 'data')]
)
def copy_input_files(
        n_clicks,
        sample_matrix,
        fragments_file,
        oligo_file,
        chr_coords_file,
        reference_file,
        sample_out_dir,
        sample_in_dir
):
    if n_clicks is None or n_clicks == 0:
        return None

    if n_clicks > 0:
        if sample_out_dir is None:
            return None

        files_to_copy = []
        if fragments_file is not None:
            files_to_copy.append(fragments_file)
        if oligo_file is not None:
            files_to_copy.append(oligo_file)
        if chr_coords_file is not None:
            files_to_copy.append(chr_coords_file)
        if sample_matrix is not None:
            files_to_copy.append(sample_matrix)

        for file in files_to_copy:
            copyfile(file, join(sample_in_dir, basename(file)))

        if reference_file is not None:
            reference_file_name = basename(reference_file)
            copyfile(reference_file, join(sample_in_dir, reference_file_name))
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


@callback(
    [Output('pp-p2f-button', 'n_clicks'),
     Output('pp-p2f-output', 'children'),
     Output('pp-p2f-dataframe-title', 'children'),
     Output('pp-p2f-dataframe', 'data'),
     Output('pp-p2f-dataframe', 'columns')],
    [Input('pp-p2f-button', 'n_clicks')],
    [Input('pp-fragments-selector', 'value'),
     Input('pp-oligo-selector', 'value')]
)
def oligo_and_fragments(n_clicks, fragments_file, oligo_file):
    if n_clicks is None or n_clicks == 0:
        return 0, None, None, [], []

    if fragments_file is None or oligo_file is None:
        return 0, "Select both fragments and capture oligos files", None, [], []

    df_oli = pd.read_csv(oligo_file, sep=core.utils.detect_delimiter(oligo_file))
    df_frag = pd.read_csv(fragments_file, sep=core.utils.detect_delimiter(fragments_file))

    title = html.H6("Oligo probes VS. Fragments ID:")
    if 'fragment' in df_oli.columns:
        data, columns = prepare_dataframe_for_output(df_oli, ["name", "fragment"])
        return 0, "Capture oligos table already contains a fragment column", title, data, columns
    else:
        core.probe2fragment.associate_probes_to_fragments(fragments_file, oligo_file)
        df_p2f = pd.read_csv(oligo_file, sep=core.utils.detect_delimiter(oligo_file))
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
        return None, [], []

    df_groups = pd.read_csv(groups_file, sep='\t')
    data, columns = prepare_dataframe_for_output(df_groups)
    title = html.H6("Probe groups :")
    return title, data, columns


@callback(
    [Output('pp-filter-button', 'n_clicks'),
     Output('pp-filter-output', 'children')],
    [Input('pp-filter-button', 'n_clicks')],
    [State('pp-current-sample-out-dir-path', 'data'),
     State('pp-current-sample-id', 'data'),
     State('pp-current-sample-file-path', 'data'),
     State('pp-fragments-selector', 'value'),
     State('pp-oligo-selector', 'value')]
)
def filter_contacts(n_clicks, output_dir, sample_id, sparse_matrix, fragments_file, oligos_file):

    if n_clicks is None or n_clicks == 0:
        return 0, dash.no_update

    if sample_id is None:
        return 0, "You have to select a sample first"

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

    core.filter.filter_contacts(oligos_file, fragments_file, sparse_matrix, output_dir)
    return 0, "Filtered contacts file created successfully"


@callback(
    [Output('pp-coverage-button', 'n_clicks'),
     Output('pp-coverage-output', 'children')],
    [Input('pp-coverage-button', 'n_clicks')],
    [State('pp-current-sample-out-dir-path', 'data'),
     State('pp-current-sample-file-path', 'data'),
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

    core.coverage.coverage(sparse_matrix, fragments_file, output_dir)
    return 0, "Coverage file created successfully"


@callback(
    [Output('pp-orga-contacts-button', 'n_clicks'),
     Output('pp-orga-contacts-output', 'children')],
    [Input('pp-orga-contacts-button', 'n_clicks'),
     State('pp-current-sample-out-dir-path', 'data'),
     State('pp-current-sample-id', 'data'),
     State('pp-oligo-selector', 'value'),
     State('pp-chr-coords', 'value'),
     State('pp-probe-groups', 'value')]
)
def fragment_contacts(n_clicks, sample_output_dir, sample_id, oligos_file, chr_coords, groups_file):
    if n_clicks is None or n_clicks == 0:
        return 0, dash.no_update

    if sample_id is None:
        return 0, "You need to select a sample first"

    output_dir = join(sample_output_dir, 'not_weighted')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    filtered_sample = join(sample_output_dir, f"{sample_id}_filtered.tsv")
    if filtered_sample is None:
        return 0, "You need to filter the sample first"
    if oligos_file is None:
        return 0, "Select a capture oligos file"
    if chr_coords is None:
        return 0, "Select a chromosome coordinates file"

    pattern = re.compile(r'.+_unbinned_contacts')
    if n_clicks == 1:
        for file in os.listdir(output_dir):
            if pattern.match(file):
                return n_clicks, "Contacts file already exists (click again to overwrite)"

    core.fragments.organize_contacts(filtered_sample, oligos_file, chr_coords, output_dir, groups_file)
    return 0, "Contacts file created successfully"


@callback(
    [Output('pp-bins-list', 'children'),
     Output('pp-binning-input-box', 'value'),
     Output('pp-stored-bins', 'data')],
    Input({'type': 'delete-button', 'index': ALL}, 'n_clicks'),
    Input('pp-binning-input-box', 'n_submit'),
    State('pp-binning-input-box', 'value'),
    State('pp-stored-bins', 'data')
)
def update_bins_list(delete_n_clicks_list, n_submit, input_value, stored_numbers):
    ctx = dash.callback_context

    if not ctx.triggered:
        return dash.no_update

    triggered_id = ctx.triggered[0]['prop_id']

    if triggered_id == 'pp-binning-input-box.n_submit' and input_value:
        if int(input_value) in stored_numbers:
            return dash.no_update
        stored_numbers.append(int(input_value))

    if 'delete-button' in triggered_id:
        clicked_button = int(json.loads(triggered_id.split('.n_clicks')[0])["index"])
        stored_numbers.pop(clicked_button)

    numbers_list = [
        html.Div([
            f"{number} kb ",
            html.Button("❌", className="custom-delete-button",
                        id={'type': 'delete-button', 'index': index})
        ], style={'display': 'inline-block', 'margin': '0px', 'width': '100px', 'height': '40px'})
        for index, number in enumerate(stored_numbers)
    ]

    return numbers_list, '', stored_numbers


@callback(
    [Output('pp-binning-button', 'n_clicks'),
     Output('pp-binning-output', 'children')],
    [Input('pp-binning-button', 'n_clicks'),
     State('pp-stored-bins', 'data'),
     State('pp-current-sample-out-dir-path', 'data'),
     State('pp-current-sample-id', 'data')],
    [State('pp-oligo-selector', 'value'),
     State('pp-chr-coords', 'value'),
     State('pp-probe-groups', 'value')]
)
def make_rebin(n_clicks, bins_list, sample_output_dir, sample_id, oligos_file, chr_coords, groups_file):

    if n_clicks is None or n_clicks == 0:
        return 0, dash.no_update
    if sample_id is None:
        return 0, "You need to select a sample first"
    if oligos_file is None:
        return 0, "Select a capture oligos file"
    if chr_coords is None:
        return 0, "Select a chromosome coordinates file"
    if len(bins_list) == 0:
        return 0, "Select at least one bin size (in kb)"

    unbinned_contacts = join(sample_output_dir, 'not_weighted', f"{sample_id}_unbinned_contacts.tsv")
    if not os.path.exists(unbinned_contacts):
        return 0, "You need to create fragment contacts tables (unbinned) first"
    output_dir = join(sample_output_dir, 'not_weighted')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if n_clicks == 1:
        for bin_kb in bins_list:
            bin_bp = int(bin_kb) * 1000
            core.binning.rebin_contacts(unbinned_contacts, chr_coords, oligos_file, bin_bp, output_dir, groups_file)

    return 0, "Binned contacts files created successfully"


@callback(
    [Output('pp-stats-button', 'n_clicks'),
     Output('pp-stats-output', 'children')],
    [Input('pp-stats-button', 'n_clicks')],
    [State('pp-current-sample-out-dir-path', 'data'),
     State('pp-current-sample-id', 'data'),
     State('pp-current-sample-file-path', 'data'),
     State('pp-oligo-selector', 'value'),
     State('pp-reference-selector', 'value'),
     State('pp-stats-cis-range-input-box', 'value')]
)
def make_statistics(n_clicks, sample_output_dir, sample_id, sample_path, oligos_file, reference, cis_range):
    if n_clicks is None or n_clicks == 0:
        return 0, dash.no_update
    if sample_id is None or sample_path is None:
        return 0, "You need to select a sample first"
    if oligos_file is None:
        return 0, "Select a capture oligos file"
    if cis_range < 0:
        return 0, "Cis range must be positive integer"

    output_dir = sample_output_dir
    sparse_matrix = sample_path
    unbinned_contacts = join(output_dir, 'not_weighted', f"{sample_id}_unbinned_contacts.tsv")

    global_stats = join(output_dir, f"{sample_id}_global_statistics.tsv")

    if n_clicks == 1:
        if global_stats in os.listdir(output_dir):
            return n_clicks, "Statistics file already exists (click again to overwrite)"

    core.statistics.get_stats(unbinned_contacts, sparse_matrix, oligos_file, output_dir, cis_range)
    if reference is not None:
        ref_name = reference.split('/')[-1].split('.')[0]
        core.statistics.compare_to_wt(global_stats, reference, ref_name)

        return 0, "Statistics files created successfully and compared to wild type reference"
    return 0, "Statistics files created successfully"


@callback(
    [Output('pp-weight-button', 'n_clicks'),
     Output('pp-weight-output', 'children')],
    [Input('pp-weight-button', 'n_clicks')],
    [State('pp-current-sample-out-dir-path', 'data'),
     State('pp-current-sample-id', 'data'),
     State('pp-reference-selector', 'value'),
     State('pp-probe-groups', 'value')]
)
def make_statistics(n_clicks, sample_output_dir, sample_id, reference, groups_file):
    if n_clicks is None or n_clicks == 0:
        return 0, dash.no_update
    if sample_id is None:
        return 0, "You need to select a sample first"
    if reference is None:
        return 0, "You need first to select a reference WT to weight with"

    global_stats = join(sample_output_dir, f"{sample_id}_global_statistics.tsv")
    ref_name = reference.split('/')[-1].split('.')[0]

    not_weighted_dir = join(sample_output_dir, 'not_weighted')
    weighted_dir = join(sample_output_dir, f'weighted_{ref_name}')
    binned_contacts_list = [f for f in os.listdir(not_weighted_dir) if '_binned_contacts' in f]
    binned_frequencies_list = [f for f in os.listdir(not_weighted_dir) if '_binned_frequencies' in f]
    unbinned_contacts = join(not_weighted_dir, f"{sample_id}_unbinned_contacts.tsv")
    unbinned_frequencies = join(not_weighted_dir, f"{sample_id}_unbinned_frequencies.tsv")

    if n_clicks == 1:
        if not os.path.exists(weighted_dir):
            os.makedirs(weighted_dir)

        core.weight.weight_mutant(
            global_stats, ref_name, unbinned_contacts, unbinned_frequencies,
            "unbinned", weighted_dir, groups_file)

        for binned_c, binned_f in zip(binned_contacts_list, binned_frequencies_list):
            bin_suffix = re.search(r'(\d+)kb', binned_c).group(1) + 'kb'
            core.weight.weight_mutant(
                global_stats, ref_name, join(not_weighted_dir, binned_c), join(not_weighted_dir, binned_f),
                f"{bin_suffix}_binned", weighted_dir, groups_file)

        return 0, "Weighted files created successfully"


@callback(
    [Output('pp-aggr-weight-selector', 'options')],
    [Input('pp-current-sample-out-dir-path', 'data')]
)
def update_aggr_weight_selector(sample_output_dir):
    if sample_output_dir is None:
        return dash.no_update

    weighted_dirs = [d for d in os.listdir(sample_output_dir) if 'weighted' in d and isdir(join(sample_output_dir, d))]
    options = [{'label': d, 'value': join(sample_output_dir, d)} for d in weighted_dirs]
    return [options]


@callback(
    [Output('pp-aggr-chr-exclusion-selector', 'options')],
    [Input('pp-chr-coords', 'value')]
)
def update_aggr_chr_exclusion_selector(chr_coords):
    if chr_coords is None:
        return dash.no_update

    df_chr_coords = pd.read_csv(chr_coords, sep=utils.detect_delimiter(chr_coords))
    chr_list = df_chr_coords['chr'].to_list()
    options = [{'label': c, 'value': c} for c in chr_list]
    return [options]

