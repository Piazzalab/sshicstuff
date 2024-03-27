import os
import base64

import pandas as pd
import numpy as np
import dash
import json
from os.path import join, dirname
from os import listdir
from dash import html
from dash import dcc
from dash import callback
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State, ALL, MATCH
import plotly.graph_objs as go

from common import generate_data_table, prepare_dataframe_for_output
import core.utils


TEMPORARY_DIRECTORY = join(dirname(dirname(os.getcwd())), "data", "__cache__")


colors = [
    'rgba(0, 0, 255, 0.8)',  # blue
    'rgba(255, 0, 0, 0.8)',  # red
    'rgba(249, 172, 37, 0.8)',  # yellow
    'rgba(245, 0, 87, 0.8)',  # pink
    'rgba(29, 233, 182, 0.8)',  # green
    'rgba(255, 234, 0, 0.8)',  # yellow 2
    'rgba(255, 11, 0, 0.8)',  # orange
    'rgba(141, 110, 99, 0.8)',  # brown
    'rgba(255, 64, 129, 0.8)',  # pink 2
    'rgba(120, 144, 156, 0.8)',  # blue grey
    'rgba(0, 131, 143, 0.8)',  # cyan
    'rgba(171, 71, 188, 0.8)',  # purple
    'rgba(255, 152, 0, 0.8)',  # amber
    'rgba(0, 150, 136, 0.8)',  # teal
    'rgba(0, 184, 212, 0.8)',  # cyan 2
    'rgba(0, 200, 83, 0.8)',  # green 2
    'rgba(229, 115, 115, 0.8)',  # red 2
    'rgba(255, 167, 38, 0.8)',  # orange 2
    'rgba(61, 90, 254, 0.8)',  # indigo
    'rgba(68, 138, 255, 0.8)',  # blue 2
    'rgba(121, 134, 203, 0.8)',  # deep purple
    'rgba(170, 102, 68, 0.8)',  # deep orange
    'rgba(255, 171, 145, 0.8)',  # pink 3
    'rgba(255, 209, 128, 0.8)'  # amber 2
]

chr_names = [f"chr{i}" for i in range(1, 17)] + ["2_micron", "mitochondrion", "chr_artificial"]
chr_pos = [230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751,
           666816, 1078177, 924431, 784333, 1091291, 948066, 6318, 85779, 7828]
chr_colors = ['#000000', '#0c090a', '#2c3e50', '#34495e', '#7f8c8d', '#8e44ad', '#2ecc71', '#2980b9',
              '#f1c40f', '#d35400', '#e74c3c', '#c0392b', '#1abc9c', '#16a085', '#bdc3c7', '#2c3e50',
              '#7f8c8d', '#f39c12', '#27ae60']


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
    dbc.Row([
        dbc.Col([
            html.H2('Probes Viewer'),
        ], width=12, style={'margin-top': '20px', 'margin-bottom': '20px'})
    ]),

    dbc.Row([
        dbc.Col([
            html.H6('Upload files'),
            dcc.Upload(
                id="pv-upload-files",
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
                id="pv-clear-list",
                className="btn btn-danger",
                children="Clear list",
            )
        ], width=2, style={'margin-top': '50px', 'margin-bottom': '20px', 'margin-left': '20px'}),
    ]),

    dbc.Row([
        dbc.Col([
            dcc.Dropdown(
                id='pv-samples-dropdown',
                options=[],
                value=[],
                placeholder="Select samples",
                multi=True
            )
        ], width=10, style={'margin-top': '10px', 'margin-bottom': '20px'})
    ]),

    dbc.Row([
        dbc.Col([
            dcc.Dropdown(id='pv-oligo-dropdown',
                         placeholder="Select capture oligos file",
                         multi=False),
        ], width=10, style={'margin-top': '10px', 'margin-bottom': '20px'}),
    ]),

    dbc.Row([
        html.Div(id='pv-probes-dropdown', children=[], style={'margin-top': '10px', 'margin-bottom': '20px'}),
    ]),

    dbc.Row([
        dbc.Col([
            html.Label("Select binning :", style={'margin-top': '10px', 'margin-bottom': '20px'}),
            dcc.Slider(
                id='pv-binning-slider',
                min=0,
                max=100,
                step=1,
                value=10,
                marks={i: str(i) for i in range(0, 101, 10)},
                included=False,
            ),
            html.Div(id='pv-slider-output-container', style={'margin-top': '20px', 'font-size': '16px'}),
            html.Br(),
        ], width=6, style={'margin-top': '0px', 'margin-bottom': '10px', 'margin-left': '0px'}),

        dbc.Col([
            dcc.Checklist(
                id="pv-sync-box",
                options=[{"label": "Sync axis", "value": "sync"}],
                value=[],
                inline=True,
                className='custom-checkbox-label',
                labelStyle={"margin": "5px"}
            )
        ], width=2, style={'margin-top': '15px', 'margin-bottom': '10px', 'margin-left': '0px'}),

        dbc.Col([
            dcc.Checklist(
                id="pv-same-graph",
                options=[{"label": "Single graph", "value": "single"}],
                value=[],
                inline=True,
                className='custom-checkbox-label',
                labelStyle={"margin": "5px"}
            )
        ], width=2, style={'margin-top': '15px', 'margin-bottom': '10px', 'margin-left': '-50px'}),

        dbc.Col([
            html.Button(id="pv-plot-buttom", className="plot-button", children="Plot"),
        ], width=2, style={'margin-top': '20px', 'margin-bottom': '0px', 'margin-left': '0px'}),

        dcc.Store(id='pv-stored-graphs-axis-range', data={}),
        html.Div(id='pv-graphs', children=[], style={'margin-top': '20px', 'margin-bottom': '20px'}),
    ]),
])


@callback(
    [Output("pv-samples-dropdown", "options"),
     Output("pv-oligo-dropdown", "options"),
     Output("pv-clear-list", "n_clicks")],
    [Input("pv-upload-files", "filename"),
     Input("pv-upload-files", "contents"),
     Input("pv-clear-list", "n_clicks")],
)
def update_file_list(uploaded_filenames, uploaded_file_contents, n_clicks):
    if uploaded_filenames is not None and uploaded_file_contents is not None:
        for name, data in zip(uploaded_filenames, uploaded_file_contents):
            save_file(name, data)

    files = uploaded_files()
    if n_clicks is not None:
        if n_clicks > 0:
            for filename in files:
                os.remove(os.path.join(TEMPORARY_DIRECTORY, filename))
            files = []

    n_clicks = 0
    if len(files) == 0:
        return files, files, n_clicks
    else:
        options = []
        for filename in files:
            options.append({'label': filename, 'value': os.path.join(TEMPORARY_DIRECTORY, filename)})
        return options, options,  n_clicks


def update_table(file_path, delim):
    if file_path and delim:
        df = pd.read_csv(file_path, sep=delim)
        data = df.to_dict('records')
        columns = [{"name": i, "id": i} for i in df.columns]
        return data, columns
    return None, None


@callback(
    [Output('pv-dataframe', 'data'),
     Output('pv-dataframe', 'columns')],
    [Input("pv-samples-dropdown", 'value')]
)
def update_dataframe(file_path, delim):
    if file_path is not None and delim is not None:
        return update_table(file_path, delim)
    return [], []


@callback(
    Output('pv-slider-output-container', 'children'),
    [Input('pv-binning-slider', 'value')])
def update_output(value):
    return f'You have selected a binning of {value} kb'


@callback(
    Output('pv-probes-dropdown', 'children'),
    Input('pv-oligo-dropdown', 'value')
)
def update_probes_dropdown(oligo_file):
    if not oligo_file:
        return None

    df = pd.read_csv(oligo_file, sep=',')
    probes = df['name'].to_list()
    return dbc.Col([
            dcc.Dropdown(
                id='pv-probes-dropdown',
                placeholder="Select probes to show",
                options=[{'label': p, 'value': p} for p in probes],
                value=None,
                multi=True
            )
        ], width=10, style={'margin-top': '10px', 'margin-bottom': '20px'}
    )


# def update_figure(
#         graph_id: int,
#         graph_dict: dict,
#         traces_colors: list,
#         binning: int,
#         chr_boundaries: list,
#         x_range=None,
#         y_range=None
# ):
#     fig = go.Figure()
#     trace_id = 0
#     for j in range(graph_dict['size']):
#         samp = graph_dict['samples'][j]
#         frag = graph_dict['fragments'][j]
#         pcr = graph_dict['pcr'][j]
#         weight = graph_dict['weight'][j]
#         filepath = graph_dict['filepaths'][j]
#         df = pd.read_csv(filepath, sep='\t')
#
#         x_col = "genome_bins" if binning > 0 else "genome_start"
#         fig.add_trace(
#             go.Scattergl(
#                 x=df[x_col],
#                 y=df[frag],
#                 name=f"{samp} - {frag} - {pcr} - {weight}",
#                 mode='lines+markers',
#                 line=dict(width=1, color=traces_colors[trace_id]),
#                 marker=dict(size=4)
#             )
#         )
#
#         fig.update_layout(
#             width=1500,
#             height=500,
#             title=f"Graphe {graph_id}",
#             xaxis=dict(domain=[0.0, 0.9], title="Genome bins"),
#             yaxis=dict(title="Contact frequency"),
#             hovermode='closest'
#         )
#         trace_id += 1
#
#     if x_range:
#         fig.update_xaxes(range=x_range)
#     if y_range:
#         fig.update_yaxes(range=y_range)
#
#     for xi, x_pos in enumerate(chr_boundaries):
#         name_pos = x_pos + 100
#         fig.add_shape(type='line',
#                       yref='paper',
#                       xref='x',
#                       x0=x_pos, x1=x_pos,
#                       y0=0, y1=1,
#                       line=dict(color='gray', width=1, dash='dot'))
#
#         fig.add_annotation(
#             go.layout.Annotation(
#                 x=name_pos,
#                 y=1.07,
#                 yref="paper",
#                 text=chr_names[xi],
#                 showarrow=False,
#                 xanchor="center",
#                 font=dict(size=11, color=chr_colors[xi]),
#                 textangle=330
#             ),
#             xref="x"
#         )
#     return fig


@callback(
    Output('pv-graphs', 'children'),
    [Input('pv-plot-buttom', 'n_clicks'),
     Input('pv-stored-graphs-axis-range', 'data')],
    [State('pv-binning-slider', 'value'),
     State({'type': 'pv-samples-dropdown', 'index': ALL}, 'value'),
     State({'type': 'pv-probes-dropdown', 'index': ALL}, 'value')],
)
def update_graphs(
        n_clicks,
        axis_range,
        binning_value,
        samples_value,
        probes_value,
        data_basedir
):
    ctx = dash.callback_context
    triggerd_input = ctx.triggered[0]['prop_id'].split('.')[0]

    if n_clicks is None or n_clicks == 0:
        return None

    pp_outputs_dir = join(data_basedir, 'outputs')
    graphs_info = {}
    nb_graphs = 0

    x_range = None
    y_range = None
    if triggerd_input == 'pv-stored-graphs-axis-range':
        if axis_range:
            x_range = axis_range['x_range']
            y_range = axis_range['y_range']

    for i, graph in enumerate(graphs_values):
        if graph is None:
            continue
        graph_id = int(graph.split(' ')[-1])
        if graph_id not in graphs_info:
            nb_graphs += 1
            graphs_info[graph_id] = {
                'samples': [],
                'fragments': [],
                'pcr': [],
                'weight': [],
                'filepaths': [],
                'size': 0,
            }
        graphs_info[graph_id]['samples'].append(samples_value[i])
        graphs_info[graph_id]['fragments'].append(probes_value[i])
        graphs_info[graph_id]['pcr'].append(pcr_value[i][-1])
        graphs_info[graph_id]['weight'].append(weight_value[i][-1])

        filedir = join(pp_outputs_dir, samples_value[i], pcr_value[i][-1], weight_value[i][-1])
        if binning_value == 0:
            filepath = join(filedir, f"{samples_value[i]}_unbinned_frequencies.tsv")
            graphs_info[graph_id]['filepaths'].append(filepath)
        else:
            filepath = join(filedir, f"{samples_value[i]}_{binning_value}kb_binned_frequencies.tsv")
            graphs_info[graph_id]['filepaths'].append(filepath)
        graphs_info[graph_id]['size'] += 1

    # TODO: use a file that stores chr data

    chr_cum_pos = list(np.cumsum(chr_pos))
    chr_boundaries = [0] + chr_cum_pos[:-1]

    figures = {}
    traces_count = 0
    for i in graphs_info:
        traces_to_add = graphs_info[i]['size']
        figures[i] = update_figure(
            graph_id=i,
            graph_dict=graphs_info[i],
            traces_colors=colors[traces_count:traces_count + traces_to_add],
            binning=binning_value,
            chr_boundaries=chr_boundaries,
            x_range=x_range,
            y_range=y_range
        )
        traces_count += traces_to_add

    graphs_layout = []
    for i in sorted(figures.keys()):
        graphs_layout.append(
            dcc.Graph(id={'type': 'graph', 'index': i},
                      config={'displayModeBar': True, 'scrollZoom': True},
                      style={'height': 'auto', 'width': '100%'},
                      figure=figures[i])
        )
    return graphs_layout


@callback(
    Output('pv-stored-graphs-axis-range', 'data'),
    Input({'type': 'graph', 'index': ALL}, 'relayoutData'),
    State('pv-sync-box', 'value')
)
def update_figure_range(relayout_data, sync_value):
    if not sync_value:
        return dash.no_update

    if len(relayout_data) == 1:
        return dash.no_update

    if not any(relayout_data):
        return [None, None]

    ctx = dash.callback_context
    input_id = json.loads(ctx.triggered[0]['prop_id'].split('.')[0])['index']

    updated_range = {'x_range': None, 'y_range': None}
    if 'xaxis.range[0]' in relayout_data[input_id]:
        updated_range['x_range'] = [relayout_data[input_id]['xaxis.range[0]'],
                                    relayout_data[input_id]['xaxis.range[1]']]

    if 'yaxis.range[0]' in relayout_data[input_id]:
        updated_range['y_range'] = [relayout_data[input_id]['yaxis.range[0]'],
                                    relayout_data[input_id]['yaxis.range[1]']]
    return updated_range


#
# def create_card(
#         index,
#         sample_options,
#         probe_options,
#         pcr_options,
#         weight_options,
#         graph_options,
#         sample_value,
#         probe_value,
#         pcr_value,
#         weight_value,
#         graph_value
# ):
#     card = dbc.Col(
#         dbc.Card([
#             dbc.CardHeader(html.Div(id={'type': 'pv-probe-card-header', 'index': index})),
#             dbc.CardBody([
#                 dbc.Row([
#                     dbc.Col([
#                         dcc.Dropdown(
#                             options=sample_options,
#                             value=sample_value,
#                             placeholder="Select sample",
#                             id={'type': 'sample-dropdown', 'index': index},
#                             multi=False,
#                         )
#                     ]),
#
#                     dbc.Col([
#                         dcc.Dropdown(
#                             options=probe_options,
#                             value=probe_value,
#                             placeholder="Select probe",
#                             id={'type': 'probe-dropdown', 'index': index},
#                             multi=False,
#                         )
#                     ]),
#                 ]),
#
#                 dbc.Row([
#                     dbc.Col([
#                         dcc.Checklist(
#                             options=pcr_options,
#                             value=pcr_value,
#                             id={'type': 'pcr-checkboxes', 'index': index},
#                             inline=True,
#                             className='custom-checkbox-label',
#                             labelStyle={"margin": "5px"}
#                         )
#                     ]),
#
#                     dbc.Col([
#                         dcc.Checklist(
#                             options=weight_options,
#                             value=weight_value,
#                             id={'type': 'weight-checkboxes', 'index': index},
#                             inline=True,
#                             className='custom-checkbox-label',
#                             labelStyle={"margin": "5px"}
#                         )
#                     ]),
#                 ]),
#
#                 dbc.Row([
#                     dbc.Col([
#                         dcc.Dropdown(
#                             id={'type': 'graph-dropdown', 'index': index},
#                             options=graph_options,
#                             value=graph_value,
#                             placeholder="Select a graph",
#                             multi=False
#                         )
#                     ])
#                 ])
#             ])
#         ])
#     )
#
#     return card
#
#
# @callback(
#     Output('pv-dynamic-probes-cards', 'children'),
#     Input('pv-number-probes', 'value'),
#     State('data-basedir', 'data'),
#     State('pv-dynamic-probes-cards', 'children')
# )
# def update_probes_cards(n_cards, data_basedir, cards_children):
#     if n_cards is None or n_cards == 0:
#         return []
#
#     pp_outputs_dir = join(data_basedir, 'outputs')
#     all_samples_items = sorted(listdir(pp_outputs_dir))
#     samples_options = [{'label': s, 'value': s} for s in all_samples_items if os.path.isdir(join(pp_outputs_dir, s))]
#     graph_options = [{'label': f'graph {x}', 'value': f'graph {x}'} for x in range(n_cards)]
#
#     existing_cards = []
#     displaying_cards = []
#     if cards_children:
#         existing_cards = cards_children[0]['props']['children']
#         for ii, item in enumerate(existing_cards):
#             cardbody = item['props']['children']['props']['children'][1]['props']['children']
#             displaying_cards.append(
#                 create_card(
#                     index=ii,
#                     sample_options=samples_options,
#                     probe_options=cardbody[0]['props']['children'][1]['props']['children'][0]['props']['options'],
#                     pcr_options=cardbody[1]['props']['children'][0]['props']['children'][0]['props']['options'],
#                     weight_options=cardbody[1]['props']['children'][1]['props']['children'][0]['props']['options'],
#                     graph_options=graph_options,
#                     sample_value=cardbody[0]['props']['children'][0]['props']['children'][0]['props']['value'],
#                     probe_value=cardbody[0]['props']['children'][1]['props']['children'][0]['props']['value'],
#                     pcr_value=cardbody[1]['props']['children'][0]['props']['children'][0]['props']['value'],
#                     weight_value=cardbody[1]['props']['children'][1]['props']['children'][0]['props']['value'],
#                     graph_value=cardbody[2]['props']['children'][0]['props']['children'][0]['props']['value']
#                 ))
#
#     if len(existing_cards) > n_cards:
#         displaying_cards = existing_cards[:n_cards]
#     if len(existing_cards) < n_cards:
#         cards_to_add = n_cards - len(existing_cards)
#         for i in range(cards_to_add):
#             displaying_cards.append(create_card(
#                 index=len(existing_cards) + i,
#                 sample_options=samples_options,
#                 probe_options=[],
#                 pcr_options=[],
#                 weight_options=[],
#                 graph_options=graph_options,
#                 sample_value=None,
#                 probe_value=None,
#                 pcr_value=None,
#                 weight_value=None,
#                 graph_value=None
#             ))
#
#     rows = []
#     for i in range(0, len(displaying_cards), 3):
#         row = dbc.Row(displaying_cards[i:i + 3], style={'margin-top': '20px', 'margin-bottom': '20px'})
#         rows.append(row)
#     return rows
#
#
# @callback(
#     Output({'type': 'pcr-checkboxes', 'index': MATCH}, 'options'),
#     Input({'type': 'sample-dropdown', 'index': MATCH}, 'value'),
#     State('data-basedir', 'data')
# )
# def update_pcr_checkboxes_options(sample_value, data_basedir):
#     pp_outputs_dir = join(data_basedir, 'outputs')
#     if sample_value is None:
#         return []
#     all_items = listdir(join(pp_outputs_dir, sample_value))
#     pcr_dirs = [item for item in all_items if os.path.isdir(join(pp_outputs_dir, sample_value, item))]
#     return [{'label': d, 'value': d} for d in pcr_dirs if 'pcr' in d.lower()]
#
#
# @callback(
#     Output({'type': 'pcr-checkboxes', 'index': MATCH}, 'value'),
#     Input({'type': 'pcr-checkboxes', 'index': MATCH}, 'value'),
# )
# def update_pcr_checkboxes(pcr_value):
#     if not pcr_value:
#         return []
#     else:
#         return [pcr_value[-1]]
#
#
# @callback(
#     Output({'type': 'weight-checkboxes', 'index': MATCH}, 'options'),
#     Input({'type': 'pcr-checkboxes', 'index': MATCH}, 'value'),
#     Input({'type': 'sample-dropdown', 'index': MATCH}, 'value'),
#     State('data-basedir', 'data')
# )
# def update_weight_checkboxes_options(pcr_value, sample_value, data_basedir):
#     ctx = dash.callback_context
#     triggerd_input = ctx.triggered[0]['prop_id'].split('.')[0]
#     if triggerd_input == '':
#         return []
#
#     pp_outputs_dir = join(data_basedir, 'outputs')
#     if not sample_value or not pcr_value or pcr_value == []:
#         return []
#
#     pcr_dir = join(pp_outputs_dir, sample_value, pcr_value[-1])
#     weight_dirs = [w for w in listdir(pcr_dir) if os.path.isdir(join(pcr_dir, w))]
#     return [{'label': d, 'value': d} for d in weight_dirs]
#
#
# @callback(
#     Output({'type': 'weight-checkboxes', 'index': MATCH}, 'value'),
#     Input({'type': 'weight-checkboxes', 'index': MATCH}, 'value')
# )
# def update_weight_checkboxes(weight_value):
#     if not weight_value:
#         return []
#     else:
#         return [weight_value[-1]]
#
#
# @callback(
#     Output({'type': 'probe-dropdown', 'index': MATCH}, 'options'),
#     Input({'type': 'weight-checkboxes', 'index': MATCH}, 'value'),
#     Input({'type': 'sample-dropdown', 'index': MATCH}, 'value'),
#     Input({'type': 'pcr-checkboxes', 'index': MATCH}, 'value'),
#     State('data-basedir', 'data')
# )
# def update_probe_dropdown_options(weight_value, sample_value, pcr_value, data_basedir):
#     ctx = dash.callback_context
#     triggerd_input = ctx.triggered[0]['prop_id'].split('.')[0]
#     if triggerd_input == '':
#         return []
#
#     pp_outputs_dir = join(data_basedir, 'outputs')
#     if sample_value is None:
#         return []
#     if pcr_value is None or pcr_value == [] or weight_value is None or weight_value == []:
#         return []
#
#     items_dir = join(pp_outputs_dir, sample_value, pcr_value[-1], weight_value[-1])
#     df = pd.read_csv(join(items_dir, f"{sample_value}_unbinned_contacts.tsv"), sep='\t')
#     probes = [c for c in df.columns if c not in ['chr', 'start', 'sizes', 'genome_start', 'end']]
#     return [{'label': f, 'value': f} for f in probes]
#
#
# @callback(
#     Output({'type': 'pv-probe-card-header', 'index': MATCH}, 'children'),
#     Input({'type': 'probe-dropdown', 'index': MATCH}, 'value'),
#     Input({'type': 'sample-dropdown', 'index': MATCH}, 'value'),
#     Input({'type': 'pcr-checkboxes', 'index': MATCH}, 'value'),
#     Input({'type': 'weight-checkboxes', 'index': MATCH}, 'value')
# )
# def update_card_header(probe_value, sample_value, pcr_value, weight_value):
#     if sample_value is None or probe_value is None:
#         return None
#     if pcr_value is None or pcr_value == [] or weight_value is None or weight_value == []:
#         return None
#
#     samp_id = sample_value.split('/')[-1]
#     return f"{samp_id} - {probe_value} - {pcr_value[-1]} - {weight_value[-1]}"
#
#
# def update_figure(
#         graph_id: int,
#         graph_dict: dict,
#         traces_colors: list,
#         binning: int,
#         chr_boundaries: list,
#         x_range=None,
#         y_range=None
# ):
#     fig = go.Figure()
#     trace_id = 0
#     for j in range(graph_dict['size']):
#         samp = graph_dict['samples'][j]
#         frag = graph_dict['fragments'][j]
#         pcr = graph_dict['pcr'][j]
#         weight = graph_dict['weight'][j]
#         filepath = graph_dict['filepaths'][j]
#         df = pd.read_csv(filepath, sep='\t')
#
#         x_col = "genome_bins" if binning > 0 else "genome_start"
#         fig.add_trace(
#             go.Scattergl(
#                 x=df[x_col],
#                 y=df[frag],
#                 name=f"{samp} - {frag} - {pcr} - {weight}",
#                 mode='lines+markers',
#                 line=dict(width=1, color=traces_colors[trace_id]),
#                 marker=dict(size=4)
#             )
#         )
#
#         fig.update_layout(
#             width=1500,
#             height=500,
#             title=f"Graphe {graph_id}",
#             xaxis=dict(domain=[0.0, 0.9], title="Genome bins"),
#             yaxis=dict(title="Contact frequency"),
#             hovermode='closest'
#         )
#         trace_id += 1
#
#     if x_range:
#         fig.update_xaxes(range=x_range)
#     if y_range:
#         fig.update_yaxes(range=y_range)
#
#     for xi, x_pos in enumerate(chr_boundaries):
#         name_pos = x_pos + 100
#         fig.add_shape(type='line',
#                       yref='paper',
#                       xref='x',
#                       x0=x_pos, x1=x_pos,
#                       y0=0, y1=1,
#                       line=dict(color='gray', width=1, dash='dot'))
#
#         fig.add_annotation(
#             go.layout.Annotation(
#                 x=name_pos,
#                 y=1.07,
#                 yref="paper",
#                 text=chr_names[xi],
#                 showarrow=False,
#                 xanchor="center",
#                 font=dict(size=11, color=chr_colors[xi]),
#                 textangle=330
#             ),
#             xref="x"
#         )
#     return fig
#
#
# @callback(
#     Output('pv-graphs', 'children'),
#     Input('pv-plot-buttom', 'n_clicks'),
#     Input('pv-stored-graphs-axis-range', 'data'),
#     State('pv-binning-slider', 'value'),
#     State({'type': 'sample-dropdown', 'index': ALL}, 'value'),
#     State({'type': 'pcr-checkboxes', 'index': ALL}, 'value'),
#     State({'type': 'weight-checkboxes', 'index': ALL}, 'value'),
#     State({'type': 'probe-dropdown', 'index': ALL}, 'value'),
#     State({'type': 'graph-dropdown', 'index': ALL}, 'value'),
#     State('data-basedir', 'data')
# )
# def update_graphs(
#         n_clicks,
#         axis_range,
#         binning_value,
#         samples_value,
#         pcr_value,
#         weight_value,
#         probes_value,
#         graphs_values,
#         data_basedir
# ):
#     ctx = dash.callback_context
#     triggerd_input = ctx.triggered[0]['prop_id'].split('.')[0]
#
#     if n_clicks is None or n_clicks == 0:
#         return None
#
#     pp_outputs_dir = join(data_basedir, 'outputs')
#     graphs_info = {}
#     nb_cards = len(samples_value)
#     nb_graphs = 0
#
#     x_range = None
#     y_range = None
#     if triggerd_input == 'pv-stored-graphs-axis-range':
#         if axis_range:
#             x_range = axis_range['x_range']
#             y_range = axis_range['y_range']
#
#     for i, graph in enumerate(graphs_values):
#         if graph is None:
#             continue
#         graph_id = int(graph.split(' ')[-1])
#         if graph_id not in graphs_info:
#             nb_graphs += 1
#             graphs_info[graph_id] = {
#                 'samples': [],
#                 'fragments': [],
#                 'pcr': [],
#                 'weight': [],
#                 'filepaths': [],
#                 'size': 0,
#             }
#         graphs_info[graph_id]['samples'].append(samples_value[i])
#         graphs_info[graph_id]['fragments'].append(probes_value[i])
#         graphs_info[graph_id]['pcr'].append(pcr_value[i][-1])
#         graphs_info[graph_id]['weight'].append(weight_value[i][-1])
#
#         filedir = join(pp_outputs_dir, samples_value[i], pcr_value[i][-1], weight_value[i][-1])
#         if binning_value == 0:
#             filepath = join(filedir, f"{samples_value[i]}_unbinned_frequencies.tsv")
#             graphs_info[graph_id]['filepaths'].append(filepath)
#         else:
#             filepath = join(filedir, f"{samples_value[i]}_{binning_value}kb_binned_frequencies.tsv")
#             graphs_info[graph_id]['filepaths'].append(filepath)
#         graphs_info[graph_id]['size'] += 1
#
#     # TODO: use a file that stores chr data
#
#     chr_cum_pos = list(np.cumsum(chr_pos))
#     chr_boundaries = [0] + chr_cum_pos[:-1]
#
#     figures = {}
#     traces_count = 0
#     for i in graphs_info:
#         traces_to_add = graphs_info[i]['size']
#         figures[i] = update_figure(
#             graph_id=i,
#             graph_dict=graphs_info[i],
#             traces_colors=colors[traces_count:traces_count + traces_to_add],
#             binning=binning_value,
#             chr_boundaries=chr_boundaries,
#             x_range=x_range,
#             y_range=y_range
#         )
#         traces_count += traces_to_add
#
#     graphs_layout = []
#     for i in sorted(figures.keys()):
#         graphs_layout.append(
#             dcc.Graph(id={'type': 'graph', 'index': i},
#                       config={'displayModeBar': True, 'scrollZoom': True},
#                       style={'height': 'auto', 'width': '100%'},
#                       figure=figures[i])
#         )
#     return graphs_layout
#
#
# @callback(
#     Output('pv-stored-graphs-axis-range', 'data'),
#     Input({'type': 'graph', 'index': ALL}, 'relayoutData'),
#     State('pv-sync-box', 'value')
# )
# def update_figure_range(relayout_data, sync_value):
#     if not sync_value:
#         return dash.no_update
#
#     if len(relayout_data) == 1:
#         return dash.no_update
#
#     if not any(relayout_data):
#         return [None, None]
#
#     ctx = dash.callback_context
#     input_id = json.loads(ctx.triggered[0]['prop_id'].split('.')[0])['index']
#
#     updated_range = {'x_range': None, 'y_range': None}
#     if 'xaxis.range[0]' in relayout_data[input_id]:
#         updated_range['x_range'] = [relayout_data[input_id]['xaxis.range[0]'],
#                                     relayout_data[input_id]['xaxis.range[1]']]
#
#     if 'yaxis.range[0]' in relayout_data[input_id]:
#         updated_range['y_range'] = [relayout_data[input_id]['yaxis.range[0]'],
#                                     relayout_data[input_id]['yaxis.range[1]']]
#     return updated_range
