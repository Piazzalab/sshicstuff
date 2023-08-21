import re
import pandas as pd
import dash
from os.path import join
from os import listdir
from dash import html
from dash import dcc
from dash import callback
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, ALL
import plotly.graph_objs as go


layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.Button(id="pv-add-card-button", className="blue-button", children="Add probe"),
            dbc.Tooltip("Click to select another probe to display",
                        target="pv-add-card-button", className="custom-tooltip", placement="right"),
        ], width=2, className="ml-auto mt-2"),
    ]),
    html.Div(id='pv-dynamic-probes-cards', children=[], style={'margin-top': '20px', 'margin-bottom': '20px'})
])


@callback(
    Output('pv-dynamic-probes-cards', 'children'),
    Output('pv-add-card-button', 'n_clicks'),
    Input('pv-add-card-button', 'n_clicks')
)
def update_probes_cards(n_clicks):
    if n_clicks is None:
        return [], dash.no_update

    probes_cards = []
    for i in range(n_clicks):
        sample_card = dbc.Col(
            dbc.Card([
                dbc.CardHeader(html.Div(id={'type': 'probe-card-header', 'index': i})),
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            dcc.Dropdown(
                                options=[],
                                value=None,
                                placeholder="Select sample",
                                id={'type': 'probe-dropdown', 'index': i},
                                multi=False,
                            )
                        ]),

                        dbc.Col([
                            dcc.Dropdown(
                                options=[],
                                value=None,
                                placeholder="Select probe",
                                id={'type': 'probe-dropdown', 'index': i},
                                multi=False,
                            )
                        ]),
                    ]),

                    dbc.Row([
                        dbc.Col([
                            dcc.Checklist(
                                options=[
                                    {'label': 'pcrfree', 'value': 'pcrfree'},
                                    {'label': 'pcrdupkept', 'value': 'pcrdupkept'}
                                ],
                                value=[],
                                id={'type': 'pcr-checkboxes', 'index': i},
                                inline=True,
                                className='custom-checkbox-label',
                                labelStyle={"margin": "5px"}
                            )
                        ]),

                        dbc.Col([
                            dcc.Checklist(
                                options=[
                                    {'label': 'weighted', 'value': 'weighted'},
                                    {'label': 'not_weighted', 'value': 'not_weighted'}
                                ],
                                value=[],
                                id={'type': 'weight-checkboxes', 'index': i},
                                inline=True,
                                className='custom-checkbox-label',
                                labelStyle={"margin": "5px"}
                            )
                        ]),
                    ]),
                ])
            ])
        )

        probes_cards.append(sample_card)

    rows = []
    for i in range(0, len(probes_cards), 3):
        row = dbc.Row(probes_cards[i:i+3], style={'margin-top': '20px', 'margin-bottom': '20px'})
        rows.append(row)

    return rows, n_clicks


@callback(
    Output({'type': 'pcr-checkboxes', 'index': ALL}, 'value'),
    Input({'type': 'pcr-checkboxes', 'index': ALL}, 'value'),
)
def update_pcr_checkboxes(pcr_values):
    for i in range(len(pcr_values)):
        if not pcr_values[i]:
            continue
        if pcr_values[i][-1] == 'pcrfree':
            pcr_values[i] = ['pcrfree']
        elif pcr_values[i][-1] == 'pcrdupkept':
            pcr_values[i] = ['pcrdupkept']
        else:
            pcr_values[i] = []
    return pcr_values


@callback(
    Output({'type': 'weight-checkboxes', 'index': ALL}, 'value'),
    Input({'type': 'weight-checkboxes', 'index': ALL}, 'value')
)
def update_weight_checkboxes(weight_values):
    for i in range(len(weight_values)):
        if not weight_values[i]:
            continue
        if weight_values[i][-1] == 'weighted':
            weight_values[i] = ['weighted']
        elif weight_values[i][-1] == 'not_weighted':
            weight_values[i] = ['not_weighted']
        else:
            weight_values[i] = []
    return weight_values





# layout = dbc.Container([
#     dbc.Row([
#         dbc.Col([
#             dbc.Card([
#                 dbc.CardHeader("Probe #1"),
#                 dbc.CardBody([
#                     dbc.Row([
#                         dbc.Col([
#                             html.Label("Select PCR mode:",
#                                        style={'margin-top': '0px', 'margin-bottom': '5px'}),
#                             dcc.Dropdown(
#                                 id='pv-pcr-selector-1',
#                                 multi=False,
#                             )
#                         ]),
#                     ]),
#
#                     dbc.Row([
#                         dbc.Col([
#                             html.Label("Select sample:",
#                                        style={'margin-top': '10px', 'margin-bottom': '5px'}),
#                             dcc.Dropdown(
#                                 id='pv-sample-selector-1',
#                                 multi=False,
#                             )
#                         ]),
#                     ]),
#
#                     dbc.Row([
#                         dbc.Col([
#                             html.Label("Select weight :",
#                                        style={'margin-top': '10px', 'margin-bottom': '5px'}),
#                             dcc.Dropdown(
#                                 id='pv-weight-selector-1',
#                                 multi=False,
#                             )
#                         ]),
#                     ]),
#
#                     dbc.Row([
#                         dbc.Col([
#                             html.Label("Select probes :",
#                                        style={'margin-top': '10px', 'margin-bottom': '5px'}),
#                             dcc.Dropdown(
#                                 id='pv-probe-selector-1',
#                                 multi=False,
#                             )
#                         ]),
#                     ]),
#                 ]),
#             ]),
#         ], width=6),
#
#         dbc.Col([
#             dbc.Card([
#                 dbc.CardHeader("Probe #2"),
#                 dbc.CardBody([
#                     dbc.Row([
#                         dbc.Col([
#                             html.Label("Select PCR mode:",
#                                        style={'margin-top': '0px', 'margin-bottom': '5px'}),
#                             dcc.Dropdown(
#                                 id='pv-pcr-selector-2',
#                                 multi=False,
#                             )
#                         ]),
#                     ]),
#
#                     dbc.Row([
#                         dbc.Col([
#                             html.Label("Select sample:",
#                                        style={'margin-top': '10px', 'margin-bottom': '5px'}),
#                             dcc.Dropdown(
#                                 id='pv-sample-selector-2',
#                                 multi=False,
#                             )
#                         ]),
#                     ]),
#
#                     dbc.Row([
#                         dbc.Col([
#                             html.Label("Select weight :",
#                                        style={'margin-top': '10px', 'margin-bottom': '5px'}),
#                             dcc.Dropdown(
#                                 id='pv-weight-selector-2',
#                                 multi=False,
#                             )
#                         ]),
#                     ]),
#
#                     dbc.Row([
#                         dbc.Col([
#                             html.Label("Select probes :",
#                                        style={'margin-top': '10px', 'margin-bottom': '5px'}),
#                             dcc.Dropdown(
#                                 id='pv-probe-selector-2',
#                                 multi=False,
#                             )
#                         ]),
#                     ]),
#                 ]),
#             ]),
#         ], width=6),
#     ], style={'margin-top': '20px', 'margin-bottom': '10px'}),
#
#
#
#     dbc.Row([
#         dbc.Col([
#             html.Label("Select binning :", style={'margin-top': '20px', 'margin-bottom': '20px'}),
#             dcc.Slider(
#                 id='pv-binning-slider',
#                 min=1,
#                 max=100,
#                 step=1,
#                 value=10,
#                 marks={1: "1"} | {i: str(i) for i in range(10, 101, 10)},
#                 included=False
#             )
#         ], width=6),
#     ]),
#
#     dbc.Row([
#         dbc.Col([
#             html.Div(id='pv-slider-output-container'),
#             html.Br(),
#         ], width=4, style={'margin-top': '10px', 'margin-bottom': '0px', 'margin-left': '20px'}),
#     ]),
#
#     dbc.Row([
#         dbc.Col([
#             dcc.Graph(id='pv-graph1',
#                       config={'displayModeBar': True, 'scrollZoom': True},
#                       style={'height': 'auto', 'width': '100%'}),
#         ], width=12, align='center'),
#     ], style={'margin-top': '300px', 'margin-left': '0px'}),
#
#     dbc.Row([
#         dbc.Col([
#             dcc.Graph(id='pv-graph2',
#                       config={'displayModeBar': True, 'scrollZoom': True},
#                       style={'height': 'auto', 'width': '100%'}),
#         ], width=12, align='center'),
#     ], style={'margin-top': '100px', 'margin-left': '0px'})
# ])
#
#
# @callback(
#     Output('pv-slider-output-container', 'children'),
#     [Input('pv-binning-slider', 'value')])
# def update_output(value):
#     return 'You have selected a binning of {} kb'.format(value)
#
#
# @callback(
#     Output('pv-probe-selector', 'options'),
#     Input('this-sample-out-dir-path', 'data')
# )
# def update_probe_selector(sample_out_dir_path):
#     pass


# @callback(
#     Output('weight-selector', 'options'),
#     Input('sample-path', 'data')
# )
# def update_weight_dir(sample_path):
#     if sample_path:
#         data_dir = sample_path
#         weighted_dir = [d for d in listdir(data_dir) if "weighted" in d]
#         return [{'label': d, 'value': d} for d in weighted_dir]
#     return dash.no_update
#
#

# def extract_kbsize(filename):
#     match = re.search(r'(\d+)kb', filename)
#     if match:
#         return match.group(1) + 'kb'
#     return None
#
#
# @callback(
#     Output('probe-selector', 'options'),
#     Output('binning-selector', 'options'),
#     Input('sample-path', 'data'),
# )
# def update_data(sample_path):
#     if sample_path:
#         # Load and set the data folders and files
#         data_dir = sample_path
#         binning_dir = join(data_dir, weight_dir)
#
#         oligos = join(data_dir, 'inputs', 'capture_oligo_positions.csv')
#         binned_files = sorted([file for file in listdir(binning_dir) if file.endswith("_binned_frequencies.tsv")],
#                               key=lambda file: int(re.search(r'(\d+)kb', file).group(1)))
#
#         binning_dict = {}
#         for file in binned_files:
#             kb_size = extract_kbsize(file)
#             if kb_size:
#                 binning_dict[kb_size] = join(binning_dir, file)
#
#         oligos_df = pd.read_csv(oligos, sep=',')
#         oligos_df['fragment'] = oligos_df['fragment'].astype(str)
#
#         # get list of probes
#         frag2probe: dict = oligos_df.set_index('fragment')['name'].to_dict()
#         fragments = list(frag2probe.keys())
#         probes_to_display = [{'label': frag2probe[frag], 'value': frag} for frag in fragments]
#         binning_to_display = [{'label': k, 'value': v} for k, v in binning_dict.items()]
#
#         return probes_to_display, binning_to_display
#     return dash.no_update, dash.no_update
#
#
# layout = dbc.Container([
#     dbc.Row([
#         dbc.Col([
#             html.Label("Select probes :"),
#             dcc.Dropdown(
#                 id='probe-selector',
#                 multi=True,
#             ),
#             html.Label("Select binning :"),
#             dcc.Dropdown(
#                 id='binning-selector',
#                 multi=False,
#             ),
#             html.Br(),
#         ], width=3, style={'position': 'absolute', 'top': '100px', 'left': '25px'}),
#     ]),
#
#     dbc.Row([
#         dbc.Col([
#             dcc.Graph(id='graph1',
#                       config={'displayModeBar': True, 'scrollZoom': True},
#                       style={'height': 'auto', 'width': '100%'}),
#         ], width=12, align='center'),
#     ], style={'margin-top': '250px', 'margin-left': '0px'}),
#
#     dbc.Row([
#         dbc.Col([
#             dcc.Graph(id='graph2',
#                       config={'displayModeBar': True, 'scrollZoom': True},
#                       style={'height': 'auto', 'width': '100%'}),
#         ], width=12, align='center'),
#     ], style={'margin-top': '100px', 'margin-left': '0px'})
# ])
#
#
# @callback(
#     Output('graph', 'figure'),
#     [Input('probe-selector', 'value'),
#      Input('binning-selector', 'value')]
# )
# def update_figure(selected_frag, selected_binning, graph_number, x_range=None, y_range=None):
#     if selected_frag and selected_binning:
#         df = pd.read_csv(selected_binning, sep='\t')
#         x_col = "genome_bins"
#         fig = go.Figure()
#
#         # Select probe to display based on graph number
#         frag_to_display = ""
#         if graph_number == 1:
#             frag_to_display = selected_frag[0]
#         elif graph_number == 2:
#             frag_to_display = selected_frag[1]
#
#         fig.add_trace(
#             go.Scattergl(
#                 x=df[x_col],
#                 y=df[frag_to_display],
#                 name=f"fragment {frag_to_display}",
#                 mode='lines+markers',
#                 line=dict(width=1, color='rgba(0,0,255,0.4)' if graph_number == 1 else 'rgba(255,0,0,0.5)'),
#                 marker=dict(size=4)
#             )
#         )
#
#         fig.update_layout(
#             width=1666,
#             height=500,
#             title=f" fragment {frag_to_display} contacts frequencies sshic binned at "
#                   f"{extract_kbsize(selected_binning)}",
#             xaxis=dict(domain=[0.0, 0.9], title="Genome bins"),
#             yaxis=dict(title="Contact frequency"),
#             hovermode='closest'
#         )
#
#         if x_range:
#             fig.update_xaxes(range=x_range)
#         if y_range:
#             fig.update_yaxes(range=y_range)
#
#         # Adding a vertical line to show chromosome
#         chr_boundaries = [0] + df.loc[df['chr'].shift(-1) != df['chr'], x_col].tolist()
#         for boundary in chr_boundaries:
#             fig.add_shape(type='line', yref='paper', xref='x', x0=boundary, x1=boundary, y0=0, y1=1,
#                           line=dict(color='gray', width=1, dash='dot'))
#
#         # Adding chromosome names between vertical lines
#         chr_names = df.loc[df['chr'].shift(-1) != df['chr'], 'chr'].tolist()
#         chr_colors = ['#000000', '#0c090a', '#2c3e50', '#34495e', '#7f8c8d', '#8e44ad', '#2ecc71', '#2980b9',
#                       '#f1c40f', '#d35400', '#e74c3c', '#c0392b', '#1abc9c', '#16a085', '#bdc3c7', '#2c3e50',
#                       '#7f8c8d', '#f39c12', '#27ae60']
#
#         for i, boundary in enumerate(chr_boundaries[:-1]):
#             name_pos = boundary + 100
#             fig.add_annotation(
#                 go.layout.Annotation(
#                     x=name_pos,
#                     y=1.07,
#                     yref="paper",
#                     text=chr_names[i],
#                     showarrow=False,
#                     xanchor="center",
#                     font=dict(size=11, color=chr_colors[i]),
#                     textangle=330
#                 ),
#                 xref="x"
#             )
#         return fig
#     return dash.no_update
#
#
# @callback(
#     Output('probe-selector', 'value'),
#     [Input('probe-selector', 'value')]
# )
# def limit_probe_selection(selected_probes):
#     if selected_probes:
#         if len(selected_probes) > 2:
#             return selected_probes[2:]
#         return selected_probes
#     return []
#
#
# @callback(
#     [Output('graph1', 'figure'),
#      Output('graph2', 'figure')],
#     [Input('probe-selector', 'value'),
#      Input('binning-selector', 'value'),
#      Input('graph1', 'relayoutData'),
#      Input('graph2', 'relayoutData')]
# )
# def sync_graphs(selected_probe, selected_binning, relayout_data1, relayout_data2):
#     ctx = dash.callback_context
#     input_id = ctx.triggered[0]['prop_id'].split('.')[0]
#
#     x_range = None
#     y_range = None
#
#     if input_id == 'graph1':
#         relayout_data = relayout_data1
#     elif input_id == 'graph2':
#         relayout_data = relayout_data2
#     else:
#         relayout_data = None
#
#     if relayout_data:
#         if 'xaxis.range[0]' in relayout_data and 'xaxis.range[1]' in relayout_data:
#             x_range = [relayout_data['xaxis.range[0]'], relayout_data['xaxis.range[1]']]
#         if 'yaxis.range[0]' in relayout_data and 'yaxis.range[1]' in relayout_data:
#             y_range = [relayout_data['yaxis.range[0]'], relayout_data['yaxis.range[1]']]
#
#     if len(selected_probe) >= 1:
#         updated_fig1 = update_figure(selected_probe, selected_binning, 1, x_range, y_range)
#     else:
#         updated_fig1 = go.Figure()
#
#     if len(selected_probe) >= 2:
#         updated_fig2 = update_figure(selected_probe, selected_binning, 2, x_range, y_range)
#     else:
#         updated_fig2 = go.Figure()
#
#     return updated_fig1, updated_fig2
