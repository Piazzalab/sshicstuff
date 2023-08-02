import os
import re
import pandas as pd
import dash
from dash import html
from dash import dcc
from dash import callback
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
import plotly.graph_objs as go

# Load and set the data folders and files
data_dir = "../../data/samples/pcrfree/AD162/not_pondered"
oligos = "../../data/inputs/capture_oligo_positions.csv"
binned_files = sorted([file for file in os.listdir(data_dir) if file.endswith("_binned_frequencies.tsv")],
                      key=lambda file: int(re.search(r'(\d+)kb', file).group(1)))

# Load binned tables as dataframe
binned_dfs = {file: pd.read_csv(os.path.join(data_dir, file), sep='\t') for file in binned_files}
oligos_df = pd.read_csv(oligos, sep=',')
oligos_df['fragment'] = oligos_df['fragment'].astype(str)

# get list of probes
frag2probe: dict = oligos_df.set_index('fragment')['name'].to_dict()
fragments = list(frag2probe.keys())

layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.Label("Select probes:"),
            dcc.Dropdown(
                id='probe-selector',
                options=[{'label': frag2probe[frag], 'value': frag} for frag in fragments],
                value=[fragments[0], fragments[1]],
                multi=True,
            ),
            html.Label("Select binning:"),
            dcc.Dropdown(
                id='binning-selector',
                options=[{'label': file, 'value': file} for file in binned_files],
                value=binned_files[0]
            ),
            html.Br(),
        ], width=3, style={'position': 'absolute', 'top': '100px', 'left': '25px'}),
    ]),

    dbc.Row([
        dbc.Col([
            dcc.Graph(id='graph1',
                      config={'displayModeBar': True, 'scrollZoom': True},
                      style={'height': 'auto', 'width': '100%'}),
        ], width=12, align='center'),
    ], style={'margin-top': '200px', 'margin-left': '0px'}),

    dbc.Row([
        dbc.Col([
            dcc.Graph(id='graph2',
                      config={'displayModeBar': True, 'scrollZoom': True},
                      style={'height': 'auto', 'width': '100%'}),
        ], width=12, align='center'),
    ], style={'margin-top': '100px', 'margin-left': '0px'})
])


def update_figure(selected_frag, selected_binning, graph_number, x_range=None, y_range=None):
    if not selected_frag:
        return go.Figure()

    df = binned_dfs[selected_binning]
    x_col = "genome_bins"
    fig = go.Figure()

    # Select probe to display based on graph number
    frag_to_display = ""
    if graph_number == 1:
        frag_to_display = selected_frag[0]
    elif graph_number == 2:
        frag_to_display = selected_frag[1]

    # get corresponding probe name for the title
    probe_to_display = frag2probe[frag_to_display]

    fig.add_trace(
        go.Scattergl(
            x=df[x_col],
            y=df[frag_to_display],
            name=probe_to_display,
            mode='lines+markers',
            line=dict(width=1, color='rgba(0,0,255,0.4)' if graph_number == 1 else 'rgba(255,0,0,0.5)'),
            marker=dict(size=4)
        )
    )

    fig.update_layout(
        width=1666,
        height=500,
        title=f"{probe_to_display} contacts frequencies sshic binned at {get_binning_from_filename(selected_binning)}",
        xaxis=dict(domain=[0.0, 0.9], title="Genome bins"),
        yaxis=dict(title="Contact frequency"),
        hovermode='closest'
    )

    if x_range:
        fig.update_xaxes(range=x_range)
    if y_range:
        fig.update_yaxes(range=y_range)

    # Adding a vertical line to show chromosome
    chr_boundaries = [0] + df.loc[df['chr'].shift(-1) != df['chr'], x_col].tolist()
    for boundary in chr_boundaries:
        fig.add_shape(type='line', yref='paper', xref='x', x0=boundary, x1=boundary, y0=0, y1=1,
                      line=dict(color='gray', width=1, dash='dot'))

    # Adding chromosome names between vertical lines
    chr_names = df.loc[df['chr'].shift(-1) != df['chr'], 'chr'].tolist()
    chr_colors = ['#000000', '#0c090a', '#2c3e50', '#34495e', '#7f8c8d', '#8e44ad', '#2ecc71', '#2980b9',
                  '#f1c40f', '#d35400', '#e74c3c', '#c0392b', '#1abc9c', '#16a085', '#bdc3c7', '#2c3e50',
                  '#7f8c8d', '#f39c12', '#27ae60']

    for i, boundary in enumerate(chr_boundaries[:-1]):
        name_pos = boundary + 100
        fig.add_annotation(
            go.layout.Annotation(
                x=name_pos,
                y=1.07,
                yref="paper",
                text=chr_names[i],
                showarrow=False,
                xanchor="center",
                font=dict(size=11, color=chr_colors[i]),
                textangle=330
            ),
            xref="x"
        )
    return fig


def get_binning_from_filename(filename: str):
    match = re.search(r"_(\d+kb)_", filename)
    if match:
        binning = match.group(1)
        return binning


@callback(
    Output('probe-selector', 'value'),
    [Input('probe-selector', 'value')]
)
def limit_probe_selection(selected_probes):
    if len(selected_probes) > 2:
        return selected_probes[2:]
    return selected_probes


@callback(
    [Output('graph1', 'figure'),
     Output('graph2', 'figure')],
    [Input('probe-selector', 'value'),
     Input('binning-selector', 'value'),
     Input('graph1', 'relayoutData'),
     Input('graph2', 'relayoutData')]
)
def sync_graphs(selected_probe, selected_binning, relayout_data1, relayout_data2):
    ctx = dash.callback_context
    input_id = ctx.triggered[0]['prop_id'].split('.')[0]

    x_range = None
    y_range = None

    if input_id == 'graph1':
        relayout_data = relayout_data1
    elif input_id == 'graph2':
        relayout_data = relayout_data2
    else:
        relayout_data = None

    if relayout_data:
        if 'xaxis.range[0]' in relayout_data and 'xaxis.range[1]' in relayout_data:
            x_range = [relayout_data['xaxis.range[0]'], relayout_data['xaxis.range[1]']]
        if 'yaxis.range[0]' in relayout_data and 'yaxis.range[1]' in relayout_data:
            y_range = [relayout_data['yaxis.range[0]'], relayout_data['yaxis.range[1]']]

    if len(selected_probe) >= 1:
        updated_fig1 = update_figure(selected_probe, selected_binning, 1, x_range, y_range)
    else:
        updated_fig1 = go.Figure()

    if len(selected_probe) >= 2:
        updated_fig2 = update_figure(selected_probe, selected_binning, 2, x_range, y_range)
    else:
        updated_fig2 = go.Figure()

    return updated_fig1, updated_fig2
