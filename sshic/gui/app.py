import os
import pandas as pd
import dash
from dash import html
from dash import dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go


# Load and set the data folders and files
data_dir = "../../test_data/sshic/AD162/"
binned_files = sorted([file for file in os.listdir(data_dir) if file.endswith("_binned_frequencies.tsv")])

# Load binned tables as dataframe
binned_dfs = {file: pd.read_csv(os.path.join(data_dir, file), sep='\t') for file in binned_files}

# get list of probes
not_a_probe_columns = ["chr", "sizes", "start", "chr_bins", "genome_bins"]
probes = [col for col in binned_dfs[binned_files[0]].columns.values if col not in not_a_probe_columns]

# Create a Dash application instance:
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])


def update_figure(selected_probe, selected_binning, fig):
    if not selected_probe:
        return go.Figure()

    df = binned_dfs[selected_binning]
    x_col = "genome_bins"
    fig = go.Figure(fig)  # Create a copy of the input figure
    fig.update_traces(visible=False)  # Hide existing traces

    probe_traces = [trace for trace in fig.data if trace['name'] == selected_probe]
    if len(probe_traces) == 0:
        # Add trace if it does not already exist
        fig.add_trace(
            go.Scattergl(
                x=df[x_col],
                y=df[selected_probe],
                name=selected_probe,
                mode='lines+markers',
                line=dict(width=1),
                marker=dict(size=4)
            )
        )

    # Show only the trace for the selected probe
    for trace in fig.data:
        if trace['name'] == selected_probe:
            trace.visible = True

    fig.update_layout(
        width=4000,
        height=600,
        xaxis=dict(title="Genome bins"),
        yaxis=dict(title="Contact frequency"),
        hovermode='closest'
    )

    # Adding a vertical line to show chromosome
    chr_boundaries = [0] + df.loc[df['chr'].shift(-1) != df['chr'], x_col].tolist()
    for boundary in chr_boundaries:
        fig.add_shape(type='line', yref='paper', xref='x', x0=boundary, x1=boundary, y0=0, y1=1,
                      line=dict(color='gray', width=1, dash='dot'))

    # Adding chromosome names between vertical lines
    chr_names = df.loc[df['chr'].shift(-1) != df['chr'], 'chr'].tolist()
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
                font=dict(size=10),
                textangle=325
            ),
            xref="x"
        )
    return fig


app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.Label("Select probe:"),
            dcc.Dropdown(
                id='probe-selector',
                options=[{'label': probe, 'value': probe} for probe in probes],
                value=probes[0]
            ),
            html.Label("Select binning:"),
            dcc.Dropdown(
                id='binning-selector',
                options=[{'label': file, 'value': file} for file in binned_files],
                value=binned_files[0]
            ),
            html.Br(),
        ], width=3, style={'position': 'absolute', 'top': '10px', 'left': '10px'}),
    ], style={'margin-top': '10px'}),

    dbc.Row([
        dbc.Col([
            dcc.Graph(id='graph1', config={'displayModeBar': True, 'scrollZoom': True}),
        ], width=12, align='center'),
    ], style={'margin-top': '200px', 'margin-left': '0px'}),

    dbc.Row([
        dbc.Col([
            dcc.Graph(id='graph2', config={'displayModeBar': True, 'scrollZoom': True}),
        ], width=12, align='center'),
    ], style={'margin-top': '200px', 'margin-left': '0px'})
])


@app.callback(
    Output('probe-selector', 'value'),
    [Input('probe-selector', 'value')]
)
def limit_probe_selection(selected_probes):
    if len(selected_probes) > 2:
        return selected_probes[:2]
    return selected_probes


@app.callback(
    Output('graph1', 'figure'),
    [Input('probe-selector', 'value'),
     Input('binning-selector', 'value')]
)
def update_graph1(selected_probe, selected_binning):
    fig = update_figure(selected_probe, selected_binning, None)
    return fig


@app.callback(
    Output('graph2', 'figure'),
    [Input('probe-selector', 'value'),
     Input('binning-selector', 'value')],
    [State('graph1', 'figure')]
)
def update_graph2(selected_probe, selected_binning, graph1_figure):
    fig = update_figure(selected_probe, selected_binning, graph1_figure)
    return fig


if __name__ == '__main__':
    app.run_server(debug=True)
