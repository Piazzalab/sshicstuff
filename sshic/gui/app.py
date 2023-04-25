import os
import re
import pandas as pd
import dash
from dash import html
from dash import dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
import plotly.graph_objs as go


# Create a Dash application instance:
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.index_string = '''
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
        <style>
            body {{
                overflow-x: scroll;
            }}
        </style>
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
'''

# Load and set the data folders and files
data_dir =  "../../test_data/AD162_classic/AD162"
binned_files = sorted([file for file in os.listdir(data_dir) if file.endswith("_binned_frequencies.tsv")],
                      key=lambda file: int(re.search(r'(\d+)kb', file).group(1)))

# Load binned tables as dataframe
binned_dfs = {file: pd.read_csv(os.path.join(data_dir, file), sep='\t') for file in binned_files}

# get list of probes
not_a_probe_columns = ["chr", "sizes", "start", "chr_bins", "genome_bins"]
probes = [col for col in binned_dfs[binned_files[0]].columns.values if col not in not_a_probe_columns]

app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.Label("Select probes:"),
            dcc.Dropdown(
                id='probe-selector',
                options=[{'label': probe, 'value': probe} for probe in probes],
                value=[probes[0], probes[1]],
                multi=True,
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


def update_figure(selected_probe, selected_binning, graph_number):
    if not selected_probe:
        return go.Figure()

    df = binned_dfs[selected_binning]
    x_col = "genome_bins"
    fig = go.Figure()

    # Select probe to display based on graph number
    probe_to_display = ""
    if graph_number == 1:
        probe_to_display = selected_probe[0]
    elif graph_number == 2:
        probe_to_display = selected_probe[1]

    fig.add_trace(
        go.Scattergl(
            x=df[x_col],
            y=df[probe_to_display],
            name=probe_to_display,
            mode='lines+markers',
            line=dict(width=1, color='rgba(0,0,255,0.4)' if graph_number == 1 else 'rgba(255,0,0,0.5)'),
            marker=dict(size=4)
        )
    )

    fig.update_layout(
        width=2000,
        height=600,
        title=f"{probe_to_display} contacts frequencies sshic binned at {get_binning_from_filename(selected_binning)}",
        xaxis=dict(domain=[0.0, 0.9], title="Genome bins"),
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


@app.callback(
    Output('probe-selector', 'value'),
    [Input('probe-selector', 'value')]
)
def limit_probe_selection(selected_probes):
    if len(selected_probes) > 2:
        return selected_probes[2:]
    return selected_probes


@app.callback(
    Output('graph1', 'figure'),
    [Input('probe-selector', 'value'),
     Input('binning-selector', 'value')]
)
def update_graph1(selected_probe, selected_binning):
    fig = update_figure(selected_probe, selected_binning, 1)
    return fig


@app.callback(
    Output('graph2', 'figure'),
    [Input('probe-selector', 'value'),
     Input('binning-selector', 'value')]
)
def update_graph2(selected_probe, selected_binning):
    fig = update_figure(selected_probe, selected_binning, 2)
    return fig


if __name__ == '__main__':
    app.run_server(debug=True)
