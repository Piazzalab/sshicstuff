import os
import re
from os.path import join
import pandas as pd
import numpy as np
import base64

# plotly
from plotly import graph_objs as go


"""
###################
VARIABLES TO IMPORT
###################
"""

__INSTALL_DIR__ = os.path.dirname(os.path.abspath(__file__))
__CACHE_DIR__ = join(__INSTALL_DIR__, "__cache__")

empty_figure = go.Figure(
    layout=go.Layout(
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        annotations=[
            dict(
                x=0.5,
                y=0.5,
                text="No data available",
                showarrow=False,
                font=dict(size=28)
            )
        ],
        hovermode='closest',
        plot_bgcolor='white',
        paper_bgcolor='white'
    )
)


"""
#################
METHODS TO IMPORT
#################
"""


def save_file(name, content):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(join(__CACHE_DIR__, name), "wb") as fp:
        fp.write(base64.decodebytes(data))


def sort_by_chr(df: pd.DataFrame, chr_list: list[str], *args: str):
    chr_list = np.unique(chr_list)
    chr_with_number = [c for c in chr_list if re.match(r'chr\d+', c)]
    chr_with_number.sort(key=lambda x: int(x[3:]))
    chr_without_number = [c for c in chr_list if c not in chr_with_number]

    order = chr_with_number + chr_without_number
    df['chr'] = df['chr'].apply(lambda x: order.index(x) if x in order else len(order))

    if args:
        df = df.sort_values(by=['chr', *args])
    else:
        df = df.sort_values(by=['chr'])

    df['chr'] = df['chr'].map(lambda x: order[x])
    df.index = range(len(df))
    return df


def transform_data(data: np.array, y_max: float, user_y_max: float, y_min: float, re_scale: bool):
    re_scale_output = ""
    if re_scale:
        if y_max <= 1.:
            # squared root transformation
            new_data = np.sqrt(data + 1e-8)
            y_max = np.sqrt(y_max) if not user_y_max else user_y_max
            y_min = np.sqrt(y_min) if y_min > 0 else 0
            re_scale_output = "sqrt"
        else:
            # log transformation
            new_data = np.log(data + 1)
            y_max = np.log(y_max) if not user_y_max else user_y_max
            y_min = 0
            re_scale_output = "log"
    else:
        new_data = data

    return new_data, y_max, y_min, re_scale_output


def uploaded_files():
    """List the files in the upload directory."""
    files = []
    for filename in os.listdir(__CACHE_DIR__):
        path = join(__CACHE_DIR__, filename)
        if os.path.isfile(path):
            files.append(filename)
    return files
