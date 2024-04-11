import os
from os.path import join, dirname
import base64

colors_rgba = [
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

colors_hex = ['#000000', '#0c090a', '#2c3e50', '#34495e', '#7f8c8d', '#8e44ad', '#2ecc71', '#2980b9',
              '#f1c40f', '#d35400', '#e74c3c', '#c0392b', '#1abc9c', '#16a085', '#bdc3c7', '#2c3e50',
              '#7f8c8d', '#f39c12', '#27ae60', '#9b59b6', '#3498db', '#e67e22', '#95a5a6', '#d35400',
              '#f1c40f', '#2980b9', '#e74c3c', '#2ecc71', '#8e44ad', '#34495e', '#1abc9c', '#c0392b',
              '#16a085', '#27ae60', '#7f8c8d', '#f39c12', '#bdc3c7', '#000000', '#0c090a', '#2c3e50',
              '#34495e', '#7f8c8d', '#8e44ad', '#2ecc71', '#2980b9', '#f1c40f', '#d35400', '#e74c3c',
              '#c0392b', '#1abc9c', '#16a085', '#bdc3c7', '#2c3e50', '#7f8c8d', '#f39c12', '#27ae60',
              '#9b59b6', '#3498db', '#e67e22', '#95a5a6', '#d35400', '#f1c40f', '#2980b9', '#e74c3c',
              '#2ecc71', '#8e44ad', '#34495e', '#1abc9c', '#c0392b', '#16a085', '#27ae60', '#7f8c8d',
              '#f39c12', '#bdc3c7', '#000000', '#0c090a', '#2c3e50', '#34495e', '#7f8c8d', '#8e44ad',
              '#2ecc71', '#2980b9', '#f1c40f', '#d35400', '#e74c3c', '#c0392b', '#1abc9c', '#16a085',
              '#bdc3c7', '#2c3e50', '#7f8c8d', '#f39c12', '#27ae60', '#9b59b6', '#3498db', '#e67e22',
              '#95a5a6', '#d35400', '#f1c40f', '#2980b9', '#e74c3c', '#2ecc71', '#8e44ad', '#34495e',
              '#1abc9c', '#c0392b', '#16a085', '#27ae60', '#7f8c8d', '#f39c12', '#bdc3c7']

TEMPORARY_DIRECTORY = join(dirname(dirname(os.getcwd())), "data", "__cache__")


def save_file(name, content):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(join(TEMPORARY_DIRECTORY, name), "wb") as fp:
        fp.write(base64.decodebytes(data))


def uploaded_files():
    """List the files in the upload directory."""
    files = []
    for filename in os.listdir(TEMPORARY_DIRECTORY):
        path = join(TEMPORARY_DIRECTORY, filename)
        if os.path.isfile(path):
            files.append(filename)
    return files
