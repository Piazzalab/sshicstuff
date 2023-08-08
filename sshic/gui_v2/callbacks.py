from os import listdir
from os.path import join, isfile
from dash import callback, no_update
from dash.dependencies import Input, Output


@callback(
    Output('sample-id-output', 'children'),
    Input('this-sample-id', 'data')
)
def display_sample_id(value):
    return f"You are working on sample {value}"


@callback(
    [Output('fragments-list-selector', 'options'),
     Output('oligo-selector', 'options')],
    Input('data-basedir', 'data')
)
def update_input_selectors(data_value):
    if data_value:
        inputs_dir = join(data_value, 'inputs')
        input_files = sorted([file for file in listdir(inputs_dir) if isfile(join(inputs_dir, file))])
        options = [{'label': file, 'value': join(inputs_dir, file)} for file in input_files]
        return options, options
    return no_update, no_update
