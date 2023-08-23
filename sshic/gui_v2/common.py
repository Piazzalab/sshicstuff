from dash import dash_table


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

