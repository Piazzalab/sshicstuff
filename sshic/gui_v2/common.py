from dash import dash_table


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

