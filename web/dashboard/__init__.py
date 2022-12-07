"""Instantiate a Dash app."""
import os

import ase.db
import dash
from dash import dash_table
from dash import dcc
from dash import html

from .data import create_dataframe
from .layout import html_layout


def init_dashboard(server):
    """Create a Plotly Dash dashboard."""
    dash_app = dash.Dash(
        server=server,
        routes_pathname_prefix="/dashapp/",
        external_stylesheets=[
            "/static/dist/css/styles.css",
            "https://fonts.googleapis.com/css?family=Lato",
        ],
    )
    # Setup DB
    db = ase.db.connect(os.getenv('DATABASE_DIR'))
    # Load DataFrame
    df = create_dataframe(db)

    # Custom HTML layout
    dash_app.index_string = html_layout

    # Create Layout
    dash_app.layout = html.Div(
        children=[
            # dcc.Graph(
            #     id="histogram-graph",
            #     figure={
            #         "data": [
            #             {
            #                 "x": df["complaint_type"],
            #                 "text": df["complaint_type"],
            #                 "customdata": df["key"],
            #                 "name": "311 Calls by region.",
            #                 "type": "histogram",
            #             }
            #         ],
            #         "layout": {
            #             "title": "NYC 311 Calls category.",
            #             "height": 500,
            #             "padding": 150,
            #         },
            #     },
            # ),
            create_data_table(df),
        ],
        id="dash-container",
    )
    return dash_app.server


def create_data_table(df):
    """Create Dash datatable from Pandas DataFrame."""
    table = dash_table.DataTable(
        id="database-table",
        columns=[{"name": i, "id": i,'presentation' : 'markdown' if i=='Structure' else 'input',"deletable": True} for i in df.columns],
        data=df.to_dict("records"),
        sort_action="native",
        sort_mode="native",
        filter_action="native",
        row_deletable=True,
        page_size=300,
        markdown_options={'link_target':'_blank',"html": True},
        style_cell_conditional=[{'if': {'column_id': 'Structure'},
         'width': '400px'},]
    )
    return table
