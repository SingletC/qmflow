"""Instantiate a Dash app."""

import dash
from dash import dcc
from dash import html

from web.dashboard.pages.utils import split_filter_part


def init_dashboard(server):
    """Create a Plotly Dash dashboard."""
    dash_app = dash.Dash(
        name=__name__,
        use_pages=True,
        server=server,
        routes_pathname_prefix="/dashapp/",
        external_stylesheets=[
            "/static/dist/css/styles.css",
            "https://fonts.googleapis.com/css?family=Lato",
        ],
    )
    # Setup DB
    dash_app.layout = html.Div([
        html.H1(''),

        html.Div(
            [
                html.Div(
                    dcc.Link(
                        f"{page['name']} - {page['path']}", href=page["relative_path"]
                    )
                )
                for page in dash.page_registry.values()
            ]
        ),

        dash.page_container
    ])
    return dash_app