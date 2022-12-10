"""Instantiate a Dash app."""
import base64
import os
import time
from argparse import ArgumentError

import ase.db
import dash
from dash import dash_table, Output, State, Input
from dash import dcc
from dash import html
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, MolFromSmiles

from calc.call import SubmitTDDFTViaAndromeda
from calc.utils import smiles_2_ase, rdkit_2_base64png
from .data import create_dataframe
from .layout import html_layout
import threading

from .utils import split_filter_part


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
            html.Div([
                html.Div(dcc.Input(id='input-on-submit', type='text')),
                html.Button('request', id='submit-val', n_clicks=0),
                html.Div(id='container-button-basic')
            ]),

            dash_table.DataTable(
                id="database-table",
                columns=[
                    {"name": i, "id": i, 'presentation': 'markdown' if i == 'Structure' else 'input', "deletable": True}
                    for
                    i in df.columns],
                data=df.to_dict("records"),
                sort_action="native",
                sort_mode="native",
                filter_action="custom",
                filter_query='',
                row_deletable=True,
                page_size=300,
                markdown_options={'link_target': '_blank', "html": True},
                style_cell_conditional=[{'if': {'column_id': 'Structure'},
                                         'width': '400px'},
                                        {'if': {'column_id': 'name'},
                                         'maxWidth': '100px'},
                                        ],

            ),
        ],
        id="dash-container",
    )
    callbacks = CallBacks(dash_app,db)
    callbacks.init_callbacks()
    return dash_app.server


class CallBacks:
    def __init__(self, app: dash.Dash, db: ase.db.core.Database):
        self.db = db
        self.app = app
        self.submit_cls = SubmitTDDFTViaAndromeda(db)
        self.df = create_dataframe(db)
        self.update_df_thead = threading.Thread(target=self.update_df)
        self.update_df_thead.start()
    def init_callbacks(self):

        @self.app.callback(
            Output('container-button-basic', 'children'),
            Input('submit-val', 'n_clicks'),
            State('input-on-submit', 'value'),

            running=[
                (Output('container-button-basic', 'children'), f'Job has been submitted', 'Job has been submitted')
            ]
        )
        def update_output(n_clicks, smiles):
            if smiles is None:
                return 'Use Chemdraw copy as smiles for input'
            else:

                return self.submit_cls.smiles_submit(smiles)

        # @self.app.callback(
        #     Output('database-table', 'data'),
        #
        #     Input('submit-val', 'n_clicks'),
        #
        # )
        # def update_table_data(n_clicks):
        #     """Create Dash datatable from Pandas DataFrame."""
        #     self.df = create_dataframe(self.db)
        #     return self.df .to_dict("records")

        @self.app.callback(
            Output('database-table', "data"),
            Input('database-table', "filter_query"))
        def update_table(filter):
            filtering_expressions = filter.split(' && ')
            dff = self.df.copy()
            for filter_part in filtering_expressions:
                if '{name}' in filter_part:
                    col_name = 'name'
                    operator = 'contains'
                    filter_value = filter_part.split(' ')[-1]
                else:
                    col_name, operator, filter_value = split_filter_part(filter_part)

                if operator in ('eq', 'ne', 'lt', 'le', 'gt', 'ge'):
                    # these operators match pandas series operator method names
                    dff = dff.loc[getattr(dff[col_name], operator)(filter_value)]
                elif operator == 'contains':
                    if col_name == 'name':
                        match = Chem.MolFromSmiles(filter_value)
                        p = Chem.AdjustQueryParameters.NoAdjustments()
                        p.makeDummiesQueries = True
                        match = Chem.AdjustQueryProperties(match, p)
                        dff = dff[dff['name'].apply(lambda x: MolFromSmiles(x).HasSubstructMatch(match))]

                    else:
                        dff = dff.loc[dff[col_name].str.contains(filter_value)]
                elif operator == 'datestartswith':
                    # this is a simplification of the front-end filtering logic,
                    # only works with complete fields in standard format
                    dff = dff.loc[dff[col_name].str.startswith(filter_value)]

            return dff.to_dict('records')

    def update_df(self):
        while True:
            self.df = create_dataframe(self.db)
            time.sleep(360)