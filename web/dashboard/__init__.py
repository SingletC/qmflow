"""Instantiate a Dash app."""
import base64
import os
from argparse import ArgumentError

import ase.db
import dash
from dash import dash_table, Output, State, Input
from dash import dcc
from dash import html
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

from calc.call import SubmitTDDFTViaAndromeda
from calc.utils import smiles_2_ase, rdkit_2_base64png
from .data import create_dataframe
from .layout import html_layout
import threading

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
                filter_action="native",
                row_deletable=True,
                page_size=300,
                markdown_options={'link_target': '_blank', "html": True},
                style_cell_conditional=[{'if': {'column_id': 'Structure'},
                                         'width': '400px'}, ]
            ),
        ],
        id="dash-container",
    )
    init_callbacks(dash_app, db)
    return dash_app.server


def init_callbacks(app, db):
    submit_cls = SubmitTDDFTViaAndromeda(db)

    @app.callback(
        Output('container-button-basic', 'children'),
        Input('submit-val', 'n_clicks'),
        State('input-on-submit', 'value'),

        running=[
            (Output('container-button-basic', 'children'),f'Job has been submitted','Job has been submitted')
        ]
    )
    def update_output(n_clicks, smiles):
        if smiles is None:
            return 'Use Chemdraw copy as smiles for input'
        try:
            rdkit_mol = Chem.MolFromSmiles(smiles)
        except ArgumentError:
            return 'Smiles wrong'
        if rdkit_mol is None:
            return 'Smiles wrong'
        AllChem.Compute2DCoords(rdkit_mol)
        img = rdkit_2_base64png(rdkit_mol)
        canonical_smiles = Chem.MolToSmiles(rdkit_mol)
        ase_atom = smiles_2_ase(smiles)
        id_ = db.reserve(name_clean=canonical_smiles)
        if id_ is not None:
            db.update(atoms=ase_atom, id=id_, data={'img': img})
            submit_cls.thread_submit(ase_atom, id_)
            return f'Smiles {smiles} has been submitted'
        else:
            return f'Structure {smiles} submitted early'

    @app.callback(
        Output('database-table', 'data'),

        Input('submit-val', 'n_clicks'),

    )
    def update_table_data(n_clicks):
        """Create Dash datatable from Pandas DataFrame."""
        df = create_dataframe(db)
        return df.to_dict("records")
