import os

import ase.db
import dash
from dash import html, dcc, Output, Input, State
from calc.call import SubmitKineticViaAndromeda
app = dash.register_page(__name__)
db = ase.db.connect(os.getenv('DATABASE_DIR'))
sub_cls = SubmitKineticViaAndromeda(db)
# Create Layout
layout = html.Div([
    html.Div(dcc.Input(id='input-on-submit', type='text')),
    html.Button('request', id='submit-val', n_clicks=0),
    html.Div(id='container-button-basic1')
])


@dash.callback(
    Output('container-button-basic1', 'children'),
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
        smiles.replace("'", '')
        return sub_cls.smiles_submit(smiles)
