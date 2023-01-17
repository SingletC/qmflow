import os

import ase.db
import dash
from dash import html, dcc, callback, Input, Output, State
from dash.dcc import send_bytes

from rdkit import Chem
from rdkit.Chem import AllChem
from ase.atoms import Atoms
from ase.io import read
import io
from ase.visualize import view
from ase.calculators.gaussian import Gaussian, GaussianOptimizer
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import uuid
import plotly.graph_objects as go # or plotly.express as px
import plotly.express as px
from plotly.subplots import make_subplots

from calc.utils import gen_uv, get_orbital_text


def smiles_2_ase(smiles: str) -> Atoms:
    a = Chem.MolFromSmiles(smiles)
    a = Chem.AddHs(a)
    AllChem.EmbedMolecule(a)
    AllChem.MMFFOptimizeMolecule(a)
    string = io.StringIO(Chem.MolToXYZBlock(a))
    ase_atoms = read(string, format='xyz')
    return ase_atoms


def quick_view(atoms: Atoms):
    return view(atoms, viewer='ngl')


def opt_pm7(atoms, steps=1000, method='PM7'):
    try:
        calc_opt = Gaussian(method=f'{method} NoSymmetry', label=f'tmp/{uuid.uuid4().hex[:6].upper()}')
        opt = GaussianOptimizer(atoms, calc_opt)
        opt.run(steps=steps)
        return atoms
    except Exception:
        return None

dash.register_page(__name__)
layout = html.Div(children=[
    html.H1(children='Compound Detail page'),
	html.Div([
        "Input Smiles",
        dcc.Input(id='input-on-submit', type='text'),
    html.Button('check result', id='submit-val', n_clicks=0),
    html.Div([html.Button("Download log file", id="download_btn",disabled=True),
                dcc.Download(id="download")],),
    ]),
	html.Br(),
    dcc.Graph(id='UV',style={'width': '40%', 'height': 400,'display': 'inline-block'}),
    dcc.Textarea(id='orbital-info',value='',style={'width': '40%', 'height': 400,'display': 'inline-block'},)
])
db = ase.db.connect(os.getenv('DATABASE_DIR'))
def smiles_2_file(smiles):
    if smiles is None:
        return None
    r = db.select(name=smiles)
    r = list(r)
    if len(r) ==1:
        return r[0].data.file_path
    else:
        return None

@dash.callback(
    Output('UV', 'figure'),
    Output('orbital-info', 'value'),
    Output("download_btn", "disabled"),
    [Input('submit-val', 'n_clicks'),
    State('input-on-submit', 'value'),],
    prevent_initial_call=True,
)
def update_output(n_clicks, smiles):
    if not smiles:
        return go.Figure() , '' , True
    file = smiles_2_file(smiles)
    if file is None:
        return go.Figure() ,'', True
    chk = file + '.chk' if '.log' not in file else file
    log = file + '.log' if '.log' not in file else file
    calc_line, calc_curve = gen_uv(log)
    figure = make_subplots(specs=[[{"secondary_y": True}]])
    figure.add_trace(go.Line(y=calc_line[1], x=calc_line.index, name='Fitted'))
    figure.add_trace(go.Line(y=calc_curve[1], x=calc_curve.index, name='Line', ), secondary_y=True, )
    # figure.update_xaxes(range=[200, 600])
    figure.update_layout(
        title=f"Exp. vs Calc UV ",
        xaxis_title="lambda (nm)",
        yaxis_title="Abs(a.u.)",
        font=dict(family="Courier New, monospace", size=18, ),
        template="plotly_white"
    )
    figure.update_yaxes(title_text="Osc Str.", secondary_y=True)
    orbital = get_orbital_text(log)
    return figure , orbital , False

@dash.callback(
    Output("download", "data"),
    Input("download_btn", "n_clicks"),
    State('input-on-submit', 'value'),
    prevent_initial_call=True,
)
def func(n_clicks,smiles):
    file = smiles_2_file(smiles)
    if file is None:
        return None
    log = file + '.log' if '.log' not in file else file
    return  dcc.send_file(log, "TD.log")

