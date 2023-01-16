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


def read_td_dft(log='TD-DFT/0.log', t=0.1):
    with open(log, mode='r') as file:
        while True:
            txt = file.readline()
            if not txt:
                return 0, 0
            if 'Excited State   ' in txt and float(txt.split()[8][2:]) > t:
                return float(txt.split()[6]), float(txt.split()[8][2:])


def gen_uv(log):
    subprocess.run(["Multiwfn", f'{log}'], input=b'''11
3
8
0.4
2''', stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    a = pd.read_csv('./spectrum_curve.txt', delimiter="   ", header=None, index_col=0)
    b = pd.read_csv('./spectrum_line.txt', delimiter="   ", header=None, index_col=0)
    return a, b

def gen_fchk(chk):
    if not os.path.exists(chk+'.fchk'):
        subprocess.run(["formchk", f'{chk}'],stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(chk+'.fchk', "rb") as fh:
        buf = io.BytesIO(fh.read())
    return buf
def get_linear_fit(df, func):
    x = df['lambda'][func]
    y = df['lambda_exp'][func]
    A = np.vstack([x, np.ones(len(x))]).T
    r = np.linalg.lstsq(A, y, rcond=None)
    m, c = r[0]
    residue = r[1]
    # plt.show()
    return r
dash.register_page(__name__)
layout = html.Div(children=[
    html.H1(children='Compound Detail page'),
	html.Div([
        "Input Smiles",
        dcc.Input(id='input-on-submit', type='text'),
    html.Button('Submit', id='submit-val', n_clicks=0),
    ]),
	html.Br(),
    dcc.Graph(id='UV'),
    dcc.Textarea(id='orbital-info',value='',style={'width': '80%', 'height': 600},),
    html.Div([html.Button("Download fchk", id="btn_txt"),
    dcc.Download(id="download")]
                     ),
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
def get_orbital_text(file):
    text = ''
    with open(file) as f:
        lines = f.readlines()
    flag = False
    for i in lines:
        if ' Excitation energies and oscillator strengths:\n' == i:
            flag = True
        if ' SavETr:  write IOETrn=' in i:
            flag = False
        if flag ==True:
            text += i

    return text
@dash.callback(
    Output('UV', 'figure'),
    Output('orbital-info', 'value'),
    [Input('submit-val', 'n_clicks'),
    State('input-on-submit', 'value'),]
)
def update_output(n_clicks, smiles):
    if not smiles:
        return go.Figure() , ''
    file = smiles_2_file(smiles)
    if file is None:
        return go.Figure()
    chk = file + '.chk'
    log = file + '.log'
    calc_line, calc_curve = gen_uv(log)
    figure = make_subplots(specs=[[{"secondary_y": True}]])
    figure.add_trace(go.Line(y=calc_line[1], x=calc_line.index, name='Fitted'))
    figure.add_trace(go.Line(y=calc_curve[1], x=calc_curve.index, name='Line', ), secondary_y=True, )
    figure.update_xaxes(range=[200, 600])
    figure.update_layout(
        title=f"Exp. vs Calc UV ",
        xaxis_title="lambda (nm)",
        yaxis_title="Abs(a.u.)",
        font=dict(family="Courier New, monospace", size=18, ),
        template="plotly_white"
    )
    figure.update_yaxes(title_text="Osc Str.", secondary_y=True)
    orbital = get_orbital_text(log)
    return figure , orbital

@dash.callback(
    Output("download", "data"),
    Input("btn_txt", "n_clicks"),
    State('input-on-submit', 'value'),
    prevent_initial_call=True,
)
def func(n_clicks,smiles):
    file = smiles_2_file(smiles)
    if file is None:
        return None
    f = gen_fchk(file)
    return send_bytes(f.read(), "orbital.fchk")

