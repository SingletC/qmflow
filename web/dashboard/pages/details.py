import os
import threading
import urllib
from typing import Union

import ase.db
import dash
import dash_bio
from dash import html, dcc, Input, Output, State
from dash_bio.utils import create_mol3d_style

from rdkit import Chem
from rdkit.Chem import AllChem
from ase.atoms import Atoms
from ase.io import read
import io
from ase.visualize import view
from ase.calculators.gaussian import Gaussian, GaussianOptimizer

import uuid
import plotly.graph_objects as go  # or plotly.express as px
from plotly.subplots import make_subplots

from calc.call import determine_nto_type
from calc.utils import gen_uv, get_orbital_text, ase_atoms_to_dash_data, smiles_2_ase, smiles_2_matched_atoms, BNcycle
from web.dashboard.pages.utils import gen_fchk, gen_cube, gen_chargeDiff, gen_NTO, ATOM_COLORS

Molecular_Orbital = 'Molecular Orbital'
# TransitionDensity= 'Transition Density of Excitation State'
Natural_Transition_Orbital = 'Natural Transition Orbital'
ChargeDiff = 'Charge Density Difference'


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
data_empty = ase_atoms_to_dash_data([])
mo_marks = {-i - 1: f'HOMO - {i}' for i in range(0, 4)}
mo_marks.update({i: f'LUMO + {i}' for i in range(0, 4)})
layout = html.Div(children=[
    html.H1(children='Compound Detail page'),
    html.Div([
        "Input Smiles",
        dcc.Input(id='input-on-submit', type='text'),
        dcc.Location(id='url', refresh=False),
        html.Button('check result', id='submit-val', n_clicks=0),
        html.Div([html.Button("Download log file", id="download_btn", disabled=True),
                  html.Br(),
                  html.A(id='downloadfhck', children='Download Gaussian formatted chk file'),
                  dcc.Download(id="download")], ),
    ]),
    html.Br(),
    dcc.Graph(id='UV', style={'width': '40%', 'height': 400, 'display': 'inline-block'}),
    dcc.Textarea(id='orbital-info', value='', style={'width': '40%', 'height': 400, 'display': 'inline-block'}, ),
    html.Div([
        dcc.Dropdown(id='drop-iso-type', options=[Molecular_Orbital, ChargeDiff, Natural_Transition_Orbital]),
        dcc.Dropdown(id='drop-iso-value'),
        dcc.Dropdown(id='drop-iso-type2'),
        dcc.Slider(-3, -1, 0.1,
                   id='slider_iso',
                   marks={i: 'iso-surface value= {}'.format(10 ** i) for i in range(-3, -1)},
                   value=-2,
                   ),
        dcc.Markdown(id='orbital_number'),
        html.Br(),
        dash_bio.Molecule3dViewer(
            style={'height': 500, 'width': '80%'},
            id='dashbio-default-molecule3d',
            modelData=data_empty,
            # orbital=
            styles=create_mol3d_style(
                data_empty['atoms'], visualization_type='stick', color_element='atom'
            )
        ),

    ])
])
db = ase.db.connect(os.getenv('DATABASE_DIR'))


def smiles_2_file(smiles):
    if smiles is None:
        return None
    rdkit_mol = Chem.MolFromSmiles(smiles)
    canonical_smiles = Chem.MolToSmiles(rdkit_mol)
    r = db.select(name=canonical_smiles)
    r = list(r)
    if len(r) == 1:
        return r[0].data.file_path
    else:
        return None
def smiles_2_db_entry(smiles) -> Union[ase.db.core.AtomsRow, None]:
    if smiles is None:
        return None
    rdkit_mol = Chem.MolFromSmiles(smiles)
    canonical_smiles = Chem.MolToSmiles(rdkit_mol)
    r = db.select(name=canonical_smiles)
    r = list(r)
    if len(r) == 1:
        return r[0]
    else:
        return None


@dash.callback(
    Output('UV', 'figure'),
    Output('orbital-info', 'value'),
    Output("download_btn", "disabled"),
    Output("dashbio-default-molecule3d", "modelData"),
    Output("dashbio-default-molecule3d", "styles"),
    Output("downloadfhck", "href"),
    [Input('submit-val', 'n_clicks'),
     State('input-on-submit', 'value'), ],
    prevent_initial_call=True,
)
def update_output(n_clicks, smiles):
    if not smiles:
        return go.Figure(), '', True, None, None, '.'
    file = smiles_2_file(smiles)
    if file is None:
        return go.Figure(), '', True, None, None, '.'
    chk = file + '.chk' if '.log' not in file else file.replace('.log', '.chk')
    log = file + '.log' if '.log' not in file else file
    calc_line, calc_curve = gen_uv(log)
    figure = make_subplots(specs=[[{"secondary_y": True}]])
    figure.add_trace(go.Line(y=calc_line[1], x=calc_line.index, name='Fitted'))
    figure.add_trace(go.Line(y=calc_curve[1], x=calc_curve.index, name='Line', ), secondary_y=True, )
    # figure.update_xaxes(range=[200, 600])
    figure.update_layout(
        title=f"Predicted UV",
        xaxis_title="lambda (nm)",
        yaxis_title="Abs(a.u.)",
        font=dict(family="Courier New, monospace", size=18, ),
        template="plotly_white"
    )
    figure.update_yaxes(title_text="Osc Str.", secondary_y=True)
    orbital = get_orbital_text(log)
    viewer_data = ase_atoms_to_dash_data(read(log))
    style = create_mol3d_style(
        viewer_data['atoms'], visualization_type='stick', color_element='atom', color_scheme=ATOM_COLORS,
    )
    t = threading.Thread(target=gen_fchk, args=[chk])
    t.start()
    return figure, orbital, False, viewer_data, style, '../downloads/' + chk.replace('.chk', '.fchk')


@dash.callback(
    Output("download", "data"),
    Input("download_btn", "n_clicks"),
    State('input-on-submit', 'value'),
    prevent_initial_call=True,
)
def func(n_clicks, smiles):
    file = smiles_2_file(smiles)
    if file is None:
        return None
    log = file + '.log' if '.log' not in file else file
    return dcc.send_file(log, "TD.log")


@dash.callback(
    Output("dashbio-default-molecule3d", "orbital"),
    Output("orbital_number", "children"),
    Output("drop-iso-type2", "options"),
    Input("drop-iso-type", "value"),
    Input("drop-iso-value", "value"),
    Input("drop-iso-type2", "value"),
    Input("slider_iso", "value"),
    State('input-on-submit', 'value'),
    prevent_initial_call=True,
)
def func(iso_type, value, type2, iso, smiles):
    if not smiles or value is None:
        return {}, '', []
    file = smiles_2_file(smiles)
    if file is None:
        return {}, '', []
    log = file + '.log' if '.log' not in file else file
    fchk = log.replace('.log', '.fchk')

    if iso_type == Molecular_Orbital:
        ase_mol = read(log)
        homo = int(sum(ase_mol.get_atomic_numbers()) / 2)
        mo = homo + value + 1  # NOTE we define HOMO is -1
        gen_cube(fchk, mo)
        with open(f'{fchk}{mo}.cube', encoding='utf-8') as f:
            orbital = f.read()
        return {'cube_file': orbital,
                'iso_val': 10 ** iso,
                'opacity': 0.95,
                'positiveVolumetricColor': 'red',
                'negativeVolumetricColor': 'blue',
                }, f'current orbital is {mo}', []
    # if iso_type == TransitionDensity:
    #     txt = gen_hole_electron_cube(fchk,value)
    #     with open(f'{fchk}{value}EH.cube',encoding='utf-8') as f:
    #         cube = f.read()
    #     return {'cube_file':  cube,
    #                     'iso_val': 10**iso,
    #                     'opacity': 0.95,
    #                     'positiveVolumetricColor': 'red',
    #                     'negativeVolumetricColor': 'blue',
    #                 }, f'Showing transition density of excitation state {value}             red(>0) blue(<0)\n'+ txt.replace('\n','\n\n')

    if iso_type == ChargeDiff:
        txt = gen_chargeDiff(fchk, value)
        with open(f'{fchk}{value}CDD.cube', encoding='utf-8') as f:
            cube = f.read()
        return {'cube_file': cube,
                'iso_val': 10 ** iso,
                'opacity': 0.95,
                'positiveVolumetricColor': 'red',
                'negativeVolumetricColor': 'blue',
                }, f'Showing charge density difference of excitation state {value}  \n' \
                   f'electron density red(increase) blue(decrease)\n' + txt.replace(
            '\n', '\n\n'), []
    if iso_type == Natural_Transition_Orbital:
        dict_, orbital_compo,atoms_percent = gen_NTO(fchk, value)
        match_idx = smiles_2_matched_atoms(smiles,BNcycle)
        opts = [{"label": f'MO {i[0]} with {float(i[1]) * 100:0.2f}% contribution', "value": i[0]} for i in
                dict_.items()]
        db_entry = smiles_2_db_entry(smiles)
        if db_entry is not None and value == 1 :
            if db_entry.get('nto_type',0) ==0:
                nto_type = determine_nto_type(orbital_compo)
                db.update(id=db_entry.id, nto_type=nto_type)
        if type2 == None:
            return {}, '', opts
        filename = f'{fchk}{value}NTO.mwfn'
        gen_cube(filename, type2)
        with open(f'{filename}{type2}.cube', encoding='utf-8') as f:
            orbital = f.read()
        return {'cube_file': orbital,
                'iso_val': 10 ** iso,
                'opacity': 0.95,
                'positiveVolumetricColor': 'red',
                'negativeVolumetricColor': 'blue',
                }, f'current orbital is NTO: excitation state {value} Orbital {type2}                   \n ' \
                   f'excited orbital_composition { {type_ : val for type_, val in zip("spdfg",orbital_compo)} } \n'\
            f'sum of NTO excited orbital_composition within BN ring:{atoms_percent[list(match_idx)].sum()} ', opts


@dash.callback(
    Output("drop-iso-value", "options"),
    Input("drop-iso-type", "value"),
    prevent_initial_call=True,
)
def func(value):
    if value == Molecular_Orbital:
        return [{"label": opt, "value": value} for value, opt in sorted(mo_marks.items())]
    if value == Natural_Transition_Orbital:
        return [{"label": f'Excited State {i}', "value": i} for i in range(1, 30)]
    if value == ChargeDiff:
        return [{"label": f'Excited State {i}', "value": i} for i in range(1, 30)]


@dash.callback(Output('input-on-submit', 'value'),
               [Input('url', 'search')])
def update_input(search):
    search = search.replace('$', '#')
    query = urllib.parse.parse_qs(search[1:])
    parameter_value = query.get('smiles', [''])[0]
    return parameter_value
