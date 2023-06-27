import base64
import io
import subprocess
from typing import Tuple, Optional

import numpy as np
import pandas as pd
from ase import neighborlist
from ase.atoms import Atoms
from ase.io import read
from madgp.geoopt.internal_coord.utils import get_conn
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

from web.dashboard.pages.utils import gen_NTO


def smiles_2_ase(smiles: str) -> Atoms:
    a = Chem.MolFromSmiles(smiles)
    a = Chem.AddHs(a)
    AllChem.EmbedMolecule(a)
    AllChem.MMFFOptimizeMolecule(a)
    string = io.StringIO(Chem.MolToXYZBlock(a))
    ase_atoms = read(string, format='xyz')
    return ase_atoms


def get_bn_idx(mol):
    conn = get_conn(mol)
    z_ls = mol.get_atomic_numbers()
    bn = [5, 7, 6, 6, 6, 6]
    bn_copy = bn.copy()
    start = True
    idx_found = []
    wrong_attep = [[] for _ in bn]
    while bn_copy:
        for i, z in enumerate(z_ls):
            # print(i)
            if start:
                pass
            elif not conn[idx_found[-1], i] or i in idx_found or i in wrong_attep[len(idx_found)] or (
                    len(idx_found) == 5 and not conn[idx_found[0], i]):
                continue
            if z == bn_copy[0]:
                if start:
                    start = False
                idx_found += [i]
                bn_copy.remove(z)
                break
        # dead loop
        else:
            wrong_pos = len(idx_found) - 1
            bn_copy.insert(0, bn[wrong_pos])
            wrong_attep[wrong_pos] += [idx_found.pop()]

    return idx_found




def determine_bncycle_index(label) -> float:
    label = label.replace('.log', '')
    fchk = f'{label}.fchk'
    orbitals, composition, atoms_percent = gen_NTO(fchk, 1)
    mol = read(f'{label}.log')
    match_idx = get_bn_idx(mol)
    return atoms_percent[list(match_idx)].sum()


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


def get_linear_fit(df, func):
    x = df['lambda'][func]
    y = df['lambda_exp'][func]
    A = np.vstack([x, np.ones(len(x))]).T
    r = np.linalg.lstsq(A, y, rcond=None)
    # m, c = r[0]
    # residue = r[1]
    # plt.show()
    return r


def smiles_2_base64png(smiles: str) -> str:
    mol_rdkit = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol_rdkit, size=(250, 250))
    buffer = io.BytesIO()
    img.save(buffer, format="png")  # Enregistre l'image dans le buffer
    myimage = buffer.getvalue()
    return base64.b64encode(myimage).decode()


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
        if flag == True:
            text += i

    return text


def get_conn(mol: Atoms, cutoff=1):
    """
    get connectivity from molecule
    Args:
        mol:
        cutoff:

    Returns:

    """
    cutOff = neighborlist.natural_cutoffs(mol, mult=cutoff)
    neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=True, skin=0.0)
    neighborList.update(mol)
    return neighborList.get_connectivity_matrix()


def generate_bonds(mol: Atoms):
    conn = get_conn(mol, cutoff=1.2)
    return conn.nonzero()


def ase_atoms_to_dash_data(mol: Atoms):
    atom = []
    for i, ato in enumerate(mol):
        xyz = ato.position
        atom.append({
            "serial": i,
            "name": ato.symbol,
            "elem": ato.symbol,
            "positions": xyz,
            "mass_magnitude": 1,
            "residue_index": 1,
            "residue_name": 'ALA',
            "chain": 1,

        })
    if mol:
        bonds = generate_bonds(mol)
    else:
        bonds = [[], []]
    return {'atoms': atom,
            'bonds': [{"atom1_index": i, "atom2_index": j, "bond_order": 1} for i, j in zip(bonds[0], bonds[1])]}


def read_gaussian_thermal(log_file):
    try:
        with open(log_file, 'r') as f:
            log_contents = f.readlines()
            for line in log_contents:
                if 'Sum of electronic and thermal Free Energies=' in line:
                    free_energy = float(line.split()[-1])
                    return free_energy
    except:
        return 0.0
