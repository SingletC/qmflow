import base64

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
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

from rdkit.Chem.rdchem import Mol


def smiles_2_ase(smiles: str) -> Atoms:
    a = Chem.MolFromSmiles(smiles)
    a = Chem.AddHs(a)
    AllChem.EmbedMolecule(a)
    AllChem.MMFFOptimizeMolecule(a)
    string = io.StringIO(Chem.MolToXYZBlock(a))
    ase_atoms = read(string, format='xyz')
    return ase_atoms


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


def rdkit_2_base64png(mol_rdkit: Mol) -> str:
    img = Draw.MolToImage(mol_rdkit)
    buffer = io.BytesIO()
    img.save(buffer, format="png")  # Enregistre l'image dans le buffer
    myimage = buffer.getvalue()
    return base64.b64encode(myimage).decode()
