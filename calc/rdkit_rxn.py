from rdkit import Chem
from rdkit.Chem import AllChem
from ase.atoms import Atoms
from ase.io import read
import io
import numpy as np


def smiles_2_ase(smiles: str) -> Atoms:
    a = Chem.MolFromSmiles(smiles)
    a = Chem.AddHs(a)
    AllChem.EmbedMolecule(a)
    AllChem.MMFFOptimizeMolecule(a)
    string = io.StringIO(Chem.MolToXYZBlock(a))
    ase_atoms = read(string, format='xyz')
    return ase_atoms


def get_r_p_from_smiles(smiles):
    reactant = AllChem.MolFromSmiles(smiles)
    reactant = Chem.AddHs(reactant)
    rxn = AllChem.ReactionFromSmarts('[C:1]1=[C:2][B:3][N:4][C:5]=[C:6]1>>[C:6]1=[C:1][C@@:2]2[B:3][N:4][C@@:5]12')
    ps = rxn.RunReactants([reactant])
    AllChem.MolToSmiles(ps[0][0])
    p = ps[0][0]
    idx_ = [i.GetPropsAsDict()['react_atom_idx'] for i in p.GetAtoms()]
    rev_idx = np.argsort(idx_)
    AllChem.EmbedMolecule(reactant)
    AllChem.UFFOptimizeMolecule(reactant)
    string = io.StringIO(Chem.MolToXYZBlock(reactant))
    r_mol = read(string, format='xyz')

    conf = reactant.GetConformer()
    init_geom = {}
    # Print the 3D coordinates of each atom
    for i, atom in enumerate(reactant.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        init_geom.update({rev_idx[i]: pos})
    Chem.SanitizeMol(p)
    conf = Chem.Conformer(p.GetNumAtoms())
    for i in init_geom:
        conf.SetAtomPosition(int(i), init_geom[i])
    p.AddConformer(conf)
    # AllChem.EmbedMolecule(p,coordMap=init_geom)
    # AllChem.EmbedMolecule(p)
    AllChem.UFFOptimizeMolecule(p)
    string = io.StringIO(Chem.MolToXYZBlock(p))
    p_mol = read(string, format='xyz')[rev_idx]
    import rmsd
    def reorient(A, B):
        A -= rmsd.centroid(A)
        B -= rmsd.centroid(B)
        U = rmsd.kabsch(A, B)
        A = np.dot(A, U)
        return A

    A = r_mol.positions
    B = p_mol.positions
    B = reorient(B, A)
    p_mol.positions = B
    return r_mol, p_mol
