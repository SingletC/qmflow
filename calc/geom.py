from ase.atoms import Atoms
import rmsd


def calc_rmsd(atoms1: Atoms, atoms2: Atoms) -> float:
    """
    calculate the rmsd between two atoms, after remove rotation and translation
    Args:
        atoms1:
        atoms2:

    Returns:

    """
    x1 = atoms1.get_positions()
    x2 = atoms2.get_positions()
    return rmsd.kabsch_rmsd(x1, x2, translate=True)