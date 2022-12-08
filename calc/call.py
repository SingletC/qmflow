from typing import Protocol

import uuid as uuid

import ase.db.core
from ase import Atoms
from ase.calculators.gaussian import Gaussian
from ase.io import read
from rdkit import Chem
from rdkit.Chem import AllChem

from calc.utils import opt_pm7, read_td_dft


class SubmitJobProtocol(Protocol):

    def submit(self, atoms: Atoms, id_) -> bool:
        ...

    def pre_submit(self, atoms: Atoms) -> bool:
        """
        if False stop submit
        """
        ...


class SubmitTDDFTViaAndromeda(SubmitJobProtocol):
    functional = 'M06'
    basis = 'jun-cc-pvtz'
    opt_calc = Gaussian(method=f'opt {functional}/{basis} IOP(2/9=2000) ', nprocshared=48, mem='100GB', )
    opt_calc.command = 'srun -n g16 48 < PREFIX.com > PREFIX.log'
    td_calc = Gaussian(method=f'{functional}/{basis} IOP(2/9=2000) scrf=(smd,solvent=cyclohexane)  TD(nstates=30) ',
                       nprocshared=48, mem='100GB', )

    def __init__(self, db: ase.db.core.Database):
        self.db = db

    def submit(self, atoms: Atoms, id_) -> bool:
        dir = f'calc/files/{uuid.uuid4().hex[:6].upper()}'
        file = dir + '.log'
        atoms = atoms.copy()
        atoms = opt_pm7(atoms)
        self.opt_calc.label = dir
        self.opt_calc.calculate(atoms=atoms)
        opt = read(file)
        self.td_calc.label = dir
        self.td_calc.calculate(opt)
        # td_mol = read(file)
        r = read_td_dft(file, t=0.05)
        lambda_ = r[0]
        osc_str = r[1]
        self.db.update(id=id_,atoms=opt,lambda_=lambda_,osc_str=osc_str,data={'file_path':file})
        return True
