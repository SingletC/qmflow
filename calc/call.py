import os
import threading
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
    functional = os.getenv('functional','M06')
    basis = os.getenv('basis', 'jun-cc-pvtz')
    opt_calc = Gaussian(method=f'opt {functional}/{basis} IOP(2/9=2000) ', nprocshared=os.getenv('GAUSSIAN_N'),
                        mem=os.getenv('GAUSSIAN_M'), )
    opt_calc.command = os.getenv('GAUSSIAN_CMD')
    td_calc = Gaussian(method=f'{functional}/{basis} IOP(2/9=2000) scrf=(smd,solvent=cyclohexane)  TD(nstates=30) '
                       , nprocshared=os.getenv('GAUSSIAN_N'),
                       mem=os.getenv('GAUSSIAN_M'), )
    td_calc.command = os.getenv('GAUSSIAN_CMD')

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
        self.db.update(id=id_, atoms=opt, lambda_=lambda_, osc_str=osc_str, data={'file_path': file})
        return True

    def thread_submit(self, atoms: Atoms, id_):
        t = threading.Thread(target=self.submit, args=[atoms, id_])
        t.start()
