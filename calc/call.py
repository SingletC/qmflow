import os
import threading
from argparse import ArgumentError
from typing import Protocol

import uuid as uuid

import ase.db.core
from ase import Atoms
from ase.calculators.gaussian import Gaussian
from ase.io import read
from rdkit import Chem
from rdkit.Chem import AllChem

from calc.utils import opt_pm7, read_td_dft, rdkit_2_base64png, smiles_2_ase
from ratelimit import limits, sleep_and_retry, RateLimitException


class SubmitJobProtocol(Protocol):

    def submit(self, atoms: Atoms, id_) -> bool:
        ...

    def pre_submit(self, atoms: Atoms) -> bool:
        """
        if False stop submit
        """
        ...


class SubmitTDDFTViaAndromeda(SubmitJobProtocol):

    def __init__(self, db: ase.db.core.Database):
        self.db = db

    @sleep_and_retry
    @limits(calls=4, period=60)
    def submit(self, atoms: Atoms, id_) -> bool:
        try:
            functional = os.getenv('functional', 'M06')
            basis = os.getenv('basis', 'jun-cc-pvtz')
            pm7_opt = Gaussian(method=f'opt PM3 IOP(2/9=2000) ', nprocshared=os.getenv('GAUSSIAN_N'),
                               mem=os.getenv('GAUSSIAN_M'), )
            pm7_opt.command = os.getenv('GAUSSIAN_CMD')
            opt_calc = Gaussian(method=f'opt {functional}/{basis} IOP(2/9=2000) ', nprocshared=os.getenv('GAUSSIAN_N'),
                                mem=os.getenv('GAUSSIAN_M'), )
            opt_calc.command = os.getenv('GAUSSIAN_CMD')
            td_calc = Gaussian(
                method=f'{functional}/{basis} IOP(2/9=2000) scrf=(smd,solvent=cyclohexane)  TD(nstates=30) '
                , nprocshared=os.getenv('GAUSSIAN_N'),
                mem=os.getenv('GAUSSIAN_M'), )
            td_calc.command = os.getenv('GAUSSIAN_CMD')
            uuid_ = uuid.uuid4()
            dir = f'calc/files/{uuid_}'
            file = dir + '.log'
            atoms = atoms.copy()
            pm7_opt.label = dir
            pm7_opt.calculate(atoms=atoms)
            pm7 = read(file)
            opt_calc.label = dir
            opt_calc.calculate(atoms=pm7)
            opt = read(file)
            td_calc.label = dir
            td_calc.calculate(atoms=opt)
            # td_mol = read(file)
            r = read_td_dft(file, t=0.05)
            lambda_ = r[0]
            osc_str = r[1]
            self.db.update(id=id_, atoms=opt, lambda_=lambda_, osc_str=osc_str, data={'file_path': file})
        except Exception as e:
            print(f'error\n {e}\n with {dir}')
        except RateLimitException:
            raise RateLimitException
        except FileNotFoundError:
            raise RateLimitException # Make srun error retry
        return True

    def thread_submit(self, atoms: Atoms, id_):
        t = threading.Thread(target=self.submit, args=[atoms, id_])
        t.start()

    def smiles_submit(self,smiles):
        if smiles is None:
            return 'Empty Smiles'
        try:
            rdkit_mol = Chem.MolFromSmiles(smiles)
        except ArgumentError:
            return 'Smiles wrong'
        if rdkit_mol is None:
            return 'Smiles wrong'
        AllChem.Compute2DCoords(rdkit_mol)
        img = rdkit_2_base64png(rdkit_mol)
        canonical_smiles = Chem.MolToSmiles(rdkit_mol)
        ase_atom = smiles_2_ase(smiles)
        id_ = self.db.reserve(name_clean=canonical_smiles)
        if id_ is not None:
            self.db.update(atoms=ase_atom, id=id_, data={'img': img})
            self.thread_submit(ase_atom, id_)
            return f'Smiles {smiles} has been submitted'
        else:
            return f'Structure {smiles} submitted early'
