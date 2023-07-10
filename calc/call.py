import os
import threading
import traceback
from argparse import ArgumentError
from typing import Protocol, List

import uuid as uuid

import ase.db.core
import pathlib
from ase import Atoms
from ase.calculators.gaussian import Gaussian
from ase.io import read
from rdkit import Chem
from rdkit.Chem import AllChem

from calc.orca import OrcaNEB
from calc.rdkit_rxn import get_r_p_from_smiles
from calc.utils import read_td_dft, smiles_2_ase, smiles_2_base64png, read_gaussian_thermal, get_bn_idx, \
    determine_bncycle_index
from ratelimit import limits, sleep_and_retry, RateLimitException

from flow.operator import ASEOperator
from flow.pipe import Pipe
from web.dashboard.pages.utils import gen_NTO

current_dir = pathlib.Path(__file__).parent.absolute()


class SubmitJobProtocol(Protocol):

    def submit(self, atoms: Atoms, id_) -> bool:
        ...

    def pre_submit(self, atoms: Atoms) -> bool:
        """
        if False stop submit
        """
        ...


class TDDFT_Ase(ASEOperator):

    def setup(self, calc: Gaussian):
        calc.parameters['chk'] = pathlib.Path(calc.label).name + '.chk'

    def __repr__(self):
        return self.__class__.__name__ + f'with calcs :{self.calcs}'


def get_random_string() -> str:
    return './files/' + uuid.uuid4().__str__()


def read_td_dft_from_ase(atoms: Atoms) -> dict:
    log = atoms.calc.label + '.log'
    lambda_, osc_str = read_td_dft(log, t=0.05)
    return {'lambda': lambda_,
            'osc_str': osc_str}


def update_db(db: ase.db.core.Database):
    def inner(id: int, atoms: Atoms, stage: int, uv: dict = None, nto_type: int = 0,
              bn_index: float = 0.0, ) -> None:
        if not uv:
            uv = {'osc_str': 0, 'lambda': 0}
        db.update(id=id, atoms=atoms, lambda_=uv['lambda'], osc_str=uv['osc_str'], nto_type=nto_type, bn_index=bn_index,
                  data={'file_path': atoms.calc.label},
                  stage=stage)

    return inner


def generate_formchk(label):
    command = f'formchk {label}.chk {label}.fchk'
    os.system(command)


def generate_nto_mwfn(label) -> List:
    fchk = f'{label}.fchk'
    orbitals, composition, atoms_percent = gen_NTO(fchk, 1)
    return composition


def determine_nto_type(composition: List) -> int:
    if composition[0] > 50:
        return -1
    if composition[1] > 50:
        return 1
    return 0


class SubmitTDDFTViaAndromeda(SubmitJobProtocol):

    def __init__(self, db: ase.db.core.Database):
        self.db = db

    @sleep_and_retry
    @limits(calls=4, period=60)
    def submit(self, atoms: Atoms, id_) -> bool:
        try:
            functional = os.getenv('functional', 'M06')
            basis = os.getenv('basis', 'jun-cc-pvtz')
            pm7_opt = Gaussian(method=f'opt PM7 IOP(2/9=2000) ', nprocshared=os.getenv('GAUSSIAN_N'), output_type='N',
                               mem=os.getenv('GAUSSIAN_M'), )
            pm3_opt = Gaussian(method=f'opt(maxstep=5 MaxCycles=999) PM3 IOP(2/9=2000) ',
                               nprocshared=os.getenv('GAUSSIAN_N'),
                               output_type='N',
                               mem=os.getenv('GAUSSIAN_M'), )
            pm3_opt.command = os.getenv('GAUSSIAN_CMD')
            pm7_opt.command = os.getenv('GAUSSIAN_CMD')
            opt_calc = Gaussian(method=f'opt {functional}/{basis} IOP(2/9=2000) ', nprocshared=os.getenv('GAUSSIAN_N'),
                                output_type='N',
                                mem=os.getenv('GAUSSIAN_M'), )
            opt_calc.command = os.getenv('GAUSSIAN_CMD')
            opt_calc2 = Gaussian(method=f'opt(maxstep=10) {functional}/{basis} IOP(2/9=2000) ',
                                 nprocshared=os.getenv('GAUSSIAN_N'),
                                 output_type='N',
                                 mem=os.getenv('GAUSSIAN_M'), )
            opt_calc2.command = os.getenv('GAUSSIAN_CMD')
            td_calc = Gaussian(
                method=f'{functional}/{basis} IOP(2/9=2000) scrf=(smd,solvent=cyclohexane)  TD(nstates=30) ',
                nprocshared=os.getenv('GAUSSIAN_N'),
                output_type='N',
                mem=os.getenv('GAUSSIAN_M'))
            td_calc.command = os.getenv('GAUSSIAN_CMD')
            update_db_func = update_db(self.db)
            stage = self.db.get(id=id_).get('stage', 0)
            label = self.db.get(id=id_).data.get('file_path') or get_random_string()
            atoms = atoms.copy() if stage == 0 else read(label + '.log')
            pipe = Pipe({ASEOperator(pm7_opt, pm3_opt).run: 'atoms'},
                        {update_db_func: None},
                        {ASEOperator(opt_calc, opt_calc2).run: 'atoms'},
                        {update_db_func: None},
                        {TDDFT_Ase(td_calc).run: 'atoms'},
                        {read_td_dft_from_ase: 'uv'},
                        # {generate_formchk: None},
                        {determine_bncycle_index: 'bn_index'},
                        {determine_nto_type: 'nto_type'},
                        {update_db_func: None},
                        update=True,
                        )
            pipe.run({'id': id_,
                      'atoms': atoms,
                      'label': label,
                      'uv': None,
                      'bn_index': 0.0,
                      'nto_type': 0,
                      'stage': stage,
                      })
        except RateLimitException:
            raise RateLimitException
        except FileNotFoundError:
            raise RateLimitException  # Make srun error retry
        except Exception as e:
            print(f'error\n {e}\n ')
            traceback.print_tb(e.__traceback__)
        return True

    def thread_submit(self, atoms: Atoms, id_):
        t = threading.Thread(target=self.submit, args=[atoms, id_])
        t.start()

    def smiles_submit(self, smiles):
        if smiles is None:
            return 'Empty Smiles'
        try:
            rdkit_mol = Chem.MolFromSmiles(smiles)
        except ArgumentError:
            return 'Smiles wrong'
        if rdkit_mol is None:
            return 'Smiles wrong'
        img = smiles_2_base64png(smiles)
        AllChem.Compute2DCoords(rdkit_mol)
        canonical_smiles = Chem.MolToSmiles(rdkit_mol)
        ase_atom = smiles_2_ase(smiles)
        id_ = self.db.reserve(name=canonical_smiles)
        if id_ is not None:
            self.db.update(atoms=ase_atom, id=id_, data={'img': img})
            self.thread_submit(ase_atom, id_)
            return f'Smiles {smiles} submitted'
        else:
            return f'Structure {smiles} submitted earlier'


class SubmitKineticViaAndromeda():
    def __init__(self, db: ase.db.core.Database):
        self.db = db

    @sleep_and_retry
    @limits(calls=4, period=60)
    def submit(self, canonical_smiles: str, id_: int) -> bool:
        try:
            method = 'M062X'
            scale = 0.97
            label = self.db.get(id=id_).get('neb_label') or get_random_string()
            self.db.update(id=id_, neb_label=label)
            r_mol, p_mol = get_r_p_from_smiles(canonical_smiles)
            pm7_opt = Gaussian(method=f'{method}/6-311++G(d,p) opt '
                                      f'IOp(2/9=2000)'
                               , nprocshared=os.getenv('GAUSSIAN_N'),
                               output_type='N',
                               mem=os.getenv('GAUSSIAN_M'), label=label,
                               )
            pm7_opt.command = os.getenv('GAUSSIAN_CMD')
            try:
                pm7_opt.calculate(r_mol)
                pm7_opt.calculate(p_mol)
            except Exception as e:
                print(e)
            neb = OrcaNEB('M062X 6-311++G(d,p)', label=label)
            try:
                ts = neb.get_ts()
            except Exception as e:
                neb.run_neb(r_mol, p_mol)
                ts = neb.get_ts()
            r_mol = neb.get_reactant()
            p_mol = neb.get_product()
            opt_calc = Gaussian(method=f'{method} opt'
                                       f' scale={scale} IOp(2/9=2000) freq'
                                , nprocshared=os.getenv('GAUSSIAN_N'),
                                output_type='N',
                                mem=os.getenv('GAUSSIAN_M'),
                                basisfile=str(current_dir / 'MG3S.gbs'))
            opt_calc.command = os.getenv('GAUSSIAN_CMD')
            opt_calc.label = label + '/r'
            opt_calc.calculate(r_mol)
            opt_calc.label = label + '/p'
            opt_calc.calculate(p_mol)
            opt_ts_calc = Gaussian(method=f'{method} opt(ts,calcfc,noeig)'
                                          f' scale={scale} IOp(2/9=2000) freq'
                                   , nprocshared=os.getenv('GAUSSIAN_N'),
                                   output_type='N',
                                   mem=os.getenv('GAUSSIAN_M'),
                                   basisfile=(str(current_dir / 'MG3S.gbs')))
            opt_ts_calc.command = os.getenv('GAUSSIAN_CMD')
            opt_ts_calc.label = label + '/ts'
            opt_ts_calc.calculate(ts)
            r_e = read_gaussian_thermal(label + '/r.log')
            p_e = read_gaussian_thermal(label + '/p.log')
            ts_e = read_gaussian_thermal(label + '/ts.log')
            delta_G = (p_e - r_e) * 627.509
            delta_G_TS = (ts_e - r_e) * 627.509
            self.db.update(id=id_, delta_G=delta_G, delta_G_TS=delta_G_TS, neb_label=label)
        except RateLimitException:
            raise RateLimitException
        except Exception as e:
            print(f'error\n {e}\n ')
            traceback.print_tb(e.__traceback__)
        return True

    def thread_submit(self, canonical_smiles: str, id_: int):
        t = threading.Thread(target=self.submit, args=[canonical_smiles, id_])
        t.start()

    def smiles_submit(self, smiles):
        if smiles is None:
            return 'Empty Smiles'
        try:
            rdkit_mol = Chem.MolFromSmiles(smiles)
        except ArgumentError:
            return 'Smiles wrong'
        if rdkit_mol is None:
            return 'Smiles wrong'
        AllChem.Compute2DCoords(rdkit_mol)
        canonical_smiles = Chem.MolToSmiles(rdkit_mol)
        id_ = self.db.get(name=canonical_smiles).id
        if id_ is None:
            return 'do TD-DFT first'
        elif self.db.get(id=id_).get('delta_G'):
            return 'already calculated'
        else:
            self.db.update(id=id_, delta_G=-1)  # label as ongoing calculation
            self.thread_submit(canonical_smiles, id_)
            return 'submited'
