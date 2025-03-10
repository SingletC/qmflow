import traceback

import ase.calculators.calculator
import pathlib
from ase import Atoms
from ase.io import read

from flow.exception import ProcessError


class ASEOperator:
    def __init__(self, *calcs: ase.calculators.calculator.Calculator):
        """

        Args:
            *calcs: calcs tuple, if first one fails, try second one
        """
        self.calcs = calcs

    def setup(self, calc: ase.calculators.calculator.Calculator):
        """
        Override to set up additional stuff before call
        Args:
            calc:

        Returns:

        """

        pass

    def post_process(self, atoms: Atoms):
        """
        Override to post calc stuff after call
        Args:
            atoms:
        Returns:
        """
        pass

    def run(self, atoms: Atoms, label: str):
        i = 0
        while i < len(self.calcs):
            try:
                atoms = read(label+'.log')
            except Exception as e :
                atoms = atoms.copy()
            calc = self.calcs[i]
            calc.label = label
            self.setup(calc)
            try:
                calc.calculate(atoms)
                atoms.calc = calc
                self.post_process(atoms)
                return atoms
            except Exception as e:
                print(f'error\n {e}\n ')
                traceback.print_tb(e.__traceback__)
                i += 1
        raise ProcessError(self.__repr__() + 'all fails')

    def __repr__(self):
        return self.__class__.__name__ + f'with calcs :{self.calcs}'
