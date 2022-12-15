import ase.calculators.calculator
from ase import Atoms


class ASEOperator:
    def __init__(self, calc: ase.calculators.calculator.Calculator):
        self.calc: ase.calculators.calculator.Calculator = calc
        self.back_calc: ase.calculators.calculator.Calculator = calc

    def run(self, atoms: Atoms, label: str):
        atoms = atoms.copy()
        self.calc.label = label
        self.calc.calculate(atoms)
        atoms.calc = self.calc
        return atoms
