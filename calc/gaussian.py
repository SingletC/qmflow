from ase import Atoms
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.gaussian import Gaussian
from ase.constraints import FixAtoms
from ase.io import write, read


class Gaussian_fix(Gaussian):
    command = 'GAUSSIAN PREFIX.com'  # avoid redirecting output, buffer may fail

    def write_input(self, atoms, properties=None, system_changes=None):
        if self.parameters.get('extra') is None:
            self.parameters['extra'] = 'IOp(2/9=2000)'
        else:
            self.parameters['extra'] += 'IOp(2/9=2000)'
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        write(self.label + '.com', atoms, properties=properties,
              format='gaussian-in', parallel=False, **self.parameters)
        # read .com file
        if atoms.constraints:
            if atoms.constraints[0].__class__.__name__ == 'FixAtoms':
                with open(self.label + '.com', 'r') as f:
                    # remove tag line
                    lines = f.readlines()
                    flag = 0
                    for idx, line in enumerate(lines):
                        if len(line.split()) ==4:
                            flag +=1
                        elif line =="\n" and flag > 1 and idx > 4:
                            for i in (atoms.constraints[0].get_indices()):
                                lines.insert(idx +1 , f"X {i+1} F\n")
                            break
                # write .com file back
                with open(self.label + '.com', 'w') as f:
                    f.writelines(lines)
if __name__ == "__main__":
    atoms = Atoms('H2O', [[0, 0, 0], [0, 0, 0.74], [0, 0, 1.64]])
    atoms.constraints = [FixAtoms(indices=[0])]
    atoms.calc = Gaussian_fix(method='opt=modredundant')
    atoms.get_potential_energy()

if __name__ == "__main__":
    atoms = read('../Gaussian.log')