import os
import re

import numpy as np

from ase.units import Hartree, Bohr
from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError


def write_molpro(atoms, **params):
    """Function to write Molpro input file
    """
    charge = params['charge']
    mult = params['mult']
    label = params['label']
    mem = params['memory']
    with open(label + '.inp', 'w') as fd:
        fd.write(f'memory,{mem},m\n')
        fd.write('geomtyp=xyz\n')
        fd.write('geometry={\n')
        fd.write(f'{len(atoms)}\n')
        fd.write('GeomXYZ\n')
        for atom in atoms:
            symbol = atom.symbol + '   '
            fd.write(symbol +
                     str(atom.position[0]) + ' ' +
                     str(atom.position[1]) + ' ' +
                     str(atom.position[2]) + '\n')
        fd.write('}\n')
        fd.write(params['input'])


class Molpro(FileIOCalculator):
    implemented_properties = ['energy', 'forces']
    default_parameters = dict(
        charge=0, mult=1,
        task='gradient',
        input='''basis=vdz
rks, b3lyp

basis=vdz
rhf

basis=VTZ-F12
rhf
CCSD(t)-F12A
FORCE''',
        n=1,
        memory=12800,)

    def __init__(self, restart=None, ignore_bad_restart_file=FileIOCalculator._deprecated, label='molpro', atoms=None,
                 n=1,**kwargs):

        restart = restart
        ignore_bad_restart_file = False
        command = os.getenv('MOLPRO_CMD') or f'molpro -t {n} --ga-impl GA PREFIX.inp > PREFIX.log'
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms,command=command, **kwargs)


    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters
        # p.write(self.label + '.ase')
        p['label'] = self.label
        write_molpro(atoms, **p)

    def read(self, label):
        FileIOCalculator.read(self, label)
        if not os.path.isfile(self.label + '.out'):
            raise ReadError
        self.read_results()

    def read_results(self):
        self.read_energy()
        if self.parameters.task.find('gradient') > -1:
            self.read_forces()

    def read_energy(self):
        """Read Energy from ORCA output file."""
        with open(self.label + '.out', mode='r', encoding='utf-8') as fd:
            text = fd.read()
        # Energy:
        re_energy = re.compile("energy= .*\n")

        self.results['energy'] = float(re_energy.search(text).group().split()[-1]) * Hartree

    def read_forces(self):
        """Read Forces from ORCA output file."""
        with open(self.label + '.out', mode='r', encoding='utf-8') as fd:
            gradients = []
            tempgrad = []
            while True:
                line = fd.readline()
                if not line:
                    break
                if ' Atom          dE/dx               dE/dy               dE/dz' in line:
                    fd.readline()
                    while True:
                        line = fd.readline()
                        if line == '\n':
                            break
                        grad = line.split()
                        tempgrad.append(grad[1:])

            np.array(tempgrad).astype(float)
            self.results['forces'] = -np.array(tempgrad).astype(float) * Hartree/Bohr