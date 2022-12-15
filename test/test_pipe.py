import unittest
import uuid

import numpy as np
from ase.calculators.gaussian import Gaussian

from flow.operator import ASEOperator
from flow.pipe import Out, Process, Pipe
from ase.atoms import Atoms

atoms = Atoms('C')
symbols = 'C'


def atoms_read(atoms: Atoms) -> None:
    return None


def atoms_read2(atoms: Atoms, atoms2: Atoms) -> None:
    return None


def string_read(key: str) -> None:
    return None


def null_2_atoms() -> Atoms:
    return atoms.copy()


def string_2_atoms(symbols: str) -> Atoms:
    return Atoms(symbols)


def atoms_2_string(atoms: Atoms) -> str:
    return str(atoms.symbols)


def atoms_2_atoms(atoms: Atoms) -> Atoms:
    """ important test function. This is ASE operator like
    """
    return atoms


def atoms_string_2_string(atoms: Atoms, symbols: str) -> str:
    return f'symbols string is {symbols} and atoms symbol is ' + str(atoms.symbols)


def get_ase_forces(atoms: Atoms) -> np.ndarray:
    return atoms.get_forces()


def get_ase_energy(atoms: Atoms) -> float:
    return atoms.get_potential_energy()


def get_random_string() -> str:
    return './tmp/'+uuid.uuid4().__str__()


class MyTestCase(unittest.TestCase):

    def test_result_read_type_check_pass(self):
        result = Out({'atoms': Atoms})
        pipe_read = Process({atoms_read: None})
        self.assertTrue(pipe_read.type_check(result))

        result = Out({'atoms': Atoms, 'atoms2': Atoms})
        pipe_read = Process({atoms_read2: None})
        self.assertTrue(pipe_read.type_check(result))
        # merge results
        result = Out({'key': str, 'atoms': Atoms})
        pipe_read = Process({atoms_read: None,
                             string_read: None})
        self.assertTrue(pipe_read.type_check(result))
        # redundant key
        result = Out({'key': str, 'atoms': Atoms})
        pipe_read = Process({atoms_read: None, })
        self.assertTrue(pipe_read.type_check(result))

    def test_result_read_type_check_fail(self):
        # Wrong result type
        result = Out({'atoms': str})
        pipe_read = Process({atoms_read: None})
        self.assertFalse(pipe_read.type_check(result))
        # Wrong result amount
        result = Out({'atoms': Atoms})
        pipe_read = Process({atoms_read2: None})
        self.assertFalse(pipe_read.type_check(result))
        # Empty result
        result = Out({})
        pipe_read = Process({atoms_read2: None})
        self.assertFalse(pipe_read.type_check(result))
        # Wrong result keys
        result = Out({'result': Atoms})
        pipe_read = Process({atoms_read2: None})
        self.assertFalse(pipe_read.type_check(result))

    def test_process_output(self):
        inp = Out({'symbols': symbols})
        process1 = Process({string_2_atoms: 'atoms'})
        result1 = process1.process(inp)
        process2 = Process({atoms_2_string: 'symbols'})
        self.assertEqual(process2.process(result1).results['symbols'],
                         symbols)

        inp = Out({'symbols': symbols, 'atoms': atoms})
        process = Process({atoms_string_2_string: 'r'})
        print(process)
        self.assertEqual(process.process(inp).results['r'],
                         atoms_string_2_string(atoms, symbols))

    def test_pipe_type_check(self):
        process1 = Process({string_2_atoms: 'atoms'})
        process2 = Process({atoms_2_string: 'symbols'})
        pipe = Pipe(process1, process2)
        self.assertTrue(pipe.type_check())

        process1 = Process({atoms_2_atoms: 'atoms'})
        process2 = Process({atoms_2_atoms: 'atoms'})
        pipe = Pipe(process1, process2)
        self.assertTrue(pipe.type_check())
        # Null-> Any
        process1 = Process({null_2_atoms: 'atoms'})
        process2 = Process({atoms_2_atoms: 'atoms'})
        pipe = Pipe(process1, process2)
        self.assertTrue(pipe.type_check())
        """     Atoms > Atoms
        Atoms->                -> Pipe.Out     # Fork type pipe
                Atoms > str
        """
        process1 = Process({atoms_2_atoms: 'atoms'})
        process2 = Process({atoms_2_string: 'symbols',
                            atoms_2_atoms: 'atoms', })
        pipe = Pipe(process1, process2)
        self.assertTrue(pipe.type_check())

    def test_pipe_type_check_fail(self):
        # Wrong args name
        process1 = Process({string_2_atoms: 'xx'})
        process2 = Process({atoms_2_string: 'symbols'})
        pipe = Pipe(process1, process2)
        self.assertFalse(pipe.type_check())
        # Wrong args type
        process1 = Process({string_2_atoms: 'xx'})
        process2 = Process({string_2_atoms: 'symbols'})
        pipe = Pipe(process1, process2)
        self.assertFalse(pipe.type_check())

        process1 = Process({atoms_2_atoms: 'atoms'})
        process2 = Process({string_2_atoms: 'symbols',  # Wrong type
                            atoms_2_atoms: 'atoms', })
        pipe = Pipe(process1, process2)
        self.assertFalse(pipe.type_check())

        process1 = Process({atoms_2_atoms: 'atoms'})
        process2 = Process({atoms_2_atoms: 'symbols',
                            string_2_atoms: 'atoms', })  # Wrong type
        pipe = Pipe(process1, process2)
        self.assertFalse(pipe.type_check())

    def test_pipe_run(self):
        """     Atoms > Atoms
        Atoms->                -> Pipe.Out     # Fork type pipe
                Atoms > str
        """
        input_ = {'atoms': atoms}
        process1 = Process({atoms_2_atoms: 'atoms'})
        process2 = Process({atoms_2_string: 'symbols',
                            atoms_2_atoms: 'atoms', })
        pipe = Pipe(process1, process2)
        results = pipe.run(dct=input_).results
        self.assertEqual(results,
                         {'atoms': atoms,
                          'symbols': symbols})

        # Null-> Any
        process1 = Process({null_2_atoms: 'atoms'})
        process2 = Process({atoms_2_atoms: 'atoms'})
        pipe = Pipe(process1, process2)
        results = pipe.run().results
        self.assertEqual(results, {'atoms': atoms})

    def test_ase_process(self):
        calc = Gaussian(method='pm7 opt',chk='td-chk')
        process1 = Process({ASEOperator(calc).run: 'atoms'})
        process2 = Process({get_ase_forces: 'force'})
        input_ = {'atoms': Atoms('C'),
                  'label': './tmp/gaussian.log'}
        # pipe = Pipe(process1, process2)
        # print(pipe.run(input_).results['force'])
        #
        # # Direct use Dict to define Pipe
        # pipe = Pipe({},
        #             {ASEOperator(calc).run: 'atoms'},
        #             {get_ase_forces: 'force', get_ase_energy: 'energy'})
        # print(pipe.run(input_).results)

        # Update type process
        pipe = Pipe(Process({get_random_string: 'label'}, update=True),
                    Process({ASEOperator(calc).run: 'atoms'}, update=True),
                    Process({get_ase_forces: 'force', get_ase_energy: 'energy'}, update=True))
        print(pipe.run(input_).results)


if __name__ == '__main__':
    unittest.main()
