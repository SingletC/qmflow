import inspect
from abc import ABC, abstractmethod

from ase import Atoms
from typing import Dict, Any, Callable, List, Union

ProcessDict = Dict[Callable, Union[str, None]]


class Out:
    def __init__(self, *args: Dict[str, Any]):
        self.results: Dict[str, Any] = {}
        for i in args:
            self.results.update(i)

    @property
    def type(self):
        return {key: value if type(value) is type else type(value) for key, value in self.results.items()}

    def get_kwargs_for_in(self, sig: inspect.Signature):
        r = {}
        for var, val in sig.parameters.items():
            r.update({var: self.results[var]})
        return r

    def __repr__(self):
        r = 'Pipe.Out obj:'
        for key, val in self.results.items():
            r += f'{key}   |   {val}' if type(val) is type else f'{key}   |   {type(val)}   |   {val}\n'
        return r


class Process:
    def __init__(self, process_dict: ProcessDict, update=False):
        """

        Args:
            process_dict: next callable
        """
        self.process_dict: Dict[Callable, Union[str, None]] = process_dict
        self.update = update

    def type_check(self, in_: Out):
        incoming_type = in_.type
        for func in self.process_dict:
            sig = inspect.signature(func)
            for var, val in sig.parameters.items():
                if incoming_type.get(var) != val.annotation:
                    return False
        return True

    def process(self, in_: Union[dict, Out, None], type_check=False) -> Out:
        """
        Pipe process method:
        pipe.Out -> [Pipe.Process] ->pipe.OUTPUT
        Beware of thread safe, if parallelize
        """
        if type(in_) is dict:
            in_ = Out(in_)
        results = in_.results if self.update else {}
        for func, maps in self.process_dict.items():
            sig = inspect.signature(func)
            # TODO exception control
            if type_check:
                r = inspect.signature(func).return_annotation
            else:
                if in_:
                    r = func(**in_.get_kwargs_for_in(sig))
                else:
                    r = func()
            if maps is None:
                continue
            else:
                results.update({maps: r})
        return Out(results)

    def __repr__(self):
        r = 'Pipe.Process obj:'
        for key, val in self.process_dict.items():
            r += f'{key}   |   {inspect.signature(key)}'
        return r


class Pipe:
    def __init__(self,
                 *args: Union[Process, ProcessDict],
                 update=False):
        """

        Args:
            update: force all process to be update type
            *args:
        """
        processes = list(args)
        self.processes = [process if type(process) is Process else Process(process) for process in processes]
        if update:
            for process in self.processes:
                process.update = True
        self.counter = 0

    def type_check(self):
        for i, process in enumerate(self.processes[::-1]):
            if not process.type_check(self.processes[i - 2].process(None, True)):
                # raise f'{process} input does not match {self.processes[i - 2]} output'
                return False
        else:
            return True

    def run(self, in_: Union[Out, None, dict] = None, dct: dict = None):
        if type(in_) is dict:
            in_ = Out(in_)
        stage = in_.results.get('stage', 0)
        self.counter = stage
        if in_ is None and type(dct) is dict:
            in_ = Out(dct)
        for process in self.processes[stage:]:
            in_ = process.process(in_)
            self.counter += 1
            in_.results['stage'] = self.counter
        return in_
