r'''
    Script to generate data on measuring time for the Grover problem

    This is used for generating data on Table 1 on the Quantum draft (see :arxiv:`2308.09510`).
'''
from __future__ import annotations
import sys, os

SCRIPT_DIR = os.path.dirname(__file__) if __name__ != "__main__" else "./"
sys.path.insert(0, os.path.join(SCRIPT_DIR, "..", "..")) # clue is here
SCRIPT_NAME = os.path.splitext(os.path.basename(__file__))[0]

from clue.linalg import CC, SparseRowMatrix, SparseVector
from itertools import product
from numpy import pi
from qiskit import AncillaRegister, QuantumCircuit, QuantumRegister
from qiskit.circuit.library import GroverOperator
from random import randint

## Imports from the local folder
from misc import *

class QuantumSearch(Experiment):
    def __init__(self, size, oracle, *, _success: list[int]=None):
        self.__size = size
        self.__oracle = oracle
        self.__success = _success

    @staticmethod
    def random(size) -> QuantumSearch:
        success = set(sorted([randint(0, 2**size-1) for _ in range(randint(1,2*size))])) # sparse success search, but random
        return QuantumSearch(size, lambda *p : 1 if int("".join(str(el) for el in p), 2) in success else 0, _success = success)

    @property
    def qbits(self) -> int: return self.__size

    @property
    def oracle(self): return self.__oracle

    def grover_matrix(self): 
        M = SparseRowMatrix(2**self.qbits, CC)
        for i,x in enumerate(product(range(2), repeat=self.qbits)):
            for j in range(2**self.qbits):
                M.increment(j, i, (-1)**(1+int(i==j)+int(self.oracle(*x)))*0.5)

        return M
    
    def grover_preparation(self) -> QuantumCircuit:
        raise NotImplementedError

    def grover_quantum(self):
        raise NotImplementedError
    
    def __repr__(self) -> str: return f"Grover for a search with {self.qbits} q-bits" + (f" {self.__success}" if self.__success != None else "")

    def size(self) -> int: return self.qbits
    def correct_size(self) -> int: return 2
    def matrix(self) -> SparseRowMatrix: return self.grover_matrix()
    def quantum(self) -> tuple[QuantumCircuit, Parameter]: return self.grover_quantum(), None
    def data(self): return []

class ToffoliSearch(QuantumSearch):
    def __init__(self, size):
        super().__init__(size, lambda *p : all(p))

    @staticmethod
    def random(size:int):
        return ToffoliSearch(size)
    
    def grover_preparation(self) -> QuantumCircuit:
        q = QuantumRegister(self.size(), "q")
        flag = AncillaRegister(1, "flag")
        state_preparation = QuantumCircuit(q, flag)
        state_preparation.h(q)
        state_preparation.x(flag)

        return state_preparation

    def grover_quantum(self):
        q = QuantumRegister(self.size(), "q")
        flag = AncillaRegister(1, "flag")

        oracle = QuantumCircuit(q, flag)
        oracle.mcp(pi, q, flag)

        grover = GroverOperator(oracle, mcx_mode="noancilla")

        return grover
    
    def __repr__(self) -> str: return f"Grover over Toffoli with {self.size()} q-bits"

def generate_example(_: str, size: int) -> ToffoliSearch:
    return ToffoliSearch.random(size)

def generate_header(csv_writer, ttype):
    if ttype in ("clue", "ddsim"):
        csv_writer.writerow(["size", "red. ratio", "time_lumping", "memory (MB)", "gate"])
    elif ttype == "full_clue":
        csv_writer.writerow(["size", "time_lumping", "kappa", "time_iteration", "memory (MB)", "gate"])
    elif ttype == "full_ddsim":
        csv_writer.writerow(["size", "kappa", "time_iteration", "memory (MB)", "gate"])
    else:
        raise NotImplementedError(f"Type of file {ttype} not recognized")

## METHODS TO GENERATE OBSERVABLES
def generate_observable_clue(search: QuantumSearch, *_) -> tuple[SparseVector]:
    return tuple([SparseVector.from_list(2**search.size()*[1], field=CC)])

def generate_observable_ddsim(search: QuantumSearch, *_) -> bool:
    return search.grover_preparation()

if __name__ == "__main__":
    ## Processing the arguments
    ## Generic part
    ttype, script = get_method(*sys.argv)
    m, M = get_size_bounds(*sys.argv)
    timeout = get_timeout(*sys.argv)
    repeats = get_repeats(*sys.argv)
    name = "grover"; obs = None

    main_script(SCRIPT_DIR, SCRIPT_NAME, name,
                [generate_header, generate_example, generate_observable_clue, generate_observable_ddsim],
                ttype, script,               
                m, M,                        
                timeout, repeats,
                obs)