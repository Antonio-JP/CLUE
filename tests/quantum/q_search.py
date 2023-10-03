r'''
    Script to generate data on measuring time for the Grover problem

    This is used for generating data on Table 1 on the Quantum draft (see :arxiv:`2308.09510`).
'''
from __future__ import annotations
import sys, os

SCRIPT_DIR = os.path.dirname(__file__) if __name__ != "__main__" else "./"
sys.path.insert(0, os.path.join(SCRIPT_DIR, "..", "..")) # clue is here

import tracemalloc
from clue import FODESystem
from clue.linalg import CC, NumericalSubspace, SparseRowMatrix, SparseVector
from csv import writer
from itertools import product
from math import ceil, sqrt
from mqt import ddsim #pylint: disable=no-name-in-module
from numpy import cdouble, pi, Inf
from numpy.linalg import matrix_power
from qiskit import AncillaRegister, QuantumCircuit, QuantumRegister, execute
from qiskit.circuit.library import GroverOperator
from random import randint
from time import time

## Imports from the local folder
from misc import *

class QuantumSearch:
    def __init__(self, size, oracle, *, _success: list[int]=None):
        self.__size = size
        self.__oracle = oracle
        self.__success = _success

    @staticmethod
    def random(size) -> QuantumSearch:
        success = set(sorted([randint(0, 2**size-1) for _ in range(randint(1,2*size))])) # sparse success search, but random
        return QuantumSearch(size, lambda *p : 1 if int("".join(str(el) for el in p), 2) in success else 0, _success = success)

    @property
    def size(self) -> int: return self.__size

    @property
    def oracle(self): return self.__oracle

    def grover_matrix(self): 
        M = SparseRowMatrix(2**self.size, CC)
        for i,x in enumerate(product(range(2), repeat=self.size)):
            for j in range(2**self.size):
                M.increment(j, i, (-1)**(1+int(i==j)+int(self.oracle(*x)))*0.5)

        return M
    
    def grover_quantum(self):
        raise NotImplementedError
    
    def __repr__(self) -> str: return f"Grover for a search with {self.size} q-bits" + (f" {self.__success}" if self.__success != None else "")

class ToffoliSearch(QuantumSearch):
    def __init__(self, size):
        super().__init__(size, lambda *p : all(p))

    @staticmethod
    def random(size:int):
        return ToffoliSearch(size)
    
    def grover_quantum(self, iterations: int):
        q = QuantumRegister(self.size, "q")
        flag = AncillaRegister(1, "flag")

        state_preparation = QuantumCircuit(q, flag)
        state_preparation.h(q)
        state_preparation.x(flag)

        oracle = QuantumCircuit(q, flag)
        oracle.mcp(pi, q, flag)

        operator = GroverOperator(oracle, mcx_mode="noancilla")

        q2 = QuantumRegister(self.size, "q")
        grover = QuantumCircuit(q2, flag, name="grover")
        grover.compose(state_preparation, inplace=True)

        grover.compose(operator.power(iterations), inplace=True)
        grover.measure_all()

        return grover
    
    def __repr__(self) -> str: return f"Grover over Toffoli with {self.size} q-bits"

def generate_valid_example(size: int) -> ToffoliSearch:
    print(f"%%% [GEN] Generating a valid Toffoli Search with {size} nodes...")
    return ToffoliSearch.random(size)

def gen_header(csv_writer, ttype):
    if ttype in ("clue", "ddsim"):
        csv_writer.writerow(["size", "red. ratio", "time_lumping", "memory (MB)", "gate"])
    elif ttype == "full_clue":
        csv_writer.writerow(["size", "time_lumping", "kappa", "time_iteration", "memory (MB)", "gate"])
    elif ttype == "full_ddsim":
        csv_writer.writerow(["size", "kappa", "time_iteration", "memory (MB)", "gate"])
    else:
        raise NotImplementedError(f"Type of file {ttype} not recognized")

def clue_reduction(size: int, result_file, timeout=0): 
    r'''
        This method computes the CLUE lumping

        This method generates a Toffoli Grover gate and computes the lumping of the problem matrix associated
        to it with CLUE.
         
        It stores the execution time of the lumping plus the reduction ratio. It stores the result on ``result_file``.
        ["size", "red. ratio", "time_lumping", "memory", "gate"]
    '''
    grover = generate_valid_example(size)

    print(f"%%% [clue] Creating the full system with the problem matrix")
    system = FODESystem.LinearSystem(grover.grover_matrix(), lumping_subspace=NumericalSubspace)
    obs = tuple([SparseVector.from_list(system.size*[1], field=CC)])
    
    print(f"%%% [clue] Computing the lumped system...")
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            ctime = time()
            lumped = system.lumping(obs, print_reduction=False, print_system=False)
            ctime = time()-ctime
    except TimeoutError:
        print(f"%%% [clue] Timeout reached for execution")
        ctime = Inf
    memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(f"%%% [clue] Storing the data...")
    if 2 != lumped.size:
        print(f"%%% [clue] ERROR!! Found weird dimension in lumping -- \n%%% \t* Expected: {2}\n%%% \t* Got: {lumped.size}\n%%% \t* Search: {grover}")
    result_file.writerow([size, lumped.size/system.size, ctime, memory, repr(grover)])

def ddsim_reduction(size: int, result_file, timeout=0): 
    r'''
        This method computes the DDSIM lumping

        This method generates a Toffoli Grover gate and computes the lumping of the problem matrix associated
        to it with DDSIM. Since we can not use linear algebra on Decision Diagrams, we compute the simulation time
        of a quantum circuit that repeats the circuit as many times as the bound of the lumping size (i.e., 2).
         
        It stores the execution time of the lumping plus the reduction ratio. It stores the result on ``result_file``.
        ["size", "red. ratio", "time_lumping", "memory", "gate"]
    '''
    grover = generate_valid_example(size)

    print(f"%%% [ddsim] Creating the full circuit and job to simulate with DDSIM")
    U_P = grover.grover_quantum(2)
    backend = ddsim.DDSIMProvider().get_backend("qasm_simulator")
    
    print(f"%%% [ddsim] Computing the simulation of the circuit...")
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            ctime = time()
            ## Executing the circuit one time
            job = execute(U_P, backend, shots=1)
            job.result()
            ctime = time()-ctime
    except TimeoutError:
        print(f"%%% [ddsim] Timeout reached for execution")
        ctime = Inf
    memory = tracemalloc.get_traced_memory()[1] / (2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(f"%%% [ddsim] Storing the data...")
    result_file.writerow([size, 2/2**size, ctime, memory, repr(grover)])

def clue_iteration(size: int, iterations, result_file, timeout=0): 
    r'''
        This method computes the CLUE iteration

        This method generates a Toffoli Grover gate and computes the lumping of the problem matrix associated
        to it with CLUE. 
        
        Then it computes the application of the reduced matrix up to ``iterations`` times. 
         
        It stores the execution time of the lumping, the number of iterations and the time for the computed iteration. 
        It stores the result on ``result_file``.
        ["size", "time_lumping", "kappa", "time_iteration", "memory", "gate"]
    '''
    grover = generate_valid_example(size)

    print(f"%%% [full-clue] Creating the full system with the problem matrix")
    system = FODESystem.LinearSystem(grover.grover_matrix(), lumping_subspace=NumericalSubspace)
    obs = tuple([SparseVector.from_list(system.size*[1], field=CC)])
    
    print(f"%%% [full-clue] Computing the lumped system...")
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            lump_time = time()
            ## Executing the circuit one time
            lumped = system.lumping(obs, print_reduction=False, print_system=False)
            lump_time = time()-lump_time
    except TimeoutError:
        print(f"%%% [full-clue] Timeout reached for execution")
        lump_time = Inf
        memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
        tracemalloc.stop()
        result_file.writerow([size, lump_time, iterations, Inf, memory, repr(grover)])
        return

    print(f"%%% [full-clue] Checking size...")
    if 2 != lumped.size:
        print(f"%%% [full-clue] ERROR!! Found weird dimension in lumping -- \n%%% \t* Expected: {2}\n%%% \t* Got: {lumped.size}\n%%% \t* Gate: {grover}")

    print(f"%%% [full-clue] Computing the iteration (U)^iterations")
    it_time = time()
    _ = matrix_power(lumped.construct_matrices("polynomial")[0].to_numpy(dtype=cdouble), iterations)
    it_time = time() - it_time
    memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
    tracemalloc.stop()

    ## We check if the matrix is diagonal
    result_file.writerow([size, lump_time, iterations, it_time, memory, repr(grover)])

def ddsim_iteration(size: int, iterations, result_file, timeout=0): 
    r'''
        This method computes the DDSIM iteration

        This method generates a Toffoli Grover gate and computes application of ``iterations`` times the Grover search gate.
         
        It stores the number of iterations and the time for the computed iteration. 
        It stores the result on ``result_file``.
        ["size", "kappa", "time_iteration", "memory", "gate"]
    '''
    grover = generate_valid_example(size)

    print(f"%%% [ddsim] Creating the full circuit and job to simulate with DDSIM")
    circuit = grover.grover_quantum(iterations)
    backend = ddsim.DDSIMProvider().get_backend("qasm_simulator")
    
    print(f"%%% [ddsim] Computing the simulation of the circuit...")
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            ctime = time()
            ## Executing the circuit one time
            job = execute(circuit, backend, shots=1)
            job.result()
            ctime = time()-ctime
    except TimeoutError:
        print(f"%%% [ddsim] Timeout reached for execution")
        ctime = Inf
    memory = tracemalloc.get_traced_memory()[1] / (2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(f"%%% [ddsim] Storing the data...")
    result_file.writerow([size, iterations, ctime, memory, repr(grover)])

if __name__ == "__main__":
    n = 1; m = 5; M=10; ttype="clue"; repeats=100; timeout=0
    ## Processing arguments
    while n < len(sys.argv):
        if sys.argv[n].startswith("-"):
            if sys.argv[n].endswith("m"):
                m = int(sys.argv[n+1]); n+=2
            elif sys.argv[n].endswith("M"):
                M = int(sys.argv[n+1]); n+=2
            elif sys.argv[n].endswith("t"):
                ttype = sys.argv[n+1] if sys.argv[n+1] in ("clue", "ddsim", "full_clue", "full_ddsim") else ttype
                n += 2
            elif sys.argv[n].endswith("to"):
                timeout = int(sys.argv[n+1]); n+=2
            elif sys.argv[n].endswith("r"):
                repeats = int(sys.argv[n+1]); n+=2
        else:
            n += 1

    methods = [clue_reduction, ddsim_reduction, clue_iteration, ddsim_iteration]
    method = methods[["clue", "ddsim", "full_clue", "full_ddsim"].index(ttype)]
    existed = os.path.exists(os.path.join(SCRIPT_DIR, f"[result]q_search_{ttype}.csv"))
    with open(os.path.join(SCRIPT_DIR, f"[result]q_search_{ttype}.csv"), "at" if existed else "wt") as result_file:
        csv_writer = writer(result_file)
        if not existed:
            gen_header(csv_writer, ttype)
        print(f"##################################################################################")
        print(f"### EXECUTION ON GROVER [{m=}, {M=}, {repeats=}, method={ttype}]")
        print(f"##################################################################################")
        for size in range(m, M+1):
            for execution in range(1,repeats+1):
                print(f"### Starting execution {execution}/{repeats} ({size=})")
                if ttype in ("clue", "ddsim"):
                    method(size, csv_writer, timeout=timeout)
                else:
                    #for it in (1,10,100):#,1000):#,10000)
                    it = ceil(sqrt(2**size))
                    print(f"------ Case with {it} iterations")
                    method(size, it, csv_writer)
                print(f"### Finished execution {execution}/{repeats}")
                result_file.flush()