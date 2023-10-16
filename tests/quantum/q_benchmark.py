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
from math import ceil, sqrt
from mqt import ddsim #pylint: disable=no-name-in-module
from mqt.bench.benchmarks import (ae, dj, ghz, graphstate, pricingput, pricingcall, portfolioqaoa, portfoliovqe, qft, 
                                  qpeexact, qpeinexact, qwalk, tsp, qnn, vqe, wstate)
from numpy import asarray, cdouble, ndarray, Inf
from numpy.linalg import matrix_power
from qiskit import Aer, QuantumCircuit, execute
from time import process_time

## Imports from the local folder
from misc import *

VALID_BENCHMARKS = {"ae": ae, "dj" : dj, "ghz": ghz, "graphstate": graphstate, "pricingput": pricingput, "pricingcall": pricingcall, 
                    "portfolioqaoa": portfolioqaoa, "portfoliovqe": portfoliovqe, "qft":qft, "qpeexact": qpeexact, "qpeinexact": qpeinexact,
                    "qwalk": qwalk, "tsp": tsp, "qnn": qnn, "vqe": vqe, "wstate": wstate}
FULL_NAMES = {"ae": "Amplitude Estimation", "dj" : "Deitsch-Jozsa", "ghz": "Greenberger-Horne_zeilinger", "graphstate": "Graph State", 
              "pricingput": "Pricing Put Option", "pricingcall": "Princing Put Call", "portfolioqaoa": "Portfolio Optimization", "portfoliovqe": "VQE Portfolio Optimizaiton", 
              "qft":"Quantum Fourier Transform", "qpeexact": "Exact Quantum Phase Estimation", "qpeinexact": "Inexact Qunatum Phase Estimation", "qwalk": "Quantum Walk", 
              "tsp": "Travelling Salesman", "qnn": "Quantum Neural Network", "vqe": "Variational Quantum Eigensolver ", "wstate": "W State"}

class QuantumBenchmark:
    BACKEND = Aer.get_backend("unitary_simulator")
    DDSIM = ddsim.DDSIMProvider().get_backend("qasm_simulator")

    def __init__(self, name, size):
        if not name in VALID_BENCHMARKS:
            raise TypeError(f"The given benchmark ({name}) not recognized.")
        self.__name = name; self.__size = size
        ## We try to build the circuit
        self.__circuit = VALID_BENCHMARKS[name].create_circuit(size)
        ## Cache for other derived attributes
        self.__unitary = None

    @staticmethod
    def random(name, size) -> QuantumBenchmark:
        r'''There is no randomness in this class'''
        return QuantumBenchmark(name, size)

    @property
    def name(self) -> str: return self.__name
    @property
    def full_name(self) -> str: return FULL_NAMES[self.name]
    @property
    def size(self) -> int: return self.__size
    @property
    def circuit(self) -> QuantumCircuit:return self.__circuit
    @property
    def unitary(self) -> ndarray:
        if self.__unitary is None:
            no_measured = self.circuit.remove_final_measurements(False)
            job = execute(no_measured, QuantumBenchmark.BACKEND, shots=8192)
            unitary = job.result().get_unitary(no_measured)
            self.__unitary = asarray(unitary)
        return self.__unitary
    
    def unitary_matrix(self) -> SparseRowMatrix: 
        return SparseRowMatrix.from_vectors([SparseVector.from_list(row, CC) for row in self.unitary])
    
    def quantum_ddsim(self) -> QuantumCircuit:
        raise NotImplementedError
    
    def __repr__(self) -> str: return f"{self.name} with {self.size} qubits." 

def generate_valid_example(name: str, size: int) -> QuantumBenchmark:
    return QuantumBenchmark.random(name, size)

def gen_header(csv_writer, ttype):
    if ttype in ("clue", "ddsim"):
        csv_writer.writerow(["size", "name", "obs", "red. ratio", "time_lumping", "memory (MB)", "problem"])
    elif ttype == "full_clue":
        csv_writer.writerow(["size", "name", "obs", "time_lumping", "kappa", "time_iteration", "memory (MB)", "problem"])
    elif ttype == "full_ddsim":
        csv_writer.writerow(["size", "name", "obs", "kappa", "time_iteration", "memory (MB)", "problem"])
    else:
        raise NotImplementedError(f"Type of file {ttype} not recognized")

def clue_reduction(name: str, size: int, result_file, observable: int | str = 0, timeout=0) -> int: 
    r'''
        This method computes the CLUE lumping

        This method generates a problem with given name and size and compute hte CLUE lumping of its unitary matrix.
         
        It stores the execution time of the lumping plus the reduction ratio. It stores the result on ``result_file``.
        ["size", "name", "obs", "red. ratio", "time_lumping", "memory", "problem"]
    '''
    benchmark = generate_valid_example(name, size)
    size = sum(len(reg) for reg in benchmark.circuit.qregs) # adjusting just in case
    if isinstance(observable, int) and observable >= 0 and observable < 2**size:
        list_to_obs = (2**size)*[0]; list_to_obs[observable] = 1
        obs = SparseVector.from_list(list_to_obs)
    elif isinstance(observable, str) and "H" in observable:
        obs = SparseVector.from_list((2**size)*[1])
    else:
        raise ValueError(f"%%% [clue @ {name}] The observable (given {observable=}) must be ain integer between 0 and {2**size-1} or a string containing 'H'")

    print(f"%%% [clue @ {name}] Creating the full system to be lumped...")
    system = FODESystem.LinearSystem(benchmark.unitary_matrix(), lumping_subspace=NumericalSubspace)
    obs = tuple([obs])
    
    print(f"%%% [clue @ {name}] Computing the lumped system...")
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            ctime = process_time()
            lumped = system.lumping(obs, print_reduction=False, print_system=False)
            ctime = process_time()-ctime
    except TimeoutError:
        print(f"%%% [clue @ {name}] Timeout reached for execution")
        ctime = Inf
    memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(f"%%% [clue @ {name}] Storing the data...")
    result_file.writerow([size, benchmark.name, observable, "unknown" if ctime == Inf else lumped.size/system.size, ctime, memory, repr(benchmark)])

    return ctime

def ddsim_reduction(name:str, size: int, result_file, observable: int|str = 0, timeout=0) -> int: 
    r'''
        This method computes the DDSIM lumping

        This method generates a problem from a benchmark and computes the lumping of the problem matrix associated
        to it with DDSIM. 

        It stores the execution time of the lumping plus the reduction ratio. It stores the result on ``result_file``.
        ["size", "name", "red. ratio", "time_lumping", "memory", "problem"]
    '''
    raise NotImplementedError("DDSIM for Benchmarks not implemented")

def clue_iteration(name: str, size: int, iterations, result_file, observable: int | str = 0, timeout=0) -> int: 
    r'''
        This method computes the CLUE iteration

        This method generates a problem with given name and size and compute hte CLUE lumping of its unitary matrix.
        
        Then it computes the application of the reduced matrix up to ``iterations`` times. 
         
        It stores the execution time of the lumping, the number of iterations and the time for the computed iteration. 
        It stores the result on ``result_file``.
        ["size", "name", "time_lumping", "kappa", "time_iteration", "memory", "gate"]
    '''
    benchmark = generate_valid_example(name, size)
    size = len(benchmark.circuit.qregs[0]) # adjusting just in case
    if isinstance(observable, int) and observable >= 0 and observable < 2**size:
        list_to_obs = (2**size)*[0]; list_to_obs[observable] = 1
        obs = SparseVector.from_list(list_to_obs)
    elif isinstance(observable, str) and "H" in observable:
        obs = SparseVector.from_list((2**size)*[1])
    else:
        raise ValueError(f"%%% [clue] The observable (given {observable=}) must be ain integer between 0 and {2**size-1} or a string containing 'H'")

    print(f"%%% [full-clue @ {name}] Creating the full system to be lumped...")
    system = FODESystem.LinearSystem(benchmark.unitary_matrix(), lumping_subspace=NumericalSubspace)
    obs = tuple([obs])
    
    print(f"%%% [full-clue @ {name}] Computing the lumped system...")
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            lump_time = process_time()
            ## Executing the circuit one time
            lumped = system.lumping(obs, print_reduction=False, print_system=False)
            lump_time = process_time()-lump_time
    except TimeoutError:
        print(f"%%% [full-clue @ {name}] Timeout reached for execution")
        lump_time = Inf
        memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
        tracemalloc.stop()
        result_file.writerow([size, benchmark.name, observable, lump_time, iterations, Inf, memory, repr(benchmark)])
        return

    print(f"%%% [full-clue @ {name}] Computing the iteration (U)^iterations")
    it_time = process_time()
    _ = matrix_power(lumped.construct_matrices("polynomial")[0].to_numpy(dtype=cdouble), iterations)
    it_time = process_time() - it_time
    memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
    tracemalloc.stop()

    ## We check if the matrix is diagonal
    result_file.writerow([size, benchmark.name, observable, lump_time, iterations, it_time, memory, repr(benchmark)])

    return lump_time + it_time

def ddsim_iteration(name: str, size: int, iterations, result_file, observable: int | str = 0, timeout=0) -> int: 
    r'''
        This method computes the DDSIM iteration

        This method generates a benchmark problem and computes application of ``iterations`` times the its circuit.
         
        It stores the number of iterations and the time for the computed iteration. 
        It stores the result on ``result_file``.
        ["size", "name", "kappa", "time_iteration", "memory", "gate"]
    '''
    raise NotImplementedError("DDSIM for Benchmarks not implemented")

if __name__ == "__main__":
    n = 1; m = 5; M=10; ttype="clue"; repeats=100; timeout=None; name=None; obs=list()
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
            elif sys.argv[n].endswith("n"):
                name = sys.argv[n+1]; n+=2
            elif sys.argv[n].endswith("obs"):
                try:
                    observable = int(sys.argv[n+1])
                except:
                    observable = sys.argv[n+1]
                obs.append(observable); n += 2
        else:
            n += 1

    if name is None:
        raise ValueError(f"At least a name must be given for the tests (with command -n)")
    methods = [clue_reduction, ddsim_reduction, clue_iteration, ddsim_iteration]
    method = methods[["clue", "ddsim", "full_clue", "full_ddsim"].index(ttype)]
    existed = os.path.exists(os.path.join(SCRIPT_DIR, "results", f"[result]q_benchmark_{ttype}.csv"))
    with open(os.path.join(SCRIPT_DIR, "results", f"[result]q_benchmark_{ttype}.csv"), "at" if existed else "wt") as result_file:
        csv_writer = writer(result_file)
        if not existed:
            gen_header(csv_writer, ttype)
        print(f"##################################################################################")
        print(f"### EXECUTION ON BENCHMARK {name} [{m=}, {M=}, {repeats=}, method={ttype}]")
        print(f"##################################################################################")
        for size in range(m, M+1):
            for execution in range(1,repeats+1):
                my_obs = list(range(2**size)) + ["H"] if len(obs) == 0 else obs
                print(f"### Starting execution {execution}/{repeats} ({size=})")
                for i,observable in enumerate(my_obs):
                    if ttype in ("clue", "ddsim"):
                        timeout -= method(name, size, csv_writer, observable=observable, timeout=timeout if timeout != None else 0)
                    else:
                        #for it in (1,10,100):#,1000):#,10000)
                        it = ceil(sqrt(2**size))
                        print(f"------ Case with {it} iterations")
                        timeout -= method(name, size, it, csv_writer, observable=observable, timeout=timeout if timeout != None else 0)
                    print(f"### -- Finished execution with {observable=} ({i+1}/{len(my_obs)})")
                    result_file.flush()
                    if timeout != None and timeout <= 0:
                        break
                    elif timeout != None:
                        timeout = ceil(timeout)
                print(f"### Finished execution {execution}/{repeats}")