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
FULL_NAMES = {"ae": "Amplitude Estimation", "dj" : "Deutsch-Jozsa", "ghz": "Greenberger-Horne-Zeilinger", "graphstate": "Graph State", 
              "pricingput": "Pricing Put Option", "pricingcall": "Princing Call Option", "portfolioqaoa": "Portfolio Optimization", "portfoliovqe": "VQE Portfolio Optimizaiton", 
              "qft":"Quantum Fourier Transform", "qpeexact": "Exact Quantum Phase Estimation", "qpeinexact": "Inexact Qunatum Phase Estimation", "qwalk": "Quantum Walk", 
              "tsp": "Travelling Salesman", "qnn": "Quantum Neural Network", "vqe": "Variational Quantum Eigensolver ", "wstate": "W State"}

class QuantumBenchmark(Experiment):
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
    def qbits(self) -> int: return self.__size
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
    
    def __repr__(self) -> str: return f"{self.name} with {self.qbits} qubits." 

    ## NECESSARY METHOD FOR EXPERIMENT
    def size(self) -> int: return sum(len(reg) for reg in self.circuit.qregs)
    def correct_size(self) -> int: return None
    def matrix(self) -> SparseRowMatrix: return self.unitary_matrix()
    def quantum(self) -> tuple[QuantumCircuit, Parameter]: return self.circuit, None

def generate_example(name: str, size: int, obs) -> QuantumBenchmark:
    return QuantumBenchmark.random(name, size)

def generate_header(csv_writer, ttype):
    if ttype in ("clue", "ddsim"):
        csv_writer.writerow(["size", "name", "obs", "red. ratio", "time_lumping", "memory (MB)", "problem"])
    elif ttype == "full_clue":
        csv_writer.writerow(["size", "name", "obs", "time_lumping", "kappa", "time_iteration", "memory (MB)", "problem"])
    elif ttype == "full_ddsim":
        csv_writer.writerow(["size", "name", "obs", "kappa", "time_iteration", "memory (MB)", "problem"])
    else:
        raise NotImplementedError(f"Type of file {ttype} not recognized")

## METHODS TO GENERATE OBSERVABLES
def generate_observable_clue(example: QuantumBenchmark, size: int, obs) -> tuple[SparseVector]:
    size = example.size()
    if isinstance(obs, int) and obs >= 0 and obs < 2**size:
        list_to_obs = (2**size)*[0]; list_to_obs[obs] = 1
        obs = SparseVector.from_list(list_to_obs)
    elif isinstance(obs, str) and "H" in obs:
        obs = SparseVector.from_list((2**size)*[1])
    else:
        raise ValueError(f"%%% [clue @ {name}] The observable (given {obs=}) must be ain integer between 0 and {2**size-1} or a string containing 'H'")
    return tuple([obs])

def generate_observable_ddsim(_: QuantumBenchmark, size: int, obs) -> bool:
    return obs == "H"

### METHODS TO GENERATE THE DATA
def generate_data(example: QuantumBenchmark, *args) -> list:
    return [size, example.name] + [*args] + [repr(example)]

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
                ttype = sys.argv[n+1] if sys.argv[n+1] in ("clue", "ddsim", "direct", "full_clue", "full_ddsim", "full_direct") else ttype
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
    methods = [clue_reduction, ddsim_reduction, direct_reduction, clue_iteration, ddsim_iteration, direct_iteration]
    method = methods[["clue", "ddsim", "direct", "full_clue", "full_ddsim", "full_direct"].index(ttype)]
    existed = os.path.exists(os.path.join(SCRIPT_DIR, "results", f"[result]q_benchmark_{ttype}.csv"))
    with open(os.path.join(SCRIPT_DIR, "results", f"[result]q_benchmark_{ttype}.csv"), "at" if existed else "wt") as result_file:
        csv_writer = writer(result_file)
        if not existed:
            generate_header(csv_writer, ttype)
        print(f"##################################################################################")
        print(f"### EXECUTION ON BENCHMARK {name} [{m=}, {M=}, {repeats=}, method={ttype}]")
        print(f"##################################################################################")
        for size in range(m, M+1):
            benchmark = generate_example(name, size, 0)
            b_size = len(benchmark.circuit.qregs[0]) # adjusting just in case
            for execution in range(1,repeats+1):
                my_obs = ([0] + ["H"] + list(range(1,2**b_size))) if len(obs) == 0 else obs
                print(f"### Starting execution {execution}/{repeats} ({size=})")
                rem_timeout = timeout
                for i,observable in enumerate(my_obs):
                    if ttype in ("clue", "ddsim"):
                        rem_timeout -= method(
                            name, 
                            generate_example, 
                            generate_observable_clue if ttype != "ddsim" else generate_observable_ddsim,
                            generate_data,
                            csv_writer, 
                            size, observable, #args
                            timeout=rem_timeout
                        )
                    else:
                        for it in (1,ceil(sqrt(2**size))):#,1000):#,10000)
                            print(f"------ Case with {it} iterations")
                            rem_timeout -= method(
                                name, 
                                generate_example, 
                                generate_observable_clue if ttype != "full_ddsim" else generate_observable_ddsim,
                                generate_data,
                                it,
                                csv_writer, 
                                size, observable, #args
                                timeout=ceil(rem_timeout)
                            )
                    print(f"### -- Finished execution with {observable=} ({i+1}/{len(my_obs)})")
                    result_file.flush()
                    if rem_timeout != None and rem_timeout <= 0:
                        break
                    elif rem_timeout != None:
                        rem_timeout = ceil(rem_timeout)
                print(f"### Finished execution {execution}/{repeats}")