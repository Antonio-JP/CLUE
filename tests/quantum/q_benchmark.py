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
from mqt.bench.benchmarks import (ae, dj, ghz, graphstate, pricingput, pricingcall, portfolioqaoa, portfoliovqe, qft, 
                                  qpeexact, qpeinexact, qwalk, tsp, qnn, vqe, wstate)
from numpy import asarray, ndarray
from qiskit import Aer, QuantumCircuit, execute

## Imports from the local folder
from misc import *

VALID_BENCHMARKS = {"ae": ae, "dj" : dj, "ghz": ghz, "graphstate": graphstate, "hhl": None, "pricingput": pricingput, "pricingcall": pricingcall, 
                    "portfolioqaoa": portfolioqaoa, "portfoliovqe": portfoliovqe, "qft":qft, "qpeexact": qpeexact, "qpeinexact": qpeinexact,
                    "qwalk": qwalk, "tsp": tsp, "qnn": qnn, "vqe": vqe, "wstate": wstate}
FULL_NAMES = {"ae": "Amplitude Estimation", "dj" : "Deutsch-Jozsa", "ghz": "Greenberger-Horne-Zeilinger", "graphstate": "Graph State", "hhl": "HHL Algorithm",
              "pricingput": "Pricing Put Option", "pricingcall": "Princing Call Option", "portfolioqaoa": "Portfolio Optimization", "portfoliovqe": "VQE Portfolio Optimizaiton", 
              "qft":"Quantum Fourier Transform", "qpeexact": "Exact Quantum Phase Estimation", "qpeinexact": "Inexact Quantum Phase Estimation", "qwalk": "Quantum Walk", 
              "tsp": "Travelling Salesman", "qnn": "Quantum Neural Network", "vqe": "Variational Quantum Eigensolver ", "wstate": "W State"}

class QuantumBenchmark(Experiment):
    BACKEND = Aer.get_backend("unitary_simulator")

    def __init__(self, name, size):
        if not name in VALID_BENCHMARKS:
            raise TypeError(f"The given benchmark ({name}) not recognized.")
        self.__name = name; self.__size = size
        ## We try to build the circuit
        if name == "hhl": # special case for HHL
            self.__circuit = QuantumCircuit.from_qasm_file(f"./circuits/hhl_indep_qiskit_{3*size-1}.qasm")
        else:
            self.__circuit = VALID_BENCHMARKS[name].create_circuit(size)
        ## Cache for other derived attributes
        self.__unitary = None
        self.observable = None

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
    
    def __repr__(self) -> str: return f"{self.name} with {self.qbits} qubits." 

    ## NECESSARY METHOD FOR EXPERIMENT
    def size(self) -> int: return self.circuit.num_qubits
    def correct_size(self) -> int: return None
    def matrix(self) -> SparseRowMatrix: return self.unitary_matrix()
    def quantum(self) -> tuple[QuantumCircuit, Parameter]: return self.circuit, None
    def data(self): return [self.full_name, self.observable]

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
        example.observable = obs
        list_to_obs = (2**size)*[0]; list_to_obs[obs] = 1
        obs = SparseVector.from_list(list_to_obs)
    elif isinstance(obs, str) and "H" in obs:
        example.observable = obs
        obs = SparseVector.from_list((2**size)*[1])
    else:
        raise ValueError(f"%%% [clue @ {example.name}] The observable (given {obs=}) must be an integer between 0 and {2**size-1} or a string containing 'H'")
    return tuple([obs])

def generate_observable_ddsim(example: QuantumBenchmark, size: int, obs) -> bool:
    size = example.size()
    if isinstance(obs, int) and obs >= 0 and obs < 2**size:
        example.observable = obs
    elif isinstance(obs, str) and "H" in obs:
        example.observable = obs
        obs = True
    else:
        raise ValueError(f"%%% [ddsim @ {example.name}] The observable (given {obs=}) must be an integer between 0 and {2**size-1} or a string containing 'H'")
    return obs

### SCRIPT ARGUMENT PROCCESS METHOD
def get_benchmark_args(*argv) -> tuple[str, list[Any]]:
    name = None; obs = list()
    for (i, arg) in enumerate(argv):
        if arg == "-n": # name
            name = argv[i+1]
        elif arg == "-obs":
            try:
                observable = int(argv[i+1])
            except:
                observable = argv[i+1]
            obs.append(observable)
    
    ## Checking correctness of arguments
    if name is None:
        raise ValueError(f"Script argument [-n] not well used: we required a follow-up name.")

    return name, obs

if __name__ == "__main__":
    ## Processing the arguments
    ## Generic part
    ttype, script = get_method(*sys.argv)
    m, M = get_size_bounds(*sys.argv)
    timeout = get_timeout(*sys.argv)
    repeats = get_repeats(*sys.argv)
    ## Specific part
    name, obs = get_benchmark_args(*sys.argv)    

    main_script(SCRIPT_DIR, SCRIPT_NAME, name,
                [generate_header, generate_example, generate_observable_clue, generate_observable_ddsim],
                ttype, script,               
                m, M,                        
                timeout, repeats,
                obs)