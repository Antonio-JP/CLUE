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
from functools import lru_cache
from itertools import product
from math import ceil, gcd, log2, sqrt
from mqt import ddsim #pylint: disable=no-name-in-module
from numpy import cdouble, pi, Inf
from numpy.linalg import matrix_power
from qiskit import AncillaRegister, QuantumCircuit, QuantumRegister, execute
from qiskit.circuit.library import GroverOperator
from random import randint, choice
from time import time, process_time

## Imports from the local folder
from misc import *

def __is_prime(n):
    return all(n%i != 0 for i in primes(int(ceil(sqrt(n)+1))))

__CACHED_PRIMES = [2,3,5,7,11,13]
@lru_cache(maxsize=None)
def primes(bound):
    if bound > __CACHED_PRIMES[-1]:
        for n in range(__CACHED_PRIMES[-1]+1, bound):
            if __is_prime(n): __CACHED_PRIMES.append(n)
            
    li = 0; ui = len(__CACHED_PRIMES) - 1
    if __CACHED_PRIMES[-1] <= bound: return __CACHED_PRIMES[:len(__CACHED_PRIMES)]
    while (ui-li) > 1:
        ci = (ui+li)//2
        if __CACHED_PRIMES[ci] < bound: li = ci
        elif __CACHED_PRIMES[ci] > bound: ui = ci
        else: ui = ci; break # found the element
    return __CACHED_PRIMES[:ui]

class QuantumOrder:
    def __init__(self, N, x, size=None):
        if gcd(N, x) != 1:
            raise ValueError(f"The mod number {N=} and the inner number {x=} must be coprime")
        self.__N = N
        self.__x = x
        if size == None: size = ceil(log2(self.N))
        self.__size = size

    @staticmethod
    def random(size) -> QuantumOrder:
        print(f"[random @ {size}] Trying to get the product of two different primes")
        candidates = [(p,q) for p,q in product(primes(ceil(2**size / 3))[1:], repeat=2) if (p < q and p*q >= 2**(size-1) and p*q < 2**size)]
        if len(candidates) == 0:
            raise ValueError(f"Impossible thing happened: no composed number between {2**(size-1)} and {2**size}")
        p,q = choice(candidates)
        N = p*q
                    
        print(f"[random @ {size}] {N} = {p} * {q}")
        ## We generate a coprime number ´x´
        x = randint(2, N-1)
        while gcd(N, x) != 1: x = randint(2, N-1)

        return QuantumOrder(N, x, size)

    @property
    def N(self) -> int: return self.__N
    @property
    def x(self) -> int: return self.__x
    @property
    def size(self) -> int: return self.__size

    def order_matrix(self): 
        f = [(self.x*y)%self.N if y < self.N else y for y in range(2**self.size)]

        M = SparseRowMatrix(2**self.size, CC)
        for i in range(2**self.size):
            M.increment(f[i], i, 1)

        return M
    
    def order_quantum(self):
        raise NotImplementedError
    
    def __repr__(self) -> str: return f"Order finding in {self.size} q-bits: mult. by {self.x} modulo {self.N}" 

def generate_valid_example(size: int) -> QuantumOrder:
    print(f"%%% [GEN] Generating a valid Order finding instance...")
    return QuantumOrder.random(size)

def gen_header(csv_writer, ttype):
    if ttype in ("clue", "ddsim"):
        csv_writer.writerow(["size", "N", "red. ratio", "time_lumping", "memory (MB)", "problem"])
    elif ttype == "full_clue":
        csv_writer.writerow(["size", "N", "time_lumping", "kappa", "time_iteration", "memory (MB)", "problem"])
    elif ttype == "full_ddsim":
        csv_writer.writerow(["size", "N", "kappa", "time_iteration", "memory (MB)", "problem"])
    else:
        raise NotImplementedError(f"Type of file {ttype} not recognized")

def clue_reduction(size: int, result_file, timeout=0): 
    r'''
        This method computes the CLUE lumping

        This method generates a Order finding problem and computes the lumping of the problem matrix associated
        to it with CLUE.
         
        It stores the execution time of the lumping plus the reduction ratio. It stores the result on ``result_file``.
        ["size", "N", "red. ratio", "time_lumping", "memory", "problem"]
    '''
    order = generate_valid_example(size)

    print(f"%%% [clue] Creating the full system with the problem matrix")
    system = FODESystem.LinearSystem(order.order_matrix(), lumping_subspace=NumericalSubspace)
    obs = tuple([SparseVector.from_list([0,1] + (system.size-2)*[0], field=CC)])
    
    print(f"%%% [clue] Computing the lumped system...")
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            ctime = process_time()
            lumped = system.lumping(obs, print_reduction=False, print_system=False)
            ctime = process_time()-ctime
    except TimeoutError:
        print(f"%%% [clue] Timeout reached for execution")
        ctime = Inf
    memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(f"%%% [clue] Storing the data...")
    result_file.writerow([size, order.N, lumped.size/system.size, ctime, memory, repr(order)])

def ddsim_reduction(size: int, result_file, timeout=0): 
    r'''
        This method computes the DDSIM lumping

        This method generates a order finding problem and computes the lumping of the problem matrix associated
        to it with DDSIM. Since we can not use linear algebra on Decision Diagrams, we compute the simulation time
        of a quantum circuit that repeats the circuit as many times as the bound of the lumping size (i.e., 2).
         
        It stores the execution time of the lumping plus the reduction ratio. It stores the result on ``result_file``.
        ["size", "red. ratio", "time_lumping", "memory", "gate"]
    '''
    raise NotImplementedError("DDSIM for ORder finding not implemented")

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
    order = generate_valid_example(size)

    print(f"%%% [full-clue] Creating the full system with the problem matrix")
    system = FODESystem.LinearSystem(order.order_matrix(), lumping_subspace=NumericalSubspace)
    obs = tuple([SparseVector.from_list([0,1] + (system.size-2)*[0], field=CC)])
    
    print(f"%%% [full-clue] Computing the lumped system...")
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            lump_time = process_time()
            ## Executing the circuit one time
            lumped = system.lumping(obs, print_reduction=False, print_system=False)
            lump_time = process_time()-lump_time
    except TimeoutError:
        print(f"%%% [full-clue] Timeout reached for execution")
        lump_time = Inf
        memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
        tracemalloc.stop()
        result_file.writerow([size, order.N, lump_time, iterations, Inf, memory, repr(order)])
        return

    print(f"%%% [full-clue] Computing the iteration (U)^iterations")
    it_time = process_time()
    _ = matrix_power(lumped.construct_matrices("polynomial")[0].to_numpy(dtype=cdouble), iterations)
    it_time = process_time() - it_time
    memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
    tracemalloc.stop()

    ## We check if the matrix is diagonal
    result_file.writerow([size, order.N, lump_time, iterations, it_time, memory, repr(order)])

def ddsim_iteration(size: int, iterations, result_file, timeout=0): 
    r'''
        This method computes the DDSIM iteration

        This method generates a Toffoli Grover gate and computes application of ``iterations`` times the Grover search gate.
         
        It stores the number of iterations and the time for the computed iteration. 
        It stores the result on ``result_file``.
        ["size", "kappa", "time_iteration", "memory", "gate"]
    '''
    raise NotImplementedError("DDSIM for Order finding not implemented")

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
    existed = os.path.exists(os.path.join(SCRIPT_DIR, f"[result]q_order_{ttype}.csv"))
    with open(os.path.join(SCRIPT_DIR, f"[result]q_order_{ttype}.csv"), "at" if existed else "wt") as result_file:
        csv_writer = writer(result_file)
        if not existed:
            gen_header(csv_writer, ttype)
        print(f"##################################################################################")
        print(f"### EXECUTION ON ORDER [{m=}, {M=}, {repeats=}, method={ttype}]")
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