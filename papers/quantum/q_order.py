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
from functools import lru_cache
from itertools import product
from math import ceil, gcd, log2, sqrt
from random import randint, choice

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

class QuantumOrder(Experiment):
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
    def qbits(self) -> int: return self.__size

    def order_matrix(self): 
        f = [(self.x*y)%self.N if y < self.N else y for y in range(2**self.size())]

        M = SparseRowMatrix(2**self.size(), CC)
        for i in range(2**self.size()):
            M.increment(f[i], i, 1)

        return M
    
    def order_quantum(self):
        raise NotImplementedError
    
    def __repr__(self) -> str: return f"Order finding in {self.size()} q-bits: mult. by {self.x} modulo {self.N}" 

    
    def size(self) -> int: return self.qbits
    def matrix(self) -> SparseRowMatrix: return self.order_matrix()
    def data(self): return [self.N]

    @staticmethod
    def generate_example(_: str, size: int) -> QuantumOrder:
        print(f"%%% [GEN] Generating a valid Order finding instance...")
        return QuantumOrder.random(size)

    @staticmethod
    def generate_header(csv_writer, ttype):
        if ttype in ("clue", "ddsim"):
            csv_writer.writerow(["size", "N", "red. ratio", "time_lumping", "memory (MB)", "problem"])
        elif ttype == "full_clue":
            csv_writer.writerow(["size", "N", "time_lumping", "kappa", "time_iteration", "memory (MB)", "problem"])
        elif ttype == "full_ddsim":
            csv_writer.writerow(["size", "N", "kappa", "time_iteration", "memory (MB)", "problem"])
        else:
            raise NotImplementedError(f"Type of file {ttype} not recognized")

    ## METHODS TO GENERATE OBSERVABLES
    @staticmethod
    def generate_observable_clue(order: QuantumOrder, *_) -> tuple[SparseVector]:
        return tuple([SparseVector.from_list([0,1] + (2**order.size()-2)*[1], field=CC)])

    @staticmethod
    def generate_observable_ddsim(order: QuantumOrder, *_) -> bool:
        raise NotImplementedError

def generate_example(_: str, size: int) -> QuantumOrder:
    return QuantumOrder.generate_example(_, size)

def generate_header(csv_writer, ttype):
    return QuantumOrder.generate_header(csv_writer, ttype)

## METHODS TO GENERATE OBSERVABLES
def generate_observable_clue(order: QuantumOrder, *_) -> tuple[SparseVector]:
    return QuantumOrder.generate_observable_clue(order, *_)

def generate_observable_ddsim(order: QuantumOrder, *_) -> bool:
    return QuantumOrder.generate_observable_ddsim(order, *_)

if __name__ == "__main__":
   ## Processing the arguments
    ## Generic part
    ttype, script = get_method(*sys.argv)
    m, M = get_size_bounds(*sys.argv)
    timeout = get_timeout(*sys.argv)
    repeats = get_repeats(*sys.argv)
    name = "order"; obs = None

    main_script(SCRIPT_DIR, SCRIPT_NAME, name,
                [generate_header, generate_example, generate_observable_clue, generate_observable_ddsim],
                ttype, script,               
                m, M,                        
                timeout, repeats,
                obs)