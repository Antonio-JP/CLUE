r'''
    Small set of functions to test the numerical lumping in Quantum experiments
'''
from __future__ import annotations
import sys, os

SCRIPT_DIR = os.path.dirname(__file__) if __name__ != "__main__" else "./"
sys.path.insert(0, os.path.join(SCRIPT_DIR, "..", "..")) # clue is here

import matplotlib.pyplot as plt 

from functools import lru_cache
from shutil import get_terminal_size

from clue import FODESystem, LDESystem
from clue.linalg import CC, OrthogonalSubspace, NumericalSubspace, SparseVector, SparseRowMatrix
from math import ceil
from numpy.linalg import svd
from numpy import matmul
from misc import Experiment
from q_sat import SATFormula
from q_maxcut import UndirectedGraph

import pickle
from pickle import PicklingError
from numpy import exp, pi

Analysis = dict[int, tuple[float,LDESystem]]

###################################################################################
### 
### METHODS FOR CREATING THE SYSTEMS AND THE OBSERVABLES FOR EACH EXPERIMENT
###
###################################################################################
def observable(E: Experiment) -> SparseVector:
    v = E.generate_observable_clue(E, E.size(), 0)[0]
    if v.norm() != 1.0:
        v.scale(1.0/v.norm())
    return v

__CACHED_Us = []
def quantum_matrix(E: Experiment) -> SparseRowMatrix:
    for (e, U) in __CACHED_Us:
        if e == E:
            return U
        
    ## Not cached version
    M = E.matrix()
    if isinstance(E, (SATFormula, UndirectedGraph)):
        par = 1/(3*max(M[i][i].real for i in range(M.nrows)))
        U = SparseRowMatrix(M.nrows, CC)
        for i in range(U.nrows):
            in_part = -1j*float(-par*pi*M[i][i].real)
            U.increment(i,i, CC(exp(in_part)))
    else:
        U = M

    __CACHED_Us.append((E, U))
    return U

def app_lumping(E: Experiment, delta: float=1e-4) -> LDESystem:
    U = quantum_matrix(E)
    if delta == 0: # exact 
        system = FODESystem.LinearSystem(U, lumping_subspace=OrthogonalSubspace)
        obs = (observable(E),)
        return system.lumping(obs, print_reduction=False, print_system=False, out_format="internal")
    else: # numerical
        system = FODESystem.LinearSystem(U, lumping_subspace=NumericalSubspace, lumping_subspace_kwds={"delta": delta})
        obs = (observable(E),)
        return system.lumping(obs, print_reduction=False, print_system=False, out_format="internal")
    
def max_epsilon(G, threshold=1e-10):
    U = quantum_matrix(G)
    system = FODESystem.LinearSystem(U, lumping_subspace=NumericalSubspace, lumping_subspace_kwds={"delta": threshold**2})
    return system.find_maximal_threshold((observable(G),))

###################################################################################
###
### METHODS FOR ANALYSIS OF A QUANTUM SYSTEM
###
###################################################################################
__CACHED_ANALYSIS = []
def analysis(E: Experiment, threshold=1e-10, min_epsilon=0) -> Analysis:
    tot_lumpings = 0

    ## Looking into CACHE
    for (e, th, A) in __CACHED_ANALYSIS:
        if (e == E) and (th >= threshold):
            print(f"[analysis] COMPLETED (with {tot_lumpings} lumpings) (from CACHE)".ljust(get_terminal_size()[0]))
            return A

    ## Looking in the files
    path_name = f"./analysis/{f'{E.name}_{E.size()}' if hasattr(E, 'name') else f'(({E}))'}[{threshold:.2E}]({min_epsilon:.2f}).an"
    if os.path.exists(path_name):
        with open(path_name, "rb") as file:
            print(f"[analysis] Loading analysis from file...".ljust(get_terminal_size()[0]), end="\r")
            A = pickle.load(file)

            print(f"[analysis] COMPLETED (with {tot_lumpings} lumpings) (Read from file)".ljust(get_terminal_size()[0]))
            return A
    
    ## Computing the ANALYSIS
    print(f"[analysis] Computing maximum epsilon".ljust(get_terminal_size()[0]), end="\r")
    me, Me = max(min_epsilon, 0) + threshold/2., max_epsilon(E, threshold)
    print(f"[analysis] Computing lumpings with given min-max epsilons...".ljust(get_terminal_size()[0]), end="\r")
    result = {me: app_lumping(E, me), Me: app_lumping(E, Me)}
    tot_lumpings = 2
    cache = True

    intervals = [(me,Me)] if (Me-me > threshold) else list()
    while len(intervals) > 0:
        try:
            print(f"[analysis] Remaining intervals: {len(intervals)}/{len(result)}".ljust(get_terminal_size()[0]), end="\r")
            me, Me = intervals.pop()
            ml = result[me]
            Ml = result[Me]
            
            ## Computing new lumping in the middle
            ne = (Me+me)/2.
            nl = app_lumping(E, ne)
            tot_lumpings += 1
            result[ne] = nl
            
            ## Adding new intervals if needed
            if ml.size != nl.size and (ne-me) > threshold:
                intervals.append((me, ne))
            if Ml.size != nl.size and (Me-ne) > threshold:
                intervals.append((ne, Me))
        except KeyboardInterrupt:
            print(f"[analysis] Broken analysis. We keep working but we do not cache the result")
            cache = False
    
    # We create a list with the (epsilon, size)
    print(f"[analysis] Sorting results...".ljust(get_terminal_size()[0]), end="\r")
    result = sorted((el, result[el].size) for el in result)
    sizes = set(el[1] for el in result)
    result = {s: (min(el[0] for el in result if el[1] == s), max(el[0] for el in result if el[1] == s)) for s in sizes}

    # Computing middle lumpings
    print(f"[analysis] Computing middle lumpings".ljust(get_terminal_size()[0]), end="\r")
    A = dict()
    for i, (s, v) in enumerate(result.items()):
        print(f"[analysis] Computing middle lumpings ({i+1}/{len(sizes)})...".ljust(get_terminal_size()[0]), end="\r")
        eps = (v[0]+v[1])/2.
        lum = app_lumping(E, eps)
        tot_lumpings += 1
        A[s] = (eps, lum)
        
    print(f"[analysis] COMPLETED: Saving to CACHE".ljust(get_terminal_size()[0]), end="\r")
    if cache:
        __CACHED_ANALYSIS.append((E, threshold, A))

    print(f"[analysis] COMPLETED: Saving to FILE".ljust(get_terminal_size()[0]), end="\r")
    with open(path_name, "wb") as file:
        try:
            print(f"[analysis] Saving analysis into file for future use...".ljust(get_terminal_size()[0]), end="\r")
            pickle.dump(A, file)
        except PicklingError as e:
            print(f"[analysis] COMPLETED (with {tot_lumpings} lumpings) (Error saving in file)".ljust(get_terminal_size()[0]))
            os.remove(path_name)
            return A
    print(f"[analysis] COMPLETED: Saved in file".ljust(get_terminal_size()[0]))
    return A

def epsilon_for_analysis(A: Analysis, size: int) -> float:
    if not size in A:
        raise ValueError("That size is not possible with Numerical Lumping")
    return A[size][0]

def lumping_from_analysis(A: Analysis, size: int) -> LDESystem:
    if not size in A:
        raise ValueError("That size is not possible with Numerical Lumping")
    return A[size][1]

def save_analysis(E: Experiment, A: Analysis):
    __CACHED_ANALYSIS.append((E, float("inf"), A))

def epsilon_intervals(E: Experiment | Analysis, threshold=1e-10) -> tuple[float,float]:
    if isinstance(E, Experiment):
        E = analysis(E, threshold)
    
    return dict(zip(sorted(E.keys(), reverse=True), _epsilon_interval(E)))

__CACHED_INTERVALS = dict()
def _epsilon_interval(A: Analysis) -> list[tuple[float, float]]:
    key = tuple(sorted(A.items(), reverse=True))
    if key not in __CACHED_INTERVALS:
        pairs = sorted([(k, e) for (k, (e, _)) in A.items()],reverse=True)
        c_epsilon = 0.0
        intervals = []
        for pair in pairs:
            intervals.append((c_epsilon, 2*pair[1] - c_epsilon))
            c_epsilon = 2*pair[1] - c_epsilon

        __CACHED_INTERVALS[key] = intervals
    return __CACHED_INTERVALS[key]
        
###################################################################################
###
### METHODS FOR COMPUTING THE ERROR WITH THE LUMPED SYSTEM
###
###################################################################################
def sample_ks(size=10, num_samples:int = 10):
    up_bound = 2**(int(ceil(size/2)))

    return [int(ceil(j*(up_bound-1)/(num_samples-1)))+1 for j in range(num_samples)]

@lru_cache(maxsize=512)
def matrix_power(A: SparseRowMatrix, power:int) -> SparseRowMatrix:
    r'''Method to compute "quickly" the power of a matrix'''
    if A.nrows != A.ncols:
        raise ValueError("Matrix must be square")
    
    if power < 0:
        raise ValueError(f"Cannot raise to power {power}, {str(A)}")
    if power == 0:
        return A.eye(A.nrows, A.field)
    if power == 1:
        return A
    return matrix_power(A, power // 2) * matrix_power(A, power // 2 + (power % 2))

def matrices_example(E: Experiment, size: int, threshold=1e-10) -> tuple[SparseRowMatrix, SparseRowMatrix, SparseRowMatrix, SparseRowMatrix]:
    r'''For a given size, this returns L, L_+, U, U_hat'''
    A = analysis(E, threshold) # This is only computed once
    al = lumping_from_analysis(A, size)

    U = quantum_matrix(E)
    L = al.lumping_matrix
    L_plus = SparseRowMatrix.from_list([[L[i][j].conjugate() for i in range(L.nrows)] for j in range(L.ncols)], field=L.field)
    U_hat = L.matmul(U).matmul(L_plus)
    
    return L, L_plus, U, U_hat

def closest_unitary(Uhat) -> SparseRowMatrix:
    W, _, V = svd(Uhat.to_numpy(dtype="complex"))
    return SparseRowMatrix.from_list(matmul(W, V), CC)

def direct_error(E: Experiment, size: int, iter_size: int = None, *, threshold) -> dict[int, SparseVector]:
    r'''
        When we compute a lumping `L` for a circuit `U`, we can build a reduced model by

        .. MATH::

            \hat U = L U L^+,

        where `L` projects `\mathbb{C}^N` to `\mathbb{C}^m` orthogonally and `L^+` embedded the subspace generated by `L` into the ambient space `\mathbb{C}^N`.

        In particular, `L^m \in \mathbb{C}^{n \times m}` is the right pseudo-inverse of `L`, providing `LL^+ = I_m`. Moreover, for any `|x\rangle \in \langle L\rangle`, 
        we have that `L^+L |x\rangle = |x\rangle`, so we can show that 
        
        .. MATH::
        
            U^k |x\rangle = L^+ \hat U^k L |x\rangle.

        If `L` is not a lumping, then this property does not hold. So we can measure how far we get from the true value of the lumped simulation. Namely:

        .. MATH::
        
            ||(U^k - L^+ \hat U^k L)|x\rangle||_2
    '''
    print(f"[direct @ {size}] Computing the matrices to compare the error...".ljust(get_terminal_size()[0]), end="\r")
    L, L_plus, U, Uhat = matrices_example(E, size, threshold)
    x = SparseVector.from_list(observable(E).to_list(), L.field)
    
    k_values = sample_ks(E.size() if iter_size is None else iter_size)
    differences = []
    for k in k_values:
        print(f"[direct @ {size}] Computing error after {k}/{k_values[-1]} iterations...".ljust(get_terminal_size()[0]), end="\r")
    
        Up = matrix_power(U, k)
        Uhp = matrix_power(Uhat, k)
        
        differences.append((Up - L_plus.matmul(Uhp.matmul(L))).dot(x))
        
    print(f"[direct @ {size}] COMPLETED".ljust(get_terminal_size()[0]))
    return dict(zip(k_values, differences))

def closest_unitary_error(E: Experiment, size: int, iter_size: int = None, *, threshold) -> dict[int, SparseVector]:
    r'''
        In :func:`direct_error`, we simply rolled with the approximate lumping we obtained. However, the reduced system 
        `\hat U` was not unitary, i.e., it was not a quantum circuit anymore. Since this matrix should be "close" to a 
        unitary matrix (the true lumped system) we can try to compute the closest unitary matrix.

        This can be achieve using the polar decomposition of a matrix (see [here](https://en.wikipedia.org/wiki/Polar_decomposition))
    '''
    print(f"[unitary @ {size}] Computing the matrices to compare the error...".ljust(get_terminal_size()[0]), end="\r")
    L, L_plus, U, Uhat = matrices_example(E, size, threshold)
    print("[unitary @ {size}] Computing closes unitary...".ljust(get_terminal_size()[0]), end="\r")
    nU = closest_unitary(Uhat)
    x = SparseVector.from_list(observable(E).to_list(), L.field)
    
    k_values = sample_ks(E.size() if iter_size is None else iter_size)
    differences = []
    for k in k_values:
        print(f"[unitary @ {size}] Computing error after {k}/{k_values[-1]} iterations...".ljust(get_terminal_size()[0]), end="\r")
        
        Up = matrix_power(U, k)
        nUp = matrix_power(nU, k)

        differences.append((Up - L_plus.matmul(nUp.matmul(L))).dot(x))
        
    print(f"[unitary @ {size}] COMPLETED".ljust(get_terminal_size()[0]))
    return dict(zip(k_values, differences))

def generate_error_graph(E: Experiment, method=direct_error, name="\hat{U}", 
                         yscale=None, xscale = None, 
                         bound_lump : int |tuple[int,int] = None, iter_size: int = None, *, threshold):
    print(f"[graph] Generating the Error graph for {method.__name__}".ljust(get_terminal_size()[0]), end="\r")
    A = analysis(E, threshold) # this is done just once
    if isinstance(bound_lump, (tuple,list)):
        m_bound, M_bound = bound_lump
    elif isinstance(bound_lump, int):
        m_bound, M_bound = 0, bound_lump
    elif bound_lump is None:
        m_bound, M_bound = 0, 2**E.size()

    d = [(s,method(E, s, iter_size, threshold=threshold)) for s in A if s < M_bound and s > m_bound]
    for (s,e) in d:
        x = e.keys()
        y = [v.norm() for v in e.values()]
        plt.plot(x, y, label=f"s={s}", linestyle="-") 
        
    plt.legend()
    plt.xlabel('Number of iterations (k)')
    plt.ylabel(f'$||U^k - L^+ {name}^k L||_2$')
    if yscale is not None: plt.yscale(yscale)
    if xscale is not None: plt.xscale(xscale)
    
    plt.show()
    print(f"[graph] COMPLETED".ljust(get_terminal_size()[0]))

__all__ = [
    "observable", "quantum_matrix", "app_lumping", "max_epsilon",
    "analysis", "save_analysis", "epsilon_intervals",
    "matrices_example", "closest_unitary", "direct_error", "closest_unitary_error", "generate_error_graph"
]