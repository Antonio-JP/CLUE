import sys

sys.path.insert(0, "../..") # clue is here

from clue.linalg import SparseRowMatrix as Circuit, SparseVector as State, NumericalSubspace, find_smallest_common_subspace
from clue.numerical_domains import CC
from clue.quantum_linalg import DensityOperator, DensityVector

from math import sqrt
from numpy import kron

X = Circuit(2, CC)
X.increment(1,0,1)
X.increment(0,1,1)
Y = Circuit(2, CC)
Y.increment(1,0,CC(1j))
Y.increment(0,1,CC(-1j))

I = Circuit.eye(2, CC)

plus = State(2, CC)
plus[0], plus[1] = 1/sqrt(2), 1/sqrt(2)

minus = State(2, CC)
minus[0], minus[1] = 1/sqrt(2), -1/sqrt(2)

def kronecker(A: Circuit, B: Circuit):
    return Circuit.from_list(kron(A.to_numpy(CC), B.to_numpy(CC)))

epsilon = 0.5

def run():
    ## We want to try to use the method `find_smallest_common_subspace` using these gates and states
    ## The matrices will be the Density Operator with some probability epsilon of doing nothing
    return find_smallest_common_subspace(
        (DensityOperator(circuits=[Y, I], probabilities=[1-epsilon, epsilon]),),
        (DensityVector.from_tensor(minus),),
        subspace_class=NumericalSubspace
    )
