import sys
import platform

# clue is here
if platform.system() == 'Linux':
    sys.path.insert(0, "../..")
elif platform.system() == "Windows":
    sys.path.insert(0, "..\..")


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

zero = State(8,CC)
zero[0] = 1

def kronecker(A: Circuit, B: Circuit):
    return Circuit.from_list(kron(A.to_numpy(CC), B.to_numpy(CC)))

# Hadamard Gate
# 1/sqrt2 * [1, 1]
#           [1,-1]
H = Circuit(2,CC)
H.increment(0,0,1/sqrt(2));H.increment(0,1,1/sqrt(2))
H.increment(1,0,1/sqrt(2));H.increment(1,1,-1/sqrt(2))

# CNOT Gate
# [1,0,0,0]
# [0,1,0,0]
# [0,0,0,1]
# [0,0,1,0]
CX = Circuit(4,CC)
CX.increment(0,0,1)
CX.increment(1,1,1)
CX.increment(2,3,1)
CX.increment(3,2,1)

# Matrix for the composition of each layer in GHZ
U_1 = kronecker(kronecker(H,I),I) # Goes from 4x4 (the first kronecker) to 8x8 (the second)
U_2 = kronecker(CX,I) # 4x4 -> 8x8
U_3 = kronecker(I,CX) # 2x2 -> 8x8

epsilon = 0.5
# DensityOperator for each of the layers
op1 = DensityOperator(circuits=[U_1,I], probabilities=[1-epsilon,epsilon])
op2 = DensityOperator(circuits=[U_2,I], probabilities=[1-epsilon,epsilon])
op3 = DensityOperator(circuits=[U_3,I], probabilities=[1-epsilon,epsilon])

def run():
    ## We want to try to use the method `find_smallest_common_subspace` using these gates and states
    ## The matrices will be the Density Operator with some probability epsilon of doing nothing
    return find_smallest_common_subspace(
        (DensityOperator(operators=[op1,op2,op3]),),
        (DensityVector.from_tensor(zero),),
        subspace_class=NumericalSubspace
    )