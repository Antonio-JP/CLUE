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

from math import sqrt, log10, floor
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
    return Circuit.from_list(kron(A.to_numpy(CC), B.to_numpy(CC)), CC)

def kron_pow(A: Circuit, n : int) -> Circuit:
    result = A
    for _ in range(1,n):
        result = kronecker(result, A)

    return result

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

# Identity 8
I8 = Circuit.eye(8, CC)

def CnNOT(n : int) -> Circuit:
    N = 2**n
    output = Circuit.eye(N,CC)
    output.increment(N-1, N-1, -1)
    output.increment(N-1, N-2, 1)
    
    output.increment(N-2, N-1, 1)
    output.increment(N-2, N-2, -1)

    return output

def not_CnNOT(n: int) -> Circuit:
    N = 2**n
    output = Circuit.eye(N,CC)
    output.increment(0, 0, -1)
    output.increment(0, 1, 1)
    
    output.increment(1, 0, 1)
    output.increment(1, 1, -1)

    return output

def G(n: int, epsilon: float) -> DensityOperator:
    O = CnNOT(n+1)
    P = [kronecker(kron_pow(H, n), I), kronecker(kron_pow(I, n), X), not_CnNOT(n+1), kronecker(kron_pow(H, n), I)]

    In = Circuit.eye(2**(n+1), CC)

    operators = [DensityOperator(circuits=[circ, In], probabilities=[1-epsilon, epsilon]) for circ in [O] + P]
    return DensityOperator(operators=operators)

def G_input(n: int) -> DensityVector:
    v = State(2**(n+1), CC)
    v[1] = 1

    return DensityVector.from_tensor(v.apply_matrix(kron_pow(H, n+1)))

def run(G: DensityOperator, v: DensityVector):
    return find_smallest_common_subspace(
        (G,),
        (v,),
        subspace_class=NumericalSubspace
    )

def evolution(U, v, starting:float, finishing:float, increase="log") -> tuple[tuple[float,int]]:
    epsilon = starting
    result = []
    while epsilon < finishing:
        print(f"Computing the reduction with noise={epsilon:.04f}", flush=True, end="\r")
        S = run(U(epsilon), v(epsilon))
        result.append((epsilon,S.dim()))

        if increase == "log":
            epsilon = epsilon + 10**floor(log10(epsilon))
        elif increase == "linear":
            epsilon += starting

    return tuple(result)

import matplotlib.pyplot as plt

def plot(result: tuple[tuple[float, int]], scale: str = "linear"):
    xvalues, yvalues = list(zip(*result))
    plt.plot(xvalues, yvalues, 'o', linestyle="-")
    plt.xscale(scale)
    plt.show()
