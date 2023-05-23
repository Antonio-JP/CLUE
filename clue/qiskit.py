r'''
    Module implementing a type of dynamical systems that are based on Quantum Circuits

    This module is somehow a parser that allows to create a dynamical system from 
    the unitary matrix that defines a quantum circuit. This uses the interface by 
    :mod:`qiskit` that allows to manipulate these circuits and generate the 
    corresponding unitary matrices.
'''

from __future__ import annotations

from qiskit import execute, QuantumCircuit, Aer
from numpy import array, matmul, block, cdouble, asarray, sqrt, kron, eye, outer, ones, zeros
from numpy.random import rand
from sympy import CC

from .clue import FODESystem
from .linalg import NumericalSubspace
from .rational_function import SparsePolynomial

class DS_QuantumCircuit(FODESystem):
    BACKEND = Aer.get_backend("unitary_simulator")

    def __init__(self, unitary, **kwds):
        nbits = len(unitary.input_dims())
        variables = [f"Q_{f'{i:b}'.rjust(nbits, '0')}" for i in range(2**nbits)]
        nvariables = 2**nbits
        matr = asarray(unitary)
        equations = [SparsePolynomial.from_vector(matr[i], variables, CC) for i in range(nvariables)]
        
        super().__init__(
            equations, variables=variables, 
            name=kwds.pop("name", "quantum circuit"), 
            lumping_subspace=kwds.pop("lumping_subspace", NumericalSubspace), 
            lumping_subspace_kwds=kwds
        )

    @staticmethod
    def from_qasm_file(path: str, **kwds) -> DS_QuantumCircuit:
        circuit = QuantumCircuit.from_qasm_file(path)
        job = execute(circuit, DS_QuantumCircuit.BACKEND, shots=kwds.pop("shots", 8192))
        unitary = job.result().get_unitary(circuit)

        return DS_QuantumCircuit(unitary, name=circuit.name, **kwds)

###########################################################################
### Building circuit matrices from diagrams
###
### In this part of the file we provide several methods that allow to take 
### numpy matrices that describes part of a quantum circuit and recombines
### them to get a final quantum circuit described with its unitary matrix.
###########################################################################
Hadamard = (1/sqrt(2))*array([[1,1],[1,-1]], dtype=cdouble)
PauliX = array([[0,1],[1,0]], dtype=cdouble); qbit_plus = PauliX
PauliY = array([[0,-1j],[1j,0]], dtype=cdouble)
PauliZ = array([[1,0],[0,-1]], dtype=cdouble)
Phase = array([[1,0],[0,1j]], dtype=cdouble)

numerical_threshold = 1e-10

def ControlledGate(U):
    size = U.shape
    return block([[eye(size[0]), zeros(size)], [zeros(size), U]])

def Entangle(*U):
    return kron(*U)

def inner_product(v, w):
    return matmul(v, w.transpose().conjugate())

def measure(v):
    if len(v.shape) != 1:
        raise TypeError("Only a state can be measured")
    if inner_product(v,v) - 1 > numerical_threshold:
        raise ValueError("The state is not well defined") 
    
    val = rand()
    val -= v[0]*v[0].conjugate(); i = 0
    while val > 0:
        i += 1
        val -= v[i]*v[i].conjugate()

    return i

def G(f, n):
    N = 2**n
    psi = (1/sqrt(N))*ones((N,))

    M = eye(N) - 2*outer(psi,psi)
    f = array([(-1)**f(x) for x in range(N)], dtype=cdouble)
    return M*f

def U(x,N):
    from math import ceil, log2
    if x < 2 or x >= N:
        raise TypeError(f"The defining {x=} must be in the range 2,...,{N}")
    L = ceil(log2(N))
    f = [(x*y)%N if y < N else y for y in range(2**L)]
    U = array([[0 if j != f[i] else 1 for j in range(2**L)] for i in range(2**L)], dtype=cdouble).transpose()

    return U

def K(U):
    return kron(Hadamard, eye(U.shape[0])) * ControlledGate(U) * kron(Hadamard, eye(U.shape[0]))
    

