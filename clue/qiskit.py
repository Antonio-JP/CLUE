r'''
    Module implementing a type of dynamical systems that are based on Quantum Circuits

    This module is somehow a parser that allows to create a dynamical system from 
    the unitary matrix that defines a quantum circuit. This uses the interface by 
    :mod:`qiskit` that allows to manipulate these circuits and generate the 
    corresponding unitary matrices.
'''

from __future__ import annotations

from itertools import product
from functools import lru_cache
from math import ceil, log2
from numpy import array, matmul, block, cdouble, asarray, sqrt, kron, eye, outer, ones, zeros, cos, sin, pi 
from numpy.random import rand
from qiskit import execute, QuantumCircuit, Aer
from sympy import CC

from .clue import FODESystem
from .linalg import NumericalSubspace
from .rational_function import SparsePolynomial

class DS_QuantumCircuit(FODESystem):
    BACKEND = Aer.get_backend("unitary_simulator")

    def __init__(self, unitary, **kwds):
        nbits = int(log2(unitary.shape[0]))
        if 2**nbits != unitary.shape[0]:
            raise ValueError("The unitary matrix is not of size 2^N")
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
# Gates for 1 q-bit
Hadamard = (1/sqrt(2))*array([[1,1],[1,-1]], dtype=cdouble)
PauliX = array([[0,1],[1,0]], dtype=cdouble); qbit_plus = PauliX
PauliY = array([[0,-1j],[1j,0]], dtype=cdouble)
PauliZ = array([[1,0],[0,-1]], dtype=cdouble)
Phase = array([[1,0],[0,1j]], dtype=cdouble)

@lru_cache(maxsize=10)
def R(k: int):
    r'''Rotation map'''
    return array([[1, 0], [0, cos(2*pi / 2**k) + sin(2*pi / 2**k)*1j]], dtype=cdouble)

@lru_cache(maxsize=10)
def S(n: int):
    r'''Swap q-bits'''
    result = eye(2**n, dtype=cdouble)
    for i,t in enumerate(product([0], *((n-1)*[[0,1]]))):
        v = sum(2**j * e for (j,e) in enumerate(t))
        if v != i: # non-symmetric number
            result[[i,v]] = result[[v,i]] # swapping rows i <-> v

    return result.transpose()

BitPlus = PauliX

# Numerical threshold for checking equality
numerical_threshold = 1e-10

# Methods to manipulate state and to measure a quantum state
def inner_product(v, w):
    r'''Method for inner or scalar product of two elements with complex elements'''
    return matmul(v, w.transpose().conjugate())

def measure(v):
    r'''
        Method to measure a quantum state.

        Given a vector representing a quantum state in terms of a basis of states, this method measures 
        the state within `v` where the probabilities of ending up in a basis state is given by the squared
        modulus of the coefficients of `v`.
    '''
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

def is_unitary(U):
    r'''Method that checks if a matrix is unitary or not'''
    return ((matmul(U, U.transpose().conjugate()) - eye(U.shape[0])) < numerical_threshold).all()

# Methods to create quantum gates
@lru_cache(maxsize=10)
def EntangledState(n: int, dtype=None):
    r'''Created the entangled state for `n` different states.'''
    return (1/sqrt(n))*ones((n,), dtype=dtype)

Psi = EntangledState #: alias for EntangledState

def ApplyWithControl(U, src: int, dst: int, tot: int):
    r'''Applies to the ``dst`` q-bit the gate ``U`` controlled by the q-bit ``src``'''
    if any(bit <= 0 or bit > tot for bit in (src, dst)):
        raise ValueError(f"The control and target q-bit must be between 1 and {tot}")
    
    size = 2**tot
    val = lambda b : 2**(tot-b)
    is_active = lambda p,b: (p//val(b)%2) == 1
    relative_vals = lambda p, b: (p-val(b), p) if is_active(p, b) else (p, p+val(b))
    result = eye(size, dtype=cdouble)
    for i in range(size):
        if is_active(i, src):
            dst_0, dst_1 = relative_vals(i, dst)
            is_on = 1 if is_active(i, dst) else 0
            result[dst_0][i] = U[0][is_on]; result[dst_1][i] = U[1][is_on]

    return result
            
def ControlledGate(U, n: int=1):
    r'''
        Building the controlled gate with `n` q-bits.

        A controlled gate applies a given gate `U` if a set of `q`-bits are all 1. If not, the controlled q-bits
        stays the same. The q-bits used for controlling do not change through the Controlled Gate.

        More info:

        * :wiki:`Controlled_NOT_gate`
        * :wiki:`Quantum_logic_gate`

        Input:

        * ``U``: unitary matrix that define a quantum gate.
        * ``n``: number of controlling q-bits. They will be add to the beginning of the string of q-bits.

        Output:

        The unitary matrix that defines the corresponding controlled gate of `U` using `n` controlling q-bits.

        Example::

            >>> from clue.qiskit import *
            >>> ControlledGate(PauliX)
            array([[1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j],
                   [0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j]])
            >>> ControlledGate(PauliX, 2)
            array([[1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j]])
    '''
    size = U.shape
    identity_block = eye(2**n * (size[0]-1))
    return block([[identity_block, zeros((identity_block.shape[0], size[1]))], [zeros((size[0], identity_block.shape[1])), U]])

def G(f, n):
    r'''
        Method to create the Grover gate.

        This method receives the search function (i.e., a boolean function such that `f(x) = 1` indicates 
        we succeed in our search and `f(x) = 0` means we failed in the search) and the number of q-bits
        for the gate.

        This method creates the unitary matrix that defines the Grover's gate used in the search algorithm.

        Input:

        * ``f``: the boolean function defining the oracle.
        * ``n``: number of q-bits used in the Grover's gate.

        Output:

        Unitary matrix defining the Grover's gate.

        Examples::

            >>> from clue.qiskit import *
            >>> f = lambda p : 0 if p != 1 else 1 # looking number 1
            >>> G(f, 2)
            array([[ 0.5+0.j,  0.5-0.j, -0.5+0.j, -0.5+0.j],
                   [-0.5+0.j, -0.5+0.j, -0.5+0.j, -0.5+0.j],
                   [-0.5+0.j,  0.5-0.j,  0.5+0.j, -0.5+0.j],
                   [-0.5+0.j,  0.5-0.j, -0.5+0.j,  0.5+0.j]])
            >>> f = lambda p : 1 if p != 1 else 0 # looking number 1
            >>> G(f, 2)
            array([[-0.5+0.j, -0.5+0.j,  0.5-0.j,  0.5-0.j],
                   [ 0.5-0.j,  0.5+0.j,  0.5-0.j,  0.5-0.j],
                   [ 0.5-0.j, -0.5+0.j, -0.5+0.j,  0.5-0.j],
                   [ 0.5-0.j, -0.5+0.j,  0.5-0.j, -0.5+0.j]])
            
    '''
    N = 2**n
    psi = Psi(N)

    M = eye(N) - 2*outer(psi,psi)
    f = array([(-1)**f(x) for x in range(N)], dtype=cdouble)
    return M*f

@lru_cache(maxsize=256)
def U(x,N):
    r'''
        Method to create the operator "multiplication by `x`".

        This method creates the unitary matrix that represents the operation

        .. MATH::

            f(y) = xy (mod N)

        for given values of `x` and `N`. This operation is performed using `L` q-bits where 
        `L` is the minimal value such that `N < 2^L`. The values between `N` and `2^L` are
        not changed by this gate.

        Input:

        * ``x``: number we are going to multiply in this gate.
        * ``N``: Value we are taking the final modulo.

        Output:

        Unitary matrix that represent the multiplication by `x` module `N`.

        Examples::

            >>> from clue.qiskit import *
            >>> U(2, 4)
            array([[1.+0.j, 0.+0.j, 1.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 1.+0.j, 0.+0.j, 1.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j]])
            >>> U(3, 4)
            array([[1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j],
                   [0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j],
                   [0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j]])
    '''
    if x < 2 or x >= N:
        raise TypeError(f"The defining {x=} must be in the range 2,...,{N}")
    L = ceil(log2(N))
    f = [(x*y)%N if y < N else y for y in range(2**L)]
    U = array([[0 if j != f[i] else 1 for j in range(2**L)] for i in range(2**L)], dtype=cdouble).transpose()

    return U

@lru_cache(maxsize=10)
def base_Fourier(n: int):
    r'''
        Method to create the quantum gate for QFT.
        
        It follows the representation of the circuit in page 219 of the Nielsen-Chuang book.
    '''
    if n == 1:
        return Hadamard
    else:
        aux = base_Fourier(n-1)
        last = kron(eye(2, dtype=cdouble), aux)

        first = kron(Hadamard, eye(aux.shape[0], dtype=cdouble))
        for i in range(2, n+1):
            first = matmul(ApplyWithControl(R(i), i, 1, n), first)

        return matmul(last,first)
    
@lru_cache(maxsize=10)
def FourierTransform(n: int):
    r'''
        Method to create the quantum gate for QFT.
        
        It follows the representation of the circuit in page 219 of the Nielsen-Chuang book.
    '''
    return matmul(S(n), base_Fourier(n))
    
@lru_cache(maxsize=10)
def FourierTransform_direct(n: int):
    root_of_unity = cos(2*pi /(2**n)) + sin(2*pi / (2**n))*1j
    main_vector = array([root_of_unity**i for i in range(2**n)], dtype=cdouble)

    rows = [ones(2**n, dtype=cdouble)]
    for i in range(2**n - 1):
        rows += [rows[-1]*main_vector]
    
    return (1/sqrt(2**n)) * array(rows, dtype=cdouble)

def K(U):
    r'''Method to create the Kitaev gate for phase estimation'''
    return matmul(kron(Hadamard, eye(U.shape[0])), matmul(ControlledGate(U), kron(Hadamard, eye(U.shape[0]))))

