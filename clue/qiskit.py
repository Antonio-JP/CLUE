r'''
    Module implementing a type of dynamical systems that are based on Quantum Circuits

    This module is somehow a parser that allows to create a dynamical system from 
    the unitary matrix that defines a quantum circuit. This uses the interface by 
    :mod:`qiskit` that allows to manipulate these circuits and generate the 
    corresponding unitary matrices.
'''

from __future__ import annotations

from collections.abc import Sequence
from itertools import chain, product
from functools import lru_cache, reduce
from math import ceil, log2
from numpy import array, matmul, block, cdouble, asarray, sqrt, kron, eye, outer, ones, zeros, exp, pi, arange
from numpy.linalg import inv, matrix_power
from numpy.random import rand
from qiskit import execute, QuantumCircuit, Aer
from sympy import CC

from .clue import FODESystem
from .linalg import NumericalSubspace
from .rational_function import SparsePolynomial

class DS_QuantumCircuit(FODESystem):
    r'''
        Class for representing quantum circuits as dynamical systems.
    '''
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
        circuit.remove_final_measurements(True) # we remove the final measeurements (hence, we do not need to remove them from the .qasm file)
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
# Numerical threshold for checking equality
numerical_threshold = 1e-10

def extend_to_power(U):
    r'''This method extends a unitary matrix to have a size a power of 2, so it fits into a quantum circuit'''
    to_add = 2**ceil(log2(U.shape[0])) - U.shape[0]
    return block([[U, zeros((U.shape[0], to_add))], [zeros((to_add, U.shape[1])), eye(to_add)]])

def compare(a, b, precision=10):
    r'''Method to compare up to a precision. If the difference or the "all" return an error, we return "False"'''
    try:
        return ((a-b).round(precision) == 0).all()
    except:
        return False

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
    
    return repeated_measure(v, 1)[0]

def repeated_measure(v, shots: int=1, out=list):
    r'''This method assumes `v` is a vector with module 1'''
    from random import choices
    output = choices(range(len(v)), v*v.conjugate(), k = shots)
    if out == dict:
        from collections import Counter
        output = Counter(output)
    return output

def is_unitary(U):
    r'''Method that checks if a matrix is unitary or not'''
    return ((matmul(U, U.transpose().conjugate()) - eye(U.shape[0])) < numerical_threshold).all()

def is_symplectic(U):
    r,c = U.shape
    if r != c: raise AssertionError("The matrix is not square")
    if r == 1: return True
    if r%2 != 0: raise AssertionError("The matrix is not square")
    om_r = r//2
    OM = block([
        [zeros((om_r,om_r),dtype=U.dtype), eye(om_r, om_r,dtype=U.dtype)], 
        [-eye(om_r, om_r,dtype=U.dtype), zeros((om_r,om_r),dtype=U.dtype)]
    ])
    return (matmul(matmul(U.transpose().conjugate(), OM), U) - OM < numerical_threshold).all()

# Gates for 1 q-bit
Hadamard = (1/sqrt(2))*array([[1,1],[1,-1]], dtype=cdouble)
PauliX = array([[0,1],[1,0]], dtype=cdouble); qbit_plus = PauliX
PauliY = array([[0,-1j],[1j,0]], dtype=cdouble)
PauliZ = array([[1,0],[0,-1]], dtype=cdouble)
Phase = array([[1,0],[0,1j]], dtype=cdouble)

@lru_cache(maxsize=10)
def Rot(k: int):
    r'''
        Rotation map
        
        Examples::

            >>> compare(Rot(1), PauliZ)
            True
            >>> Rot(2).round(15)
            array([[1.+0.j, 0.+0.j],
                   [0.+0.j, 0.+1.j]])


    '''
    return array([[1, 0], [0, exp(2j*pi / (2**k))]], dtype=cdouble)

@lru_cache(maxsize=10)
def Sw(n: int):
    r'''
        Swap q-bits
    
        Examples::

            >>> compare(Sw(1), eye(2))
            True
            >>> Sw(2)
            array([[1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j],
                   [0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j]])
    '''
    result = eye(2**n, dtype=cdouble)
    for i,t in enumerate(product([0], *((n-1)*[[0,1]]))):
        v = sum(2**j * e for (j,e) in enumerate(t))
        if v != i: # non-symmetric number
            result[[i,v]] = result[[v,i]] # swapping rows i <-> v

    return result.transpose()

BitPlus = PauliX #: alias for the PauliX 

# Methods to create quantum gates
@lru_cache(maxsize=10)
def EntangledState(n: int, dtype=cdouble):
    r'''Created the entangled state for `n` different states.'''
    return (1/sqrt(n))*ones((n,), dtype=dtype)

Psi = EntangledState #: alias for EntangledState

def PermuteBits(permutation: Sequence[int]):
    r'''
        Gate that permutes the q-bits in order.

        If the given permutation is `(\sigma_0,\ldots,\sigma_{n-1})`, then the quantum circuit will 
        move to `i`-th position the `\sigma_i` q-bit of the input.

        Examples::

            >>> PermuteBits([2,1,0]) # inverting the order
            array([[1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j],
                   [0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j]])

        This is exactly the same as the ``Sw`` gate::

            >>> compare(PermuteBits([2,1,0]), Sw(3))
            True

        However, this gate allows to move q-bits easier::

            >>> PermuteBits([2,0,1]) # moving (b_0,b_1,b_2) -> (b_2,b_0,b_1)
            array([[1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j],
                   [0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j]])
            >>> from numpy.linalg import matrix_power
            >>> compare(matrix_power(PermuteBits([2,0,1]),3), eye(2**3))
            True


    '''
    if any(i not in permutation for i in range(len(permutation))):
        raise TypeError("The given list is not a permutation")
    
    n = len(permutation); size = 2**n
    result = zeros((size,size), dtype=cdouble)
    def canonical(size, ind, dtype=cdouble): 
        out = zeros(size, dtype=dtype)
        out[ind] = 1
        return out

    for (k, bits) in enumerate(product([0,1], repeat=n)):
        permuted_value = sum(2**(n-1-i) * bits[permutation[i]] for i in range(n))
        result[:,k] = canonical(size,permuted_value)

    return result

def ApplyWithControl(U, src: int|Sequence[int], dst: int|Sequence[int], tot: int):
    r'''
        Applies to the ``dst`` q-bit the gate ``U`` controlled by the q-bit ``src``
    '''
    src = [src] if isinstance(src, int) else src
    dst = [dst] if isinstance(dst, int) else dst

    if any(bit < 0 or bit >= tot for bit in chain(src, dst)):
        raise ValueError(f"The control and target q-bit must be between 0 and {tot-1}")
    if U.shape != (2**len(dst),2**len(dst)):
        raise TypeError("The number of goal q-bits must match the size of the given gate")
    if any(bit in dst for bit in src):
        raise ValueError("The control and the target can not intersect")
    
    involved = src + dst
    permutation = [i for i in range(tot) if i not in involved] + src + dst
    inverse_permutation = [permutation.index(i) for i in range(tot)]

    first = PermuteBits(permutation)
    inner = kron(eye(2**(tot-len(src)-len(dst))), ControlledGate(U, len(src)))
    last = PermuteBits(inverse_permutation)

    return matmul(last, matmul(inner, first))
            
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
            array([[1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j],
                   [0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j]])
    '''
    size = U.shape
    identity_block = eye((2**n - 1)* size[0])
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
            >>> f = lambda p : 1 if p != 1 else 0 # looking anything but 1
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
        size = 2**n

        # We create the first gates
        first = kron(Hadamard, eye(size//2, dtype=cdouble))
        for i in range(2, n+1):
            first = matmul(ApplyWithControl(Rot(i), i-1, 0, n), first)

        # We then apply the base Fourier of one less q-bits
        last = kron(eye(2, dtype=cdouble), base_Fourier(n-1)) # adding one extra bit to previous unswept Fourier

        return matmul(last,first)
    
@lru_cache(maxsize=10)
def FourierTransform(n: int):
    r'''
        Method to create the quantum gate for QFT.
        
        It follows the representation of the circuit in page 219 of the Nielsen-Chuang book.
    '''
    return matmul(Sw(n), base_Fourier(n))
    
@lru_cache(maxsize=10)
def FourierTransform_direct(n: int):
    from numpy import vander
    size = 2**n
    roots_of_unity = exp(2j*pi/(size) * arange(size))
    return (1/sqrt(size)) * vander(roots_of_unity, increasing=True)

def K(U):
    r'''Method to create the Kitaev gate for phase estimation with 1 q-bit'''
    return matmul(kron(Hadamard, eye(U.shape[0])), matmul(ControlledGate(U), kron(Hadamard, eye(U.shape[0]))))

def PhaseEstimation(U, n):
    r'''Get Phase Estimation circuit with `n` q-bits as estimation'''
    if 2**(int(log2(U.shape[0]))) != U.shape[0]:
        raise TypeError("The given gate do not apply to the base case of q-bits.")
    nbits_U = int(log2(U.shape[0]))
    tot_bits = n + nbits_U

    to_apply = [kron(reduce(lambda p,q: kron(p,q), n*[Hadamard]), eye(U.shape[0]))]
    Us = [matrix_power(U, 2**i) for i in range(n)]
    to_apply.extend([ApplyWithControl(Us[i], n-1-i, [n+i for i in range(nbits_U)], tot_bits) for i in range(n)])

    to_apply.append(kron(inv(FourierTransform_direct(n)), eye(U.shape[0]))) # adding the inverse of fourier at the end

    # we now apply all matrices in inverse order
    return reduce(lambda p,q: matmul(q,p), to_apply)
