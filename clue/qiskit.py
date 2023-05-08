r'''
    Module implementing a type of dynamical systems that are based on Quantum Circuits

    This module is somehow a parser that allows to create a dynamical system from 
    the unitary matrix that defines a quantum circuit. This uses the interface by 
    :mod:`qiskit` that allows to manipulate these circuits and generate the 
    corresponding unitary matrices.
'''

from itertools import product
from qiskit import execute, QuantumCircuit, Aer
from numpy import asarray
from sympy import CC

from .clue import FODESystem
from .linalg import NumericalSubspace
from .rational_function import SparsePolynomial

class DS_QuantumCircuit(FODESystem):
    BACKEND = Aer.get_backend("unitary_simulator")

    def __init__(self, path, **kwds):
        circuit = QuantumCircuit.from_qasm_file(path)
        job = execute(circuit, DS_QuantumCircuit.BACKEND, shots=kwds.pop("shots", 8192))
        unitary = job.result().get_unitary(circuit)

        nbits = len(unitary.input_dims())
        variables = [f"Q_{f'{i:b}'.ljust(nbits, '0')}" for i in range(2**nbits)]
        nvariables = 2**nbits
        matr = asarray(unitary)
        equations = [SparsePolynomial.from_vector(matr[i], variables, CC) for i in range(nvariables)]
        
        super().__init__(
            equations, variables=variables, 
            name=circuit.name, 
            lumping_subspace=NumericalSubspace, 
            lumping_subspace_kwds=kwds
        )
