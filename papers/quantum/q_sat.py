r'''
    Script to generate data on measuring time for the 3-SAT problem

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
from math import sqrt
from numpy import count_nonzero, diag, diagonal, exp
from qiskit.circuit import QuantumCircuit, Parameter
from qiskit.circuit.library import PhaseGate
from random import randint

## Imports from the local folder
from misc import *

class Clause:
    def __init__(self, *elements):
        elements = tuple(set(tuple(el) for el in elements)) # removed repeated
        if any((i != j and elements[i][0] == elements[j][0]) for i in range(len(elements)) for j in range(len(elements))): # if two opposite variables are present we have a trivial clause 
            elements = tuple()

        self.__values = tuple([el[0] for el in elements])
        self.__pos = tuple([el[1] for el in elements])

    ### OTHER BUILDERS
    @staticmethod
    def random(max_variables: int, num_variables: int) -> Clause:
        return Clause(*tuple(set((randint(0,max_variables-1), randint(0,1)) for _ in range(num_variables))))
    
    @staticmethod
    def parse(clause: str) -> Clause:
        clause = clause.strip().removeprefix("(").removesuffix(")")
        return Clause(*tuple((int(var.split("_")[1]), int(not var.startswith("-"))) for var in (prt.strip() for prt in clause.split("v"))))
    
    def is_trivial(self) -> bool:
        return len(self) == 0
    
    ### STATIC METHODS
    @staticmethod
    def cp(circuit: QuantumCircuit, dt: float | Parameter, b: list[int]) -> QuantumCircuit:
        r'''Build the quantum gate Controlled Phase over `circuit` with parameter `dt` over the qubits in `b`'''
        c, t = b[0], b[1]
        par = dt/2

        ## Apply phase shift to both qubits
        circuit.rz(par, c)
        circuit.rz(par, t)
        ## We do a CNOT
        circuit.cx(c,t)
        ## Undo phase shift in case of need
        circuit.rz(-par, t)
        ## Undo the CNOT
        circuit.cx(c,t)

        return circuit
    
    @staticmethod
    def ccp(circuit: QuantumCircuit, dt: float | Parameter, b: list[int]) -> QuantumCircuit:
        r'''Build the quantum gate Controlled-Controlled Phase over `circuit` with parameter `dt` over the qubits in `b`'''
        s1,s2,t = b
        par = dt/4
        
        ## Apply some phase shift to all qubits
        circuit.rz(par, s1)
        circuit.rz(par, s2)
        circuit.rz(2*par, t)

        ## Apply CNOT and CCNOT
        circuit.ccx(s1, s2, t)
        circuit.cx(s1, s2)
        ## Undo phase shift in case of need
        circuit.rz(-par, s2)
        circuit.rz(-2*par, t)
        ## Undo the CNOT and CCNOT
        circuit.cx(s1, s2)
        circuit.ccx(s1, s2, t)

        return circuit

    @staticmethod
    def F(circuit: QuantumCircuit, dt: float | Parameter, b: list[int]) -> QuantumCircuit:
        r'''Build the quantum gate F over `circuit` with parameter `dt` over the qubits in `b`'''
        circuit.rz(-dt, b[0])

        return circuit
    @staticmethod
    def T(circuit: QuantumCircuit, dt: float | Parameter, b: list[int]) -> QuantumCircuit:
        r'''Build the quantum gate T over `circuit` with parameter `dt` over the qubits in `b`'''
        circuit.x(b[0])
        circuit.rz(-dt, b[0])
        circuit.x(b[0])

        return circuit
    @staticmethod
    def TT(circuit: QuantumCircuit, dt: float | Parameter, b: list[int]) -> QuantumCircuit:
        r'''Build the quantum gate TT over `circuit` with parameter `dt` over the qubits in `b`'''
        ## Negating qubit 0
        circuit.x(b[0])
        ## Apply conditioned T to qubit 2
        circuit.cx(b[0], b[1])
        Clause.cp(circuit, -dt, b)
        circuit.cx(b[0], b[1])
        ## Restoring qubit 0
        circuit.x(b[0])
        return circuit
    @staticmethod
    def FT(circuit: QuantumCircuit, dt: float | Parameter, b: list[int]) -> QuantumCircuit:
        r'''Build the quantum gate FT over `circuit` with parameter `dt` over the qubits in `b`'''
        ## Apply conditioned T to qubit 2
        circuit.cx(b[0], b[1])
        Clause.cp(circuit, -dt, b)
        circuit.cx(b[0], b[1])
        return circuit
    @staticmethod
    def FF(circuit: QuantumCircuit, dt: float | Parameter, b: list[int]) -> QuantumCircuit:
        r'''Build the quantum gate FF over `circuit` with parameter `dt` over the qubits in `b`'''
        ## Apply conditioned F to qubit 2
        Clause.cp(circuit, -dt, b)
        return circuit
    @staticmethod
    def FFF(circuit: QuantumCircuit, dt: float | Parameter, b: list[int]) -> QuantumCircuit:
        r'''Build the quantum gate FFF over `circuit` with parameter `dt` over the qubits in `b`'''
        return Clause.ccp(circuit, -dt, b)
    @staticmethod
    def FTT(circuit: QuantumCircuit, dt: float | Parameter, b: list[int]) -> QuantumCircuit:
        r'''Build the quantum gate FTT over `circuit` with parameter `dt` over the qubits in `b`'''
        circuit.cx(b[0], b[1])
        circuit.ccx(b[0], b[1], b[2])
        Clause.ccp(circuit, -dt, b)
        circuit.ccx(b[0], b[1], b[2])
        circuit.cx(b[0], b[1])
        return circuit
    @staticmethod
    def FFT(circuit: QuantumCircuit, dt: float | Parameter, b: list[int]) -> QuantumCircuit:
        r'''Build the quantum gate FFT over `circuit` with parameter `dt` over the qubits in `b`'''
        circuit.ccx(b[0],b[1],b[2])
        Clause.ccp(circuit, -dt, b)
        circuit.ccx(b[0],b[1],b[2])
        return circuit
    @staticmethod
    def TTT(circuit: QuantumCircuit, dt: float | Parameter, b: list[int]) -> QuantumCircuit:
        r'''Build the quantum gate TTT over `circuit` with parameter `dt` over the qubits in `b`'''
        circuit.x(b[0])
        circuit.cx(b[0], b[1])
        circuit.ccx(b[0], b[1], b[2])
        Clause.ccp(circuit, -dt, b)
        circuit.ccx(b[0], b[1], b[2])
        circuit.cx(b[0], b[1])
        circuit.x(b[0])
        return circuit

    ### PROPERTIES OF A CLAUSE
    @property
    def variables(self):
        return self.__values
    @property
    def terms(self):
        return zip(self.__values, self.__pos)

    ### METHODS OF A CLAUSE
    def eval(self, *values):
        return any(self.__pos[i] == values[v] for i,v in enumerate(self.__values))

    def quantum_gate(self, circuit: QuantumCircuit, dt: Parameter) -> QuantumCircuit:
        r'''Build over `circuit` the quantum gate corresponding to the clause'''
        m = len(self) # number of variables in the clause
        f = len([v for (v,p) in self.terms if p == 0]) # false variables
        b = [v for (v,p) in self.terms if p == 0] + [v for (v,p) in self.terms if p == 1]

        gate_gen = [
            None, 
            [Clause.T, Clause.F], 
            [Clause.TT,Clause.FT, Clause.FF], 
            [Clause.TTT,Clause.FTT,Clause.FFT,Clause.FFF]
        ]
        return gate_gen[m][f](circuit, dt, b)
    
    ### MAGIC METHODS
    def __len__(self) -> int:
        return len(self.variables)

    def __call__(self, *values):
        return self.eval(*values)

    def __eq__(self, other):
        if not isinstance(other, Clause):
            return False
        
        if not set(self.__values) == set(other.__values):
            return False
        
        return all(self.__pos[i] == other.__pos[other.__values.index(v)] for (i,v) in enumerate(self.__values))
    
    def __hash__(self) -> int:
        return hash(tuple([self.__values, self.__pos]))

    def __repr__(self) -> str:
        return "(" + " v ".join(f"{'-' if not cl[1] else ''}x_{cl[0]}" for cl in self.terms) + ")"

class SATFormula(set[Clause], Experiment):
    def __init__(self, *clauses: Clause, n=None):
        self.n = n
        super().__init__()
        for clause in clauses: self.add(clause)

    ### METHODS OF SET
    def add(self, __element: Clause) -> None:
        if not isinstance(__element, Clause):
            raise TypeError(f"SATFormula only allow clauses as elements")
        if self.n != None:
            if any(v >= self.n for v in __element.variables):
                raise ValueError(f"Formula bounded to {self.n} variables. Received variable outside: {__element.variables}")
        if not __element.is_trivial():
            super().add(__element)

    ### OTHER BUILDERS
    @staticmethod
    def random(max_variables: int, max_clauses: int, force_clauses: bool = True, max_per_clause: int = 3) -> SATFormula:
        formula = SATFormula(*[Clause.random(max_variables, max_per_clause) for _ in range(max_clauses)], n=max_variables)

        if force_clauses:
            while(len(formula) != max_clauses):
                formula.add(Clause.random(max_variables, max_per_clause))

        return formula
    
    @staticmethod
    def parse(formula) -> SATFormula:
        return SATFormula(
            [
                Clause.parse(clause) 
                for clause in (prt.strip() for prt in formula.strip().removeprefix("[").removesuffix("]").split("^")) 
                if len(clause) > 0
            ]
        )
    
    ### PROPERTIES OF FORMULA
    @property
    def total_size(self) -> int:
        return self.n if self.n != None else max(self.variables)+1
    @property
    def total_variables(self) -> tuple[int]:
        return tuple(range(self.total_size))
    @property
    def variables(self):
        return set(sum((list(clause.variables) for clause in self), []))
    
    ### METHODS FOR 3-SAT FORMULAS
    def value(self, *values) -> int:
        return sum(clause(*values) for clause in self)
    
    def eval(self, *values) -> bool:
        return all(clause(*values) for clause in self)
    
    def eval_matrix(self) -> SparseRowMatrix:
        r'''Create a diagonal matrix with the value of each evaluation'''
        M = SparseRowMatrix(2**self.total_size, field=CC)
        for i,values in enumerate(product(range(2), repeat=self.total_size)):
            M.increment(i,i,self.value(*values))

        return M
    
    def eval_values(self):
        M = self.eval_matrix()
        return set(M[(i,i)] for i in range(M.nrows))
    
    def eval_split(self) -> tuple[SparseRowMatrix, SparseRowMatrix]:
        r'''Return a direct lumping and reduced model from the formula'''
        eigenvalues = dict()
        for i,values in enumerate(product(range(2), repeat=self.total_size)):
            cost = self.value(*values)
            if not cost in eigenvalues:
                eigenvalues[cost] = list()
            eigenvalues[cost].append(i)

        d = len(eigenvalues)
        L = SparseRowMatrix((d,2**(self.total_size)), CC)
        U = SparseRowMatrix((d,d), CC)
        for i,eigenvalue in enumerate(eigenvalues.keys()):
            m = 1/sqrt(len(eigenvalues[eigenvalue]))
            for j in eigenvalues[eigenvalue]:
                L.increment(i,j,m)
            U.increment(i,i,exp(-1j*eigenvalue))

        return L,U

    def eval_quantum(self, circuit: QuantumCircuit = None, par: Parameter = None) -> tuple[QuantumCircuit, Parameter]:
        if circuit is None:
            par = Parameter("t")
            circuit = QuantumCircuit(len(self.total_variables))

        for clause in self:
            circuit = clause.quantum_gate(circuit, par)

        return circuit, par

    def eval_quantumB(self, circuit: QuantumCircuit = None, par: Parameter = None) -> tuple[QuantumCircuit, Parameter]:
        r'''
            From :arxiv:`0001106v1`, the begin Hamiltonian is the exponential matrix 

            .. MATH::

                H_B = e^{-itB}

            where `B = \sum_{j=1}^n d_j B_j`. First of all, the matrices `B_j` commute among themselves, since they only 
            act on one of the q-bits. Hence, we can compute the matrix `H_B` as the product of the exponential of the independent
            `B_j` with the corresponding factor `-id_j t`.
             
            Now, `B_j` diagonalizes with eigenvalues 1 and 0 on the basis `\ket{+},\ket{-}`. Hence we can deduce that 

            .. MATH::

                e^{-id_j t B_j} = H_j\begin{pmatrix}e^{-id_j t} & 0\\0 & 1\ned{pmatrix}H_j = H_j \sigma_j^x P_j(-d_j t) \sigma_j^x H_j.
        '''
        if circuit is None:
            par = Parameter("t")
            circuit = QuantumCircuit(len(self.total_variables))

        count = [sum(n in clause.variables for clause in self) for n in self.total_variables]
        for i in range(len(self.total_variables)):
            circuit.h(i)
            circuit.x(i)
            circuit.rz(-par*count[i], i)
            circuit.x(i)
            circuit.h(i)
        
        return circuit, par

    def eval_begin_matrix(self) -> SparseRowMatrix:
        count = [sum(n in clause.variables for clause in self) for n in self.total_variables]
        M = SparseRowMatrix(2**self.total_size, field=CC)
        for index in range(self.total_size):
            for values in product(range(2), repeat=(self.size()-1)):
                values = list(values)
                i = int("".join([str(v) for v in values[:index]] + ['0'] + [str(v) for v in values[index:]]), 2)
                j = int("".join([str(v) for v in values[:index]] + ['1'] + [str(v) for v in values[index:]]), 2)
                M.increment(i,i,count[index]*0.5)
                M.increment(i,j,count[index]*(-0.5))
                M.increment(j,i,count[index]*(-0.5))
                M.increment(j,j,count[index]*0.5)

        return M
    
    def eval_begin_quantum(self) -> tuple[QuantumCircuit, Parameter]:
        raise NotImplementedError("QUantum circuit for begin matrix not yet implemented")
    
    ### MAGIC METHODS
    def __call__(self, *values) -> bool:
        return self.eval(*values)
    
    def __repr__(self) -> str:
        return "["+ " ^ ".join(repr(clause) for clause in self) +"]"

    ### OTHER STATIC METHODS
    @staticmethod
    def store_circuit(formula: SATFormula, parameter, name: str="formula"):
        circuit, par = formula.eval_quantum()
        circuit = circuit.assign_parameters({par: parameter})

        n = formula.total_size; m = len(formula)
        name = f"{name}_{n}_{m}"
        final_name = name; i = 0
        while os.path.exists(os.path.join(SCRIPT_DIR, "hamiltonians", f"{final_name}.qasm")):
            final_name = f"{name}[{i}]"
            i += 1
        circuit.qasm(True, os.path.join(SCRIPT_DIR, "hamiltonians", f"{final_name}.qasm"))
        with open(os.path.join(SCRIPT_DIR, "hamiltonians", f"{final_name}.qasm"), "a+") as f:
            f.write(f"\n\n// True formula: {formula}\n")
            variables = formula.variables
            f.write(f"// Variables appearing: {variables}\n")
            f.write(f"// Appear all: {len(variables) == n}\n")
            f.write(f"\n// For creating the entangled state:\n")
            for i in range(n):
                f.write(f"//h q[{i}];\n")
            f.write(f"//barrier q;")

    ## NECESSARY METHOD FOR EXPERIMENT
    def size(self) -> int: return self.total_size
    def correct_size(self) -> int: return len(self.eval_values())
    def direct(self) -> tuple[SparseRowMatrix, SparseRowMatrix]: return self.eval_split()
    def matrix(self) -> SparseRowMatrix: return self.eval_matrix()
    def matrix_B(self, red_U: ndarray) -> ndarray: 
        if count_nonzero(red_U - diag(diagonal(red_U))): # not diagonal
            return diag([len(self)] + (red_U.shape[0]-1)*[1/len(self)])
        else:
            raise NotImplementedError(f"[full-clue] Base hamiltonian not defined when U_P is diagonal")
    def quantum(self, circuit: QuantumCircuit, par: Parameter) -> tuple[QuantumCircuit, Parameter]: return self.eval_quantum(circuit, par)
    def quantum_B(self, circuit: QuantumCircuit, par: Parameter) -> tuple[QuantumCircuit, Parameter]: return self.eval_quantumB(circuit, par)
    def data(self): return [len(self)]

    @staticmethod
    def generate_example(_:str, size: int) -> SATFormula:
        formula = SATFormula.random(size, randint(size, 3*size))
        while len(formula.variables) != formula.total_size:
            formula = SATFormula.random(size, randint(size, 3*size))
        
        return formula

    @staticmethod
    def generate_header(csv_writer, ttype):
        if ttype in ("clue", "ddsim", "direct"):
            csv_writer.writerow(["size", "clauses", "red. ratio", "time_lumping", "memory (MB)", "formula"])
        elif ttype in ("full_clue", "full_direct"):
            csv_writer.writerow(["size", "clauses", "time_lumping", "kappa", "time_iteration", "memory (MB)", "formula"])
        elif ttype == "full_ddsim":
            csv_writer.writerow(["size", "clauses", "kappa", "time_iteration", "memory (MB)", "formula"])
        elif ttype == "full_quokka#":
            csv_writer.writerow(["size", "clauses", "kappa", "time_encoding", "time_iteration", "memory (MB)", "formula"])
        else:
            raise NotImplementedError(f"Type of file {ttype} not recognized")

    ## METHODS TO GENERATE OBSERVABLES
    @staticmethod
    def generate_observable_clue(formula: SATFormula, *_) -> tuple[SparseVector]:
        return tuple([SparseVector.from_list(2**formula.size()*[1], field=CC)])

    @staticmethod
    def generate_observable_ddsim(formula: SATFormula, *_) -> bool:
        return bool(formula)

def generate_example(_:str, size: int) -> SATFormula:
    return SATFormula.generate_example(_, size)

def generate_header(csv_writer, ttype):
    return SATFormula.generate_header(csv_writer, ttype)

## METHODS TO GENERATE OBSERVABLES
def generate_observable_clue(formula: SATFormula, *_) -> tuple[SparseVector]:
    return SATFormula.generate_observable_clue(formula, *_)

def generate_observable_ddsim(formula: SATFormula, *_) -> bool:
    return SATFormula.generate_observable_ddsim(formula, *_)

if __name__ == "__main__":
    ## Processing the arguments
    ttype, script = get_method(*sys.argv)
    m, M = get_size_bounds(*sys.argv)
    timeout = get_timeout(*sys.argv)
    repeats = get_repeats(*sys.argv)
    name = "sat"; obs = None
    
    main_script(SCRIPT_DIR, SCRIPT_NAME, name,
                [generate_header, generate_example, generate_observable_clue, generate_observable_ddsim],
                ttype, script,               
                m, M,                        
                timeout, repeats,
                obs)