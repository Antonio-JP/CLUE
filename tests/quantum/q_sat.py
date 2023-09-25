r'''
    Script to generate data on measuring time for the 3-SAT problem

    This is used for generating data on Table 1 on the Quantum draft (see :arxiv:`2308.09510`).
'''
from __future__ import annotations
import sys, os

SCRIPT_DIR = os.path.dirname(__file__) if __name__ != "__main__" else "./"
sys.path.insert(0, os.path.join(SCRIPT_DIR, "..", "..")) # clue is here

import tracemalloc
from clue import FODESystem
from clue.linalg import CC, NumericalSubspace, SparseRowMatrix, SparseVector
from csv import writer
from functools import lru_cache
from itertools import product
from mqt import ddsim #pylint: disable=no-name-in-module
from numpy import cdouble, count_nonzero, diag, diagonal, matmul, Inf
from numpy.linalg import matrix_power
from qiskit import execute
from qiskit.circuit import QuantumCircuit, Parameter
from qiskit.circuit.library import PhaseGate
from random import randint
from time import time

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
    @lru_cache(maxsize=16)
    def F(dt) -> QuantumCircuit:
        C = QuantumCircuit(1, name=f"F")
        C.p(-dt,0)
        return C
    @staticmethod
    @lru_cache(maxsize=16)
    def T(dt) -> QuantumCircuit:
        C = QuantumCircuit(1, name=f"T")
        C.x(0)
        C.p(-dt,0)
        C.x(0)
        return C
    @staticmethod
    @lru_cache(maxsize=16)
    def TT(dt) -> QuantumCircuit:
        C = QuantumCircuit(2, name=f"TT")
        ## Negating qubit 0
        C.x(0)
        ## Apply conditioned T to qubit 2
        C.cx(0, 1)
        C.cp(-dt, 0, 1)
        C.cx(0, 1)
        ## Restoring qubit 0
        C.x(0)
        return C
    @staticmethod
    @lru_cache(maxsize=16)
    def FT(dt) -> QuantumCircuit:
        C = QuantumCircuit(2, name=f"FT")
        ## Apply conditioned T to qubit 2
        C.cx(0, 1)
        C.cp(-dt, 0, 1)
        C.cx(0, 1)
        return C
    @staticmethod
    @lru_cache(maxsize=16)
    def FF(dt) -> QuantumCircuit:
        C = QuantumCircuit(2, name=f"FF")
        ## Apply conditioned F to qubit 2
        C.cp(-dt, 0, 1)
        return C
    @staticmethod
    @lru_cache(maxsize=16)
    def FFF(dt) -> QuantumCircuit:
        C = QuantumCircuit(3, name=f"FFF")
        C.append(PhaseGate(-dt).control(2), [0,1,2])
        
        return C
    @staticmethod
    @lru_cache(maxsize=16)
    def FTT(dt) -> QuantumCircuit:
        C = QuantumCircuit(3, name=f"FTT")
        C.cx(0,1)
        C.ccx(0,1,2)
        C.append(PhaseGate(-dt).control(2), [0,1,2])
        C.ccx(0,1,2)
        C.cx(0,1)
        return C
    @staticmethod
    @lru_cache(maxsize=16)
    def FFT(dt) -> QuantumCircuit:
        C = QuantumCircuit(3, name=f"FFT")
        C.ccx(0,1,2)
        C.append(PhaseGate(-dt).control(2), [0,1,2])
        C.ccx(0,1,2)
        return C
    @staticmethod
    @lru_cache(maxsize=16)
    def TTT(dt) -> QuantumCircuit:
        C = QuantumCircuit(3, name=f"TTT")
        C.x(0)
        C.cx(0,1)
        C.ccx(0,1,2)
        C.append(PhaseGate(-dt).control(2), [0,1,2])
        C.ccx(0,1,2)
        C.cx(0,1)
        C.x(0)
        return C

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

    def quantum_gate(self, dt) -> tuple[QuantumCircuit, tuple[int]]:
        m = len(self) # number of variables in the clause
        f = len([v for (v,p) in self.terms if p == 0]) # false variables
        b = [v for (v,p) in self.terms if p == 0] + [v for (v,p) in self.terms if p == 1]
        
        gate_gen = [
            None, 
            [Clause.T, Clause.F], 
            [Clause.TT,Clause.FT, Clause.FF], 
            [Clause.TTT,Clause.FTT,Clause.FFT,Clause.FFF]
        ]
        return (gate_gen[m][f](dt), b)
        
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

class SATFormula(set[Clause]):
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

    def eval_quantum(self) -> tuple[QuantumCircuit, Parameter]:
        par = Parameter("t")
        circuit = QuantumCircuit(len(self.total_variables))
        gates = [clause.quantum_gate(par) for clause in self]
        trotter(circuit, gates, 2)

        return circuit, par
    
    def eval_quantumB(self) -> tuple[QuantumCircuit, Parameter]:
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
        par = Parameter("t")
        count = [sum(n in clause.variables for clause in self) for n in self.total_variables]
        circuit = QuantumCircuit(len(self.total_variables))
        for i in range(len(self.total_variables)):
            circuit.h(i)
            circuit.x(i)
            circuit.p(-par*count[i], i)
            circuit.x(i)
            circuit.h(i)
        
        return circuit, par

    def eval_begin_matrix(self) -> SparseRowMatrix:
        count = [sum(n in clause.variables for clause in self) for n in self.total_variables]
        M = SparseRowMatrix(2**self.total_size, field=CC)
        for index in range(self.total_size):
            for values in product(range(2), repeat=(n-1)):
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
        circuit = circuit.bind_parameters({par: parameter})

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

def generate_valid_example(size: int) -> SATFormula:
    print(f"%%% [GEN] Generating a valid formula of SAT3 with {size} variables...")
    formula = SATFormula.random(size, randint(size, 3*size))
    while len(formula.variables) != formula.total_size:
        print(f"%%% [GEN] Generated formula with fewer variables than expected. Repeating...")
        formula = SATFormula.random(size, randint(size, 3*size))
    
    print(f"%%% [GEN] Generated formula with {len(formula)} clauses")
    return formula

def gen_header(csv_writer, ttype):
    if ttype in ("clue", "ddsim"):
        csv_writer.writerow(["size", "clauses", "red. ratio", "time_lumping", "memory (MB)", "formula"])
    elif ttype == "full_clue":
        csv_writer.writerow(["size", "clauses", "time_lumping", "kappa", "time_iteration", "memory (MB)", "formula"])
    elif ttype == "full_ddsim":
        csv_writer.writerow(["size", "clauses", "kappa", "time_iteration", "memory (MB)", "formula"])
    else:
        raise NotImplementedError(f"Type of file {ttype} not recognized")

def clue_reduction(size: int, result_file, timeout=0): 
    r'''
        This method computes the CLUE lumping

        This method generates a random 3-SAT formula and computes the lumping of the problem matrix associated
        to it with CLUE.
         
        It stores the execution time of the lumping plus the reduction ratio. It stores the result on ``result_file``.
        ["size", "clauses", "red. ratio", "time_lumping", "memory", "formula"]
    '''
    formula = generate_valid_example(size)
    print(f"%%% [clue] Computing all true values for formula...")
    eval_values = formula.eval_values()

    print(f"%%% [clue] Creating the full system with the problem matrix")
    system = FODESystem.LinearSystem(formula.eval_matrix(), lumping_subspace=NumericalSubspace)
    obs = tuple([SparseVector.from_list(system.size*[1], field=CC)])
    
    print(f"%%% [clue] Computing the lumped system...")
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            ctime = time()
            lumped = system.lumping(obs, print_reduction=False, print_system=False)
            ctime = time()-ctime
    except TimeoutError:
        print(f"%%% [clue] Timeout reached for execution")
        ctime = Inf
    
    memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(f"%%% [clue] Storing the data...")
    if len(eval_values) != lumped.size:
        print(f"%%% [clue] ERROR!! Found weird dimension in lumping -- \n%%% \t* Expected: {len(eval_values)}\n%%% \t* Got: {lumped.size}\n%%% \t* Formula: {formula}")
    result_file.writerow([size, len(formula), lumped.size/system.size, ctime, memory, repr(formula)])

def ddsim_reduction(size: int, result_file, timeout=0): 
    r'''
        This method computes the DDSIM lumping

        This method generates a random 3-SAT formula and computes the lumping of the problem matrix associated
        to it with DDSIM. Since we can not use linear algebra on Decision Diagrams, we compute the simulation time
        of a quantum circuit that repeats the circuit as many times as the bound of the lumping size (i.e., number
        of clauses).
         
        It stores the execution time of the lumping plus the reduction ratio. It stores the result on ``result_file``.
        ["size", "clauses", "red. ratio", "time_lumping", "memory", "formula"]
    '''
    formula = generate_valid_example(size)
    print(f"%%% [ddsim] Computing all true values for formula...")
    eval_values = formula.eval_values() if size <=20 else list(range(len(formula)))

    print(f"%%% [ddsim] Creating the full circuit and job to simulate with DDSIM")
    m = len(eval_values) # number of values of true simultaneously
    U_P, par = formula.eval_quantum(); U_P = U_P.bind_parameters({par: 1/(m*1000)})
    circuit = loop(U_P, size, m, True, True)
    backend = ddsim.DDSIMProvider().get_backend("qasm_simulator")
    
    print(f"%%% [ddsim] Computing the simulation of the circuit...")
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            ctime = time()
            ## Executing the circuit one time
            job = execute(circuit, backend, shots=1)
            job.result()
            ctime = time()-ctime
    except TimeoutError:
        print(f"%%% [ddsim] Timeout reached for execution")
        ctime = Inf
    memory = tracemalloc.get_traced_memory()[1] / (2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(f"%%% [ddsim] Storing the data...")
    result_file.writerow([size, len(formula), m/2**size, ctime, memory, repr(formula)])

def clue_iteration(size: int, iterations, result_file, timeout=0): 
    r'''
        This method computes the CLUE iteration

        This method generates a random 3-SAT formula and computes the lumping of the problem matrix associated
        to it with CLUE. 
        
        Then it creates a valid begin Hamiltonian within the invariant space and compute the alternate application 
        of both the begin and problem Hamiltonian in the reduced space (i.e., we used the lumped matrix from the 
        actual lumping) up to ``iterations`` times. 
         
        It stores the execution time of the lumping, the number of iterations and the time for the computed iteration. 
        It stores the result on ``result_file``.
        ["size", "clauses", "time_lumping", "kappa", "time_iteration", "memory", "formula"]
    '''
    formula = generate_valid_example(size)
    
    print(f"%%% [full-clue] Computing all true values for formula...")
    eval_values = formula.eval_values()

    print(f"%%% [full-clue] Creating the full system with the problem matrix")
    system = FODESystem.LinearSystem(formula.eval_matrix(), lumping_subspace=NumericalSubspace)
    obs = tuple([SparseVector.from_list(system.size*[1], field=CC)])
    
    print(f"%%% [full-clue] Computing the lumped system...")
    tracemalloc.start()    
    try:
        with(Timeout(timeout)):
            lump_time = time()
            ## Executing the circuit one time
            lumped = system.lumping(obs, print_reduction=False, print_system=False)
            lump_time = time()-lump_time
    except TimeoutError:
        print(f"%%% [full-clue] Timeout reached for execution")
        lump_time = Inf
        memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
        tracemalloc.stop()
        result_file.writerow([size, len(formula), lump_time, iterations, Inf, memory, repr(formula)])
        return

    print(f"%%% [full-clue] Checking size...")
    if len(eval_values) != lumped.size:
        print(f"%%% [full-clue] ERROR!! Found weird dimension in lumping -- \n%%% \t* Expected: {len(eval_values)}\n%%% \t* Got: {lumped.size}\n%%% \t* Formula: {formula}")

    print(f"%%% [full-clue] Getting the reduced U_P")
    U_P = lumped.construct_matrices("polynomial")[0].to_numpy(dtype=cdouble)
    print(f"%%% [full-clue] Getting the reduced U_B")
    if count_nonzero(U_P - diag(diagonal(U_P))): # not diagonal
        U_B = diag([len(formula)] + (U_P.shape[0]-1)*[1/len(formula)])
    else:
        raise NotImplementedError(f"[full-clue] Base hamiltonian not defined when U_P is diagonal")
    
    print(f"%%% [full-clue] Computing the iteration (U_P*U_B)^iterations")
    it_time = time()
    _ = matrix_power(matmul(U_P, U_B), iterations)
    it_time = time() - it_time
    memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
    tracemalloc.stop()

    ## We check if the matrix is diagonal
    result_file.writerow([size, len(formula), lump_time, iterations, it_time, memory, repr(formula)])

def ddsim_iteration(size: int, iterations, result_file, timeout=0): 
    r'''
        This method computes the DDSIM iteration

        This method generates a random 3-SAT formula and computes application of ``iterations`` times the alternating 
        of the begin and problem Hamiltonian associated with the formula. 
        The method will build both circuits, alternate them and combine them on a loop of ``iterations`` times. 
         
        It stores the number of iterations and the time for the computed iteration. 
        It stores the result on ``result_file``.
        ["size", "clauses", "kappa", "time_iteration", "memory", "formula"]
    '''
    formula = generate_valid_example(size)
    print(f"%%% [full-ddsim] Computing all true values for formula...")
    eval_values = formula.eval_values() if size <=20 else list(range(len(formula)))

    print(f"%%% [full-ddsim] Creating the full circuit and job to simulate with DDSIM")
    m = len(eval_values) # number of values of true simultaneously
    U_P, par = formula.eval_quantum(); U_P = U_P.bind_parameters({par: 1/(m*10*iterations)})
    U_B, par = formula.eval_quantumB(); U_B = U_B.bind_parameters({par: 1/(m*10*iterations)})
    U_B.append(U_P, U_B.qregs[0]) # Now U_B is the alternate circuit U_B * U_P
    circuit = loop(U_B, size, iterations, True, True)
    backend = ddsim.DDSIMProvider().get_backend("qasm_simulator")
    
    print(f"%%% [full-ddsim] Computing the simulation of the circuit...")
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            ctime = time()
            ## Executing the circuit one time
            job = execute(circuit, backend, shots=1)
            job.result()
            ctime = time()-ctime
    except TimeoutError:
        print(f"%%% [full-ddsim] Timeout reached for execution")
        ctime = Inf
    memory = tracemalloc.get_traced_memory()[1] / (2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(f"%%% [full-ddsim] Storing the data...")
    result_file.writerow([size, len(formula), iterations, ctime, memory, repr(formula)])

if __name__ == "__main__":
    n = 1; m = 5; M=10; ttype="clue"; repeats=100; timeout=0
    ## Processing arguments
    while n < len(sys.argv):
        if sys.argv[n].startswith("-"):
            if sys.argv[n].endswith("m"):
                m = int(sys.argv[n+1]); n+=2
            elif sys.argv[n].endswith("M"):
                M = int(sys.argv[n+1]); n+=2
            elif sys.argv[n].endswith("t"):
                ttype = sys.argv[n+1] if sys.argv[n+1] in ("clue", "ddsim", "full_clue", "full_ddsim") else ttype
                n += 2
            elif sys.argv[n].endswith("to"):
                timeout = int(sys.argv[n+1]); n+=2
            elif sys.argv[n].endswith("r"):
                repeats = int(sys.argv[n+1]); n+=2
        else:
            n += 1

    methods = [clue_reduction, ddsim_reduction, clue_iteration, ddsim_iteration]
    method = methods[["clue", "ddsim", "full_clue", "full_ddsim"].index(ttype)]
    existed = os.path.exists(os.path.join(SCRIPT_DIR, f"[result]q_sat_{ttype}.csv"))
    with open(os.path.join(SCRIPT_DIR, f"[result]q_sat_{ttype}.csv"), "at" if existed else "wt") as result_file:
        csv_writer = writer(result_file)
        if not existed:
            gen_header(csv_writer, ttype)
        print(f"##################################################################################")
        print(f"### EXECUTION ON Q-SAT [{m=}, {M=}, {repeats=}, method={ttype}]")
        print(f"##################################################################################")
        for size in range(m, M+1):
            for execution in range(1,repeats+1):
                print(f"### Starting execution {execution}/{repeats} ({size=})")
                if ttype in ("clue", "ddsim"):
                    method(size, csv_writer, timeout=timeout)
                else:
                    for it in (1,10,100):#,1000):#,10000)
                        print(f"------ Case with {it} iterations")
                        method(size, it, csv_writer)
                print(f"### Finished execution {execution}/{repeats}")
                result_file.flush()