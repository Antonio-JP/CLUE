r'''
    Script to generate data on measuring time for the MaxCut problem

    This is used for generating data on Table 1 on the Quantum draft (see :arxiv:`2308.09510`).
'''
from __future__ import annotations
import sys, os

SCRIPT_DIR = os.path.dirname(__file__) if __name__ != "__main__" else "./"
sys.path.insert(0, os.path.join(SCRIPT_DIR, "..", "..")) # clue is here

import tracemalloc
from collections import defaultdict
from clue import FODESystem
from clue.linalg import CC, NumericalSubspace, SparseRowMatrix, SparseVector
from csv import writer
from itertools import product
from mqt import ddsim
from qiskit import execute
from qiskit.circuit import QuantumCircuit, Parameter
from random import sample, random
from time import time

## Imports from the local folder
from misc import *

class UndirectedGraph(defaultdict):
    def __init__(self):
        super().__init__()

        self.__vertices = []

    @staticmethod
    def random(vertices: int, edges:int = None, density:float = None):
        ## We create the vertices
        G = UndirectedGraph()
        for _ in range(vertices): G.add_vertex()

        possible_edges = [edge for edge in product(G, repeat=2) if edge[0] != edge[1]]
        possible_edges = list(set(tuple(sorted(list(edge))) for edge in possible_edges))

        ## If we are given number of edges, we do the pick randomly
        if edges != None:
            edges = sample(possible_edges, edges)
        else: # Otherwise we do with the density
            edges = [edge for edge in possible_edges if random() < density]

        for (s,t) in edges:
            G.add_edge(s,t)

        return G

    def __setitem__(self, __key, __value) -> None:
        raise NotImplementedError
    
    def __delitem__(self, __key) -> None:
        raise NotImplementedError

    def add_vertex(self, label=None):
        if label is None:
            label = 0
            while label in self:
                label += 1
        if label in self:
            raise ValueError(f"Repeated labels are not allowed")
        self.__vertices.append(label)
        super().__setitem__(label, set())

    def add_edge(self, src, trg):
        if src not in self:
            raise KeyError(f"Vertex {src} is nto in the graph.")
        if trg not in self:
            raise KeyError(f"Vertex {trg} is nto in the graph.")
        if trg in self[src]:
            raise ValueError(f"Repeated edges are not allow")
        
        self[src].add(trg); self[trg].add(src)

    @property
    def vertices(self):
        return self.__vertices.copy()

    @property
    def edges(self):
        return list(set(tuple(set([v,w])) for v in self for w in self[v]))

    def adjacency_matrix(self):
        M = SparseRowMatrix(len(self))
        for (i,v) in enumerate(self.__vertices):
            for w in self[v]:
                M.increment(i, self.__vertices.index(w), 1)
        return M
        
    def cut_cost(self, V_0: set, V_1 = None):
        ## Checking V_0
        if any(v not in self for v in V_0):
            raise ValueError(f"The given set {V_0=} contains elements not in the grpah")
        
        ## Checking V_1
        if V_1 != None:
            if V_0.union(V_1) != set(self.__vertices):
                raise ValueError(f"The given cut does not cover the whole graph")
        else:
            V_1 = set([v for v in self if v not in V_0])

        ## Counting edges from V_0 to V_1
        return sum(len(self[v].intersection(V_1)) for v in V_0)

    def cut_matrix(self):
        r'''Return a 2^n square matrix with the cut counting'''
        M = SparseRowMatrix(2**len(self), field=CC)
        for i,values in enumerate(product(range(2), repeat=len(self))):
            M.increment(i,i,self.cut_cost(set(v for (j,v) in enumerate(self.__vertices) if values[j] == 0)))

        return M

    def cut_values(self) -> set[int]:
        M = self.cut_matrix()
        return set(M[(i,i)] for i in range(M.nrows))

    def visualize(self, save=False):
        import networkx as nx
        import matplotlib.pyplot as plt
        G = nx.Graph()
        G.add_edges_from(self.edges)
        nx.draw_networkx(G)
        if save:
            plt.savefig(save)
        else:
            plt.show()
        plt.close()

    def quantum_cut(self):
        circuit = QuantumCircuit(len(self))
        t = Parameter("t")
        edge_gate = QuantumCircuit(2, name="E")
        edge_gate.x(0)
        edge_gate.cp(t, 0, 1)
        edge_gate.x(0); edge_gate.x(1)
        edge_gate.cp(t, 1, 0)
        edge_gate.x(1)

        for edge in self.edges:
            circuit.append(edge_gate, edge)

        return circuit, t
            
    @staticmethod
    def store_circuit(graph: UndirectedGraph, parameter, name="graph"):
        circuit, par = graph.quantum_cut()
        circuit = circuit.bind_parameters({par: parameter})

        name = f"{name}_{len(graph)}_{len(graph.edges)}"
        final_name = name; i = 0
        while os.path.exists(os.path.join(SCRIPT_DIR, "graphs", f"{final_name}.qasm")):
            final_name = f"{name}_{i}"
            i += 1
        circuit.qasm(True, os.path.join(SCRIPT_DIR, "graphs", f"{final_name}.qasm"))
        with open(os.path.join(SCRIPT_DIR, "graphs", f"{final_name}.qasm"), "a+") as f:
            f.write(f"\n\n// Description of the graph:\n")
            f.write(f"//\t * Vertices: {graph.vertices}\n")
            f.write(f"//\t * Edges: {graph.edges}\n")
            f.write(f"//\t * Number of edges: {len(graph.edges)}\n\n")
            f.write(f"\n// For creating the entangled state:\n")
            for i in range(len(graph)):
                f.write(f"//h q[{i}];\n")
            f.write(f"//barrier q;")
        graph.visualize(os.path.join(SCRIPT_DIR, "graphs", f"{final_name}.png"))

    def __repr__(self):
        return f"(Vertices: {self.vertices}) -- (Edges: {self.edges})"

def generate_valid_example(size: int) -> UndirectedGraph:
    print(f"%%% [GEN] Generating a valid graph with {size} nodes...")
    graph = UndirectedGraph.random(size, density=1/3)    
    while(len(graph.edges) == 0):
        graph = UndirectedGraph.random(size, density=1/3)    
    print(f"%%% [GEN] Generated a graph with {len(graph.edges)} edges")
    return graph

def gen_header(csv_writer, ttype):
    if ttype in ("clue", "ddsim"):
        csv_writer.writerow(["size", "edges", "red. ratio", "time_lumping", "memory (MB)", "graph"])
    elif ttype == "full_clue":
        csv_writer.writerow(["size", "edges", "time_lumping", "kappa", "time_iteration", "memory (MB)", "graph"])
    elif ttype == "full_ddsim":
        csv_writer.writerow(["size", "edges", "kappa", "time_iteration", "memory (MB)", "graph"])
    else:
        raise NotImplementedError(f"Type of file {ttype} not recognized")

def clue_reduction(size: int, result_file): 
    r'''
        This method computes the CLUE lumping

        This method generates a random 3-SAT formula and computes the lumping of the problem matrix associated
        to it with CLUE.
         
        It stores the execution time of the lumping plus the reduction ratio. It stores the result on ``result_file``.
        ["size", "edges", "red. ratio", "time_lumping", "memory", "graph"]
    '''
    graph = generate_valid_example(size)
    print(f"%%% [clue] Computing all cut values for graph...")
    eval_values = graph.cut_values()

    print(f"%%% [clue] Creating the full system with the problem matrix")
    system = FODESystem.LinearSystem(graph.cut_matrix(), lumping_subspace=NumericalSubspace)
    obs = tuple([SparseVector.from_list(system.size*[1], field=CC)])
    
    print(f"%%% [clue] Computing the lumped system...")
    tracemalloc.start()
    ctime = time()
    lumped = system.lumping(obs, print_reduction=False, print_system=False)
    ctime = time()-ctime
    memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(f"%%% [clue] Storing the data...")
    if len(eval_values) != lumped.size:
        print(f"%%% [clue] ERROR!! Found weird dimension in lumping -- \n%%% \t* Expected: {len(eval_values)}\n%%% \t* Got: {lumped.size}\n%%% \t* Graph: {graph}")
    result_file.writerow([size, len(graph.edges), lumped.size/system.size, ctime, memory, repr(graph)])

def ddsim_reduction(size: int, result_file): 
    r'''
        This method computes the DDSIM lumping

        This method generates a random 3-SAT formula and computes the lumping of the problem matrix associated
        to it with DDSIM. Since we can not use linear algebra on Decision Diagrams, we compute the simulation time
        of a quantum circuit that repeats the circuit as many times as the bound of the lumping size (i.e., number
        of clauses).
         
        It stores the execution time of the lumping plus the reduction ratio. It stores the result on ``result_file``.
        ["size", "edges", "red. ratio", "time_lumping", "memory", "graph"]
    '''
    graph = generate_valid_example(size)
    print(f"%%% [ddsim] Computing all cut values for graph...")
    eval_values = graph.cut_values() if size <=20 else list(range(len(graph)))

    print(f"%%% [ddsim] Creating the full cirtuit and job to simulate with DDSIM")
    m = len(eval_values) # number of values of true simulataneously
    U_P, par = graph.quantum_cut(); U_P = U_P.bind_parameters({par: 1/(m*1000)})
    circuit = loop(U_P, size, m, True, True)
    backend = ddsim.DDSIMProvider().get_backend("qasm_simulator")
    
    print(f"%%% [ddsim] Computing the simulation of the circuit...")
    tracemalloc.start()
    ctime = time()
    ## Executing the circuit one time
    job = execute(circuit, backend, shots=1)
    job.result()
    ctime = time()-ctime
    memory = tracemalloc.get_traced_memory()[1] / (2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(f"%%% [ddsim] Storing the data...")
    result_file.writerow([size, len(graph.edges), m/2**size, ctime, memory, repr(graph)])

def clue_iteration(size: int, iterations, result_file): 
    r'''
        This method computes the CLUE iteration

        This method generates a random 3-SAT formula and computes the lumping of the problem matrix associated
        to it with CLUE. 
        
        Then it creates a valid begin Hamiltonian within the invariant space and compute the alternate application 
        of both the begin and problem Hamiltonian in the reduced space (i.e., we used the lumped matrix from the 
        actual lumping) up to ``iterations`` times. 
         
        It stores the execution time of the lumping, the number of iterations and the time for the computed iteration. 
        It stores the result on ``result_file``.
        ["size", "edges", "time_lumping", "kappa", "time_iteration", "memory", "graph"]
    '''
    graph = generate_valid_example(size)
    raise NotImplementedError(f"Experiment 'full_clue' not yet implemented")

def ddsim_iteration(size: int, iterations, result_file): 
    r'''
        This method computes the DDSIM iteration

        This method generates a random 3-SAT formula and computes application of ``iterations`` times the alternating 
        of the begin and problem Hamiltonian associated with the formula. 
        The method will build both circuits, alternate them and combine them on a loop of ``iterations`` times. 
         
        It stores the number of iterations and the time for the computed iteration. 
        It stores the result on ``result_file``.
        ["size", "edges", "kappa", "time_iteration", "memory", "graph"]
    '''
    graph = generate_valid_example(size)
    raise NotImplementedError(f"Experiment 'full_ddsim' not yet implemented")

if __name__ == "__main__":
    n = 1; m = 3; M=10; ttype="clue"; repeats=100
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
            elif sys.argv[n].endswith("r"):
                repeats = int(sys.argv[n+1]); n+=2
        else:
            n += 1

    methods = [clue_reduction, ddsim_reduction, clue_iteration, ddsim_iteration]
    method = methods[["clue", "ddsim", "full_clue", "full_ddsim"].index(ttype)]
    existed = os.path.exists(os.path.join(SCRIPT_DIR, f"[result]q_maxcut_{ttype}.csv"))
    with open(os.path.join(SCRIPT_DIR, f"[result]q_maxcut_{ttype}.csv"), "at" if existed else "wt") as result_file:
        csv_writer = writer(result_file)
        if not existed:
            gen_header(csv_writer, ttype)
        print(f"##################################################################################")
        print(f"### EXECUTION ON MAXCUT [{m=}, {M=}, {repeats=}, method={ttype}]")
        print(f"##################################################################################")
        for size in range(m, M+1):
            for execution in range(1,repeats+1):
                print(f"### Starting execution {execution}/{repeats} ({size=})")
                if ttype in ("clue", "ddsim"):
                    method(size, csv_writer)
                else:
                    for it in (10,100,1000):#,10000)
                        print(f"------ Case with {it} iterations")
                        method(size, it, csv_writer)
                print(f"### Finished execution {execution}/{repeats}")
                result_file.flush()