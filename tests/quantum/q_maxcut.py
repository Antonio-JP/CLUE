r'''
    Script to generate data on measuring time for the MaxCut problem

    This is used for generating data on Table 1 on the Quantum draft (see :arxiv:`2308.09510`).
'''
from __future__ import annotations
import sys, os

SCRIPT_DIR = os.path.dirname(__file__) if __name__ != "__main__" else "./"
sys.path.insert(0, os.path.join(SCRIPT_DIR, "..", "..")) # clue is here
SCRIPT_NAME = os.path.splitext(os.path.basename(__file__))[0]

from clue.linalg import CC, SparseRowMatrix, SparseVector
from collections import defaultdict
from itertools import product
from math import sqrt
from numpy import  count_nonzero, diag, diagonal, exp
from qiskit.circuit import QuantumCircuit, Parameter
from random import sample, random

## Imports from the local folder
from misc import *

class UndirectedGraph(defaultdict, Experiment):
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
            raise ValueError(f"The given set {V_0=} contains elements not in the graph")
        
        ## Checking V_1
        if V_1 != None:
            if V_0.union(V_1) != set(self.__vertices):
                raise ValueError(f"The given cut does not cover the whole graph")
        else:
            V_1 = set([v for v in self if v not in V_0])

        ## Counting edges from V_0 to V_1
        return sum(len(self[v].intersection(V_1)) for v in V_0)

    def cut_matrix(self) -> SparseRowMatrix:
        r'''Return a 2^n square matrix with the cut counting'''
        M = SparseRowMatrix(2**len(self), field=CC)
        for i,values in enumerate(product(range(2), repeat=len(self))):
            M.increment(i,i,self.cut_cost(set(v for (j,v) in enumerate(self.__vertices) if values[j] == 0)))

        return M

    def cut_values(self) -> set[int]:
        M = self.cut_matrix()
        return set(M[(i,i)] for i in range(M.nrows))
    
    def cut_split(self) -> tuple[SparseRowMatrix, SparseRowMatrix]:
        r'''Return a direct lumping and reduced model from the graph'''
        eigenvalues = dict()
        for i,values in enumerate(product(range(2), repeat=len(self))):
            cost = self.cut_cost(set(v for (j,v) in enumerate(self.__vertices) if values[j] == 0))
            if not cost in eigenvalues:
                eigenvalues[cost] = list()
            eigenvalues[cost].append(i)

        d = len(eigenvalues)
        L = SparseRowMatrix((d,2**(len(self))), CC)
        U = SparseRowMatrix((d,d), CC)
        for i,eigenvalue in enumerate(eigenvalues.keys()):
            m = 1/sqrt(len(eigenvalues[eigenvalue]))
            for j in eigenvalues[eigenvalue]:
                L.increment(i,j,m)
            U.increment(i,i,exp(-1j*eigenvalue))

        return L,U

    def visualize(self, save=False):
        import networkx as nx #pylint: disable=import-error
        import matplotlib.pyplot as plt
        G = nx.Graph()
        G.add_edges_from(self.edges)
        nx.draw_networkx(G)
        if save:
            plt.savefig(save)
        else:
            plt.show()
        plt.close()

    def quantum_cut(self) -> tuple[QuantumCircuit, Parameter]:
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
    
    def quantum_cutB(self) -> tuple[QuantumCircuit, Parameter]:
        r'''
            From :arxiv:`1411.4028v1` the begin Hamiltonian is the exponential matrix

            .. MATH::

                H_B = e^{B}
            
            where `B = sum_{i=1}^n \sigma_i^x`. Since the gates `\sigma_i^x` commutes for different 
            values of `i`, then the matrix `H_B` is the product (i.e., composition) of the individual
            exponential.

            The matrix `\sigma_i^x` diagonalizes in the `\ket{+},\ket{-}` basis, meaning that:

            .. MATH::

                e^{-it\sigma_{j}^x} = H e^{-it H\sigma_j^x H} H,

            where now `H\sigma_j^x H = [(1,0),(0,-1)]^T = \sigma_j^z`. Finally, we can write the exponential
            of `-it\sigma_j^z` as the following 4-gates: `P_j(-i*t) \rightarrow \sigma_j^x \rightarrow P_j(i*t) \rightarrow \sigma_j^x`.

            Hence, to build `H_B` we need to concatenate 6 gates per q-bit:

            ``c.h(j).p(-t,j).x(j).p(t,j).x(j).h(j)``
        '''
        circuit = QuantumCircuit(len(self))
        t = Parameter("t")
        for i in range(len(self)):
            circuit.h(i)
            circuit.p(-t,i)
            circuit.x(i)
            circuit.p(t,i)
            circuit.x(i)
            circuit.h(i)

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
    
    ## NECESSARY METHOD FOR EXPERIMENT
    def size(self) -> int: return len(self)
    def correct_size(self) -> int: return len(self.cut_values())
    def direct(self) -> tuple[SparseRowMatrix, SparseRowMatrix]: return self.cut_split()
    def matrix(self) -> SparseRowMatrix: return self.cut_matrix()
    def matrix_B(self, red_U: ndarray) -> ndarray: 
        if count_nonzero(red_U - diag(diagonal(red_U))): # not diagonal
            return diag([len(self.edges)] + (red_U.shape[0]-1)*[1/len(self.edges)])
        else:
            ## We assume is because it comes from the direct algorithm
            raise NotImplementedError(f"[full-clue] Base hamiltonian not defined when U_P is diagonal")
    def quantum(self) -> tuple[QuantumCircuit, Parameter]: return self.quantum_cut()
    def quantum_B(self) -> tuple[QuantumCircuit, Parameter]: return self.quantum_cutB()
    def data(self): return [len(self.edges)]

## METHODS TO GENERATE THE EXAMPLES OR DATA HEADER
def generate_example(_: str, size: int) -> UndirectedGraph:
    graph = UndirectedGraph.random(size, density=1/3)    
    while(len(graph.edges) == 0):
        graph = UndirectedGraph.random(size, density=1/3)    
    return graph

def generate_header(csv_writer, ttype):
    if ttype in ("clue", "ddsim", "direct"):
        csv_writer.writerow(["size", "edges", "red. ratio", "time_lumping", "memory (MB)", "graph"])
    elif ttype in ("full_clue", "full_direct"):
        csv_writer.writerow(["size", "edges", "time_lumping", "kappa", "time_iteration", "memory (MB)", "graph"])
    elif ttype == "full_ddsim":
        csv_writer.writerow(["size", "edges", "kappa", "time_iteration", "memory (MB)", "graph"])
    else:
        raise NotImplementedError(f"Type of file {ttype} not recognized")

## METHODS TO GENERATE OBSERVABLES
def generate_observable_clue(graph: UndirectedGraph, *_) -> tuple[SparseVector]:
    return tuple([SparseVector.from_list(2**graph.size()*[1], field=CC)])

def generate_observable_ddsim(graph: UndirectedGraph, *_) -> bool:
    return bool(graph)

if __name__ == "__main__":
    ## Processing the arguments
    ttype, script = get_method(*sys.argv)
    m, M = get_size_bounds(*sys.argv)
    timeout = get_timeout(*sys.argv)
    repeats = get_repeats(*sys.argv)
    name = "maxcut"; obs = None
    
    main_script(SCRIPT_DIR, SCRIPT_NAME, name,
                [generate_header, generate_example, generate_observable_clue, generate_observable_ddsim],
                ttype, script,               
                m, M,                        
                timeout, repeats,
                obs)