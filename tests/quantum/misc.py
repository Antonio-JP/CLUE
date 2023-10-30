r'''
    Some auxiliary methods
'''
from __future__ import annotations

from clue import FODESystem, NumericalSubspace, SparseRowMatrix, SparseVector
from csv import writer
from datetime import datetime
from functools import cached_property, lru_cache, reduce
from itertools import product
from math import ceil,inf,sqrt
from mqt.ddsim import DDSIMProvider, CircuitSimulator
from pygraphviz import AGraph, Edge
from numpy import array, cdouble, eye, matmul, ndarray, sqrt
from numpy.linalg import matrix_power
from qiskit import execute
from qiskit.circuit import Parameter, QuantumCircuit
from sympy.parsing.sympy_parser import parse_expr
from time import process_time
from typing import Any, Callable
import os, re, signal, tracemalloc

## Utilities for script
class Timeout(object):
    def __init__(self, seconds):
        self.seconds = 0 if seconds == None else seconds
        self.old = None
    def __enter__(self):
        self.old = signal.signal(signal.SIGALRM, Timeout.alarm_handler)
        signal.alarm(self.seconds)
        return self
    def __exit__(self, type, value, traceback):
        signal.alarm(0)
        signal.signal(signal.SIGALRM, self.old)

    @staticmethod
    def alarm_handler(sgn, _):
        if(sgn == signal.SIGALRM):
            raise TimeoutError
        else:
            raise RuntimeError

## Utilities for different functionalities
def loop(circuit: QuantumCircuit, iterations: int, state_preparation: bool | QuantumCircuit = True, measure: bool = False):
    r'''
        Creates a quantum circuit as a loop of a circuit with fixed number of iterations.

        If ``state_preparation`` is ``bool``: we prepend H to all bits if ``True`` else nothing.
        If ``state_preparation`` is a circuit: we use it as a state preparation circuit. Must fit ``circuit``
    '''
    circuit = circuit.remove_final_measurements(inplace=False)
    final = QuantumCircuit(*circuit.qregs, *circuit.cregs)

    if isinstance(state_preparation, QuantumCircuit):
        final.append(state_preparation, final.qbit_argument_conversion(range(final.num_qubits)), final.cbit_argument_conversion(range(final.num_clbits)))
    elif isinstance(state_preparation, bool):
        for qreg in final.qregs:
            final.h(qreg)
    elif isinstance(state_preparation, int):
        for qbit in final.qbit_argument_conversion([i for i,b in enumerate(bin(state_preparation)[::-1][:-2]) if b == "1"]):
            final.x(qbit) # we set a 1 in each qubit that is needed to set the appropriate input
    
    final.append(circuit.power(iterations), final.qbit_argument_conversion(range(final.num_qubits)), final.cbit_argument_conversion(range(final.num_clbits)))

    if measure:
        final.measure_all()

    return final

def trotter(circuit: QuantumCircuit, gates: list[tuple[QuantumCircuit, list[int]]], order=2):
    r'''
        Computes the Trotter approximation of a list of clauses

        If an operation can be expressed as `e^{t\sum H_i}`, then we can approximate this 
        operation by the product of the matrices:
        
        .. MATH::

            e^{t\sum H_i} = \prod e^{tH_i} + \mathcal{O}(t^2)

        This is called Suzuki-Trotter decomposition. It has several formulas for 
        different orders of approximation. 

        INPUT:

        * ``circuit``: the circuit on which we will apply the Trotter decomposition,
        * ``gates``: the set of operations `e^{tH_i}` in form `(gate, bits)`. Hence, multiplying 
          by `e^{tH_i}` can e done with ``circuit.append(*clauses[i])``.
        * ``order``: order of the Trotter decomposition.

        WARNING:

        The circuit is modified in-place.
    '''
    if order == 2:
        for clause in gates:
            circuit.append(*clause)
    else:
        raise NotImplementedError(f"Trotter decomposition for order {order} not yet implemented")
    
    return circuit

__CACHE_DDSIM_GRAPH = list()
def ddsim_graph(circuit: QuantumCircuit, iterations: int, state_preparation: bool | QuantumCircuit = True, store: str|None = None) -> AGraph:
    r'''
        Method that generates a DiGraph for the execution of a circuit ``iterations`` times.

        This method uses DDSIM to compute the iteration of a circuit a given amount of times and then proceed
        to extract the Directed Graph that represents the final state. With this, we can then proceed to 
        apply linear algebra operations.
    '''
    ## We look in the CACHE
    where = None
    key = (iterations, state_preparation if (not isinstance(state_preparation, QuantumCircuit)) else None)
    for (circ, graphs) in __CACHE_DDSIM_GRAPH:
        if circ is circuit:
            if key in graphs:
                return graphs[key]
            where = graphs
            break
    else:
        where = dict()
        __CACHE_DDSIM_GRAPH.append((circuit, where))

    ## Creating the iterated circuit
    circuit = loop(circuit, iterations, state_preparation, False)
    simulator = CircuitSimulator(circuit)
    simulator.simulate(shots=1) # we run the circuit

    ## We create the Graph structure in networkx
    result = AGraph(string=simulator.export_dd_to_graphviz_str())
    if store != None:
        result.draw(f"./graphs/{store}_{iterations}.png", prog="dot")
    where[key] = result

    return result

def CLUE_circuit(circuit: QuantumCircuit, state_preparation: bool | QuantumCircuit = False, epsilon : float = 1e-6, bound: None|int = None, store: str|None = None) -> tuple[GraphVector]:
    r'''
        Perform approximate CLUE for a quantum circuit using the Orthognal projection approach
    '''
    B: list[GraphVector] = []; go_on = True
    while go_on and (True if bound is None else len(B) < bound):
        print(datetime.now(), f"[CLUE-circuit] Starting case with {len(B)} iterations")
        nextG = GraphVector(ddsim_graph(circuit, len(B), state_preparation, store))
        print(datetime.now(), f"[CLUE-circuit] \t- Computed vector with {len(B)} iterations")
        ## We compute the projection over the already computed basis
        candidate = [nextG]
        for (i,V) in enumerate(B):
            print(datetime.now(), f"[CLUE-circuit] \t- Computing scalar product ({i+1}/{len(B)})...")
            candidate.append((V, -(nextG*V)))

        print(datetime.now(), f"[CLUE-circuit] \t- Building the new candidate...")
        candidate = GraphVector(*candidate)

        print(datetime.now(), f"[CLUE-circuit] \t- Computing its norm...")
        print(datetime.now(), f"[CLUE-circuit] \t- Candidate with norm {float(abs(candidate.norm))}")

        if float(abs(candidate.norm)) < epsilon:
            go_on = False
        else:
            B.append(candidate.normalized_vector)
    print(datetime.now(), f"[CLUE-circuit] Finished CLUE. Dimension of lumping: {len(B)}")
    return tuple(B)

def CLUE_reduced(circuit: QuantumCircuit, lumping : list[GraphVector]) -> ndarray:
    graphs = [ddsim_graph(circuit, i) for i in range(len(lumping))]
    C = array([([cdouble(0)] + [lumping[i].coeff(graphs[j]) for j in range(len(lumping)-1)]) for i in range(len(lumping))])
    last = GraphVector(ddsim_graph(circuit, len(lumping)+1))
    last_in_basis = array([last*v for v in lumping])
    U = array([row+last_in_basis for row in C])

    return U.transpose()

class GraphVector:
    r'''
        Class that represent a linear combination of decision diagrams.

        Input:

        - ``graphs``: a list defining the linear combination. It is a list of `(G, c)` where `c` is a complex
          coefficient and `G` is the graph that we are combining. We allow to only provide the graph `G` when 
          the coefficient is 1.

        This class is immutable, i.e., it does not have any operation that changes the structure inplace.

        This class provides enough method to perform linear algebra operations. More precisely, we provide an 
        implementation for all the operations required for a Gram-Schmitd scheme, namely:

        * We can compute the addition of two :class:`GraphVector`.
        * We can scale a :class: `GraphVector` using a complex scalar.
        * We can compute the scalar product (see :func:`dot`) of two :class:`GraphVector`.
    '''
    def __init__(self, *graphs : tuple[cdouble, AGraph] | AGraph | GraphVector):
        self.__data = dict()
        for graph in graphs:
            if isinstance(graph, (tuple, list)):
                if len(graph) != 2: raise ValueError(f"[GraphVector] Bad format of input: we require that tuples/lists have length 2. Found {graph}")
                G,c = graph
            else:
                c = cdouble(1)
                G = graph

            # Case when
            if isinstance(G, GraphVector):
                for (g,d) in G.to_list():
                    if not g in self.__data:
                        self.__data[g] = cdouble(0)
                    self.__data[g] += c*d
            elif isinstance(G, AGraph):
                if not G in self.__data:
                    self.__data[G] = cdouble(0)
                self.__data[G] += c
            else:
                raise TypeError(f"[GraphVector] Wrong type for graph: either class :AGraph: or class :GraphVector:. Got {G.__class__}")

        ## We clean the zeros
        self.__data = {G: c for (G,c) in self.__data.items() if c != 0}

    def to_list(self) -> list[tuple[AGraph, cdouble]]:
        return list(self.__data.items())
    
    def coeff(self, graph: AGraph):
        return self.__data.get(graph, cdouble(0))

    def __len__(self):
        return len(self.__data)
    
    ##############################################################################################
    ### Vector methods
    ##############################################################################################
    def scalar(self, scale: cdouble) -> GraphVector:
        return GraphVector((self, scale))
    
    def scalar_dot(self, other: GraphVector) -> cdouble:
        ## Scalar product is linear, so we can do a loop until the scalar product of two AGraph
        return sum((c1*c2.conjugate()*GraphVector.graph_scalar(G1,G2) for ((G1,c1),(G2,c2)) in product(self.to_list(), other.to_list())), start=cdouble(0))
    

    __GRAPH_SCALAR_CACHE : dict[set[GraphVector], tuple[tuple[GraphVector,GraphVector],cdouble]] = dict()
    __GRAPH_SCALAR_CACHE_hits : int = 0
    __GRAPH_SCALAR_CACHE_calls : int = 0
    @staticmethod
    # @lru_cache(maxsize=256)
    def graph_scalar(G1: AGraph, G2: AGraph) -> cdouble:
        GraphVector.__GRAPH_SCALAR_CACHE_calls += 1
        key = frozenset((G1, G2))
        if not key in GraphVector.__GRAPH_SCALAR_CACHE:
            print(datetime.now(), f"[graph-scalar] Computing the scalar product of two graphs: ({len(G1.nodes())}, {len(G1.edges())}) * ({len(G2.nodes())}, {len(G2.edges())})")
            if GraphVector.graph_depth(G1) != GraphVector.graph_depth(G2):
                raise TypeError(f"Different depth of two graphs. Can not compute scalar product")
            ## The step with "root" is done first and we start recursion with the first non-root vertex
            edge1 = G1.out_edges("root")[0]; edge2 = G2.out_edges("root")[0]
            val1 = GraphVector.t2c(edge1); val2 = GraphVector.t2c(edge2, True)
            result = val1*val2*GraphVector._graph_scalar(G1, edge1[1], G2, edge2[1])

            GraphVector.__GRAPH_SCALAR_CACHE[key] = ((G1,G2), result)
        else:
            GraphVector.__GRAPH_SCALAR_CACHE_hits += 1

        (g1,g2), result = GraphVector.__GRAPH_SCALAR_CACHE[key]
        if (G1,G2) == (g1,g2): return result
        else: return result.conjugate()
    
    @staticmethod
    def t2c(edge, conjugate=False) -> cdouble:
        to_parse = edge.attr["tooltip"]
        transformations = [
            ("π", "pi"),
            ("^ℯ", "E**"),
            ("ℯ", "*E**"),
            ("ipi ", "I*pi*"),
            ("ipi", "I*pi"),
            ("√(\d+|pi)", lambda M: str(sqrt(float(parse_expr(M.groups()[0]).n()))))
        ]
        output = cdouble(parse_expr(reduce(lambda string, trans : re.sub(trans[0], trans[1], string), [to_parse] + transformations)))
        if conjugate: output = output.conjugate()
        return output

    @staticmethod
    def _graph_scalar(G1: AGraph, root1: Edge, G2: AGraph, root2: Edge):
        ## We first check if we are in the base case
        if root1 == G1.get_node("t") and root2 == G2.get_node("t"):
            return cdouble(1)
        elif root1 == G1.get_node("t") or root2 == G2.get_node("t"):
            raise ValueError("Found the end of a graph but not the other: same depth?")
        
        ## If not, we proceed to an recursion step
        ## Getting left/right for node 1
        edges1 = G1.out_edges(root1)
        left1 = None; right1 = None
        for edge in edges1:
            if "tailport" in edge.attr:
                if "0" in edge.attr["tailport"]:
                    if left1 != None: raise TypeError("Repeated left edge for a vertex?")
                    left1 = (edge[1], GraphVector.t2c(edge))
                elif "1" in edge.attr["tailport"]:
                    if right1 != None: raise TypeError("Repeated right edge for a vertex?")
                    right1 = (edge[1], GraphVector.t2c(edge))
            else:
                raise TypeError("No 'tailport' attribute outside the root?")
        ## Getting left/right for node 1
        edges2 = G2.out_edges(root2)
        left2 = None; right2 = None
        for edge in edges2:
            if "tailport" in edge.attr:
                if "0" in edge.attr["tailport"]:
                    if left2 != None: raise TypeError("Repeated left edge for a vertex?")
                    left2 = (edge[1], GraphVector.t2c(edge,True))
                elif "1" in edge.attr["tailport"]:
                    if right2 != None: raise TypeError("Repeated right edge for a vertex?")
                    right2 = (edge[1], GraphVector.t2c(edge,True))
            else:
                raise TypeError("No 'tailport' attribute outside the root?")
        
        ## Computing the resursion
        left = cdouble(0) if (left1 is None or left2 is None) else left1[1]*left2[1]*GraphVector._graph_scalar(G1, left1[0], G2, left2[0])
        right = cdouble(0) if (right1 is None or right2 is None) else right1[1]*right2[1]*GraphVector._graph_scalar(G1, right1[0], G2, right2[0])
        
        return left + right
        
    @staticmethod
    def graph_depth(G: AGraph) -> int:
        r'''Computes the depth of a Graph. We use a Depth-Search of node "t" from node "root"'''
        goal = G.get_node("t")
        queue = [(G.get_node("root"), 0)]
        while len(queue) > 0:
            C, D = queue.pop()
            if C == goal:
                return D
            for N in G.itersucc(C):
                queue.append((N, D+1))
    
    @cached_property
    def norm(self) -> cdouble:
        return sqrt(self * self)
    
    @cached_property
    def normalized_vector(self):
        return self * (1/self.norm)

    ##############################################################################################
    ### Arithmetic methods
    ##############################################################################################
    def __add__(self, other: GraphVector) -> GraphVector:
        if not isinstance(other, GraphVector): 
            if other == 0:
                return self
            return NotImplemented
        return GraphVector(self, other)
    
    def __radd__(self, other: GraphVector) -> GraphVector:
        return self.__add__(other)
    
    def __neg__(self) -> GraphVector:
        return GraphVector(*[(G,-c) for (G,c) in self.to_list()])
    
    def __sub__(self, other: GraphVector) -> GraphVector:
        return self + (-other)
    
    def __rsub__(self, other: GraphVector) -> GraphVector:
        return (-self).__add__(other)
    
    def __mul__(self, other: cdouble | GraphVector) -> GraphVector | cdouble:
        if isinstance(other, GraphVector):
            return self.scalar_dot(other)
        else:
            try:
                other = cdouble(other)
            except:
                return NotImplemented
            return self.scalar(other)
    
    def __rmul__(self, other: cdouble | GraphVector) -> GraphVector | cdouble:
        if isinstance(other, GraphVector):
            return other.scalar_dot(self)
        return self.__mul__(other)

## General Experiment interface
class Experiment:
    r'''Interface for experiments to run other methods of this module'''
    def size(self) -> int: raise NotImplementedError(f"Method for getting 'qbits size' not implemented")
    def correct_size(self) -> int: raise NotImplementedError(f"Method for getting 'correct_size' not implemented")
    def direct(self) -> tuple[SparseRowMatrix, SparseRowMatrix]: raise NotImplementedError(f"Method for getting 'direct lumping' not implemented")
    def matrix(self) -> SparseRowMatrix: raise NotImplementedError(f"Method for getting 'matrix' not implemented")
    def matrix_B(self, red_U: ndarray) -> ndarray: raise NotImplementedError(f"Method for getting 'matrix begin' not implemented")
    def quantum(self) -> tuple[QuantumCircuit, Parameter]: raise NotImplementedError(f"Method for getting 'quantum circuit' not implemented")
    def quantum_B(self) -> tuple[QuantumCircuit, Parameter]: raise NotImplementedError(f"Method for getting 'quantum begin' not implemented")
    def data(self): pass

## Different execution methods
def clue_reduction(name: str, 
                   generate_example: Callable[[Any],Experiment], 
                   generate_observable: Callable[[Experiment, Any], tuple[SparseVector]], 
                   generate_data: Callable[[Experiment,Any], tuple], 
                   result_file, *args, timeout:float=0, **kwds): 
    r'''
        This method computes the CLUE lumping.

        This method compute the CLUE lumping of an example. The arguments are as follows:

        * ``generate_example``: from [name, *args, **kwds] generates a valid example.
        * ``generate_observable``: from [system, *args, **kwds] generates a valid observable.
        * ``generate_data``: from [size, time, ratio, memory, experiment] generates the output row for CSV
        * ``result_fle``: CSV writer to put the result
        * ``args``: arguments for the generate functions.
        * ``timeout``: timeout used during the lumping.
        * ``kwds``: named arguments for the generate functions.

        This method generates aninstance of a problem, create a valid observable for the lumping and compute the lumping w.r.t.
        that observable. Then it stores the reduction ratio, time used, and memory spent in the execution.

        It stores the execution time of the lumping plus the reduction ratio. It stores the result on ``result_file``. It uses
        method ``generate_data`` to format the CSV output.
    '''
    print(datetime.now(), f"%%% [clue @ {name}] Computing CLUE reduction for {name} and arguments {args} and {kwds}...", flush=True)
    experiment = generate_example(name, *args, **kwds)
    try: 
        true_size = experiment.correct_size()
    except:
        true_size = None

    print(datetime.now(), f"%%% [clue @ {name}] Creating the full system to apply CLUE", flush=True)
    system = FODESystem.LinearSystem(experiment.matrix(), lumping_subspace=NumericalSubspace)
    obs = generate_observable(experiment, *args, **kwds)
    
    print(datetime.now(), f"%%% [clue @ {name}] Computing the lumped system...", flush=True)
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            ctime = process_time()
            lumped = system.lumping(obs, print_reduction=False, print_system=False)
            ctime = process_time()-ctime
    except TimeoutError:
        print(datetime.now(), f"%%% [clue @ {name}] Timeout reached for execution", flush=True)
        ctime = inf
    memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(datetime.now(), f"%%% [clue @ {name}] Checking correct size (if possible)", flush=True)
    if true_size != None and ctime < inf:
        if true_size != lumped.size:
            print(datetime.now(), f"%%% [clue @ {name}] ERROR!! Found weird dimension in lumping -- \n%%% \t* Expected: {true_size}\n%%% \t* Got: {lumped.size}\n%%% \t* Experiment: {experiment}", flush=True)
    result_file.writerow(generate_data(experiment, lumped.size/system.size, ctime, memory))

    return ctime

def ddsim_reduction(name: str, 
                   generate_example: Callable[[Any],Experiment], 
                   generate_observable: Callable[[Experiment, Any], bool|QuantumCircuit], 
                   generate_data: Callable[[Experiment,Any], tuple], 
                   result_file, *args, timeout:float=0, **kwds) -> float: 
    r'''
        This method computes the DDSIM lumping

        This method compute the DDSIM iteration required for lumping of an example. The arguments are as follows:

        * ``generate_example``: from [name, *args, **kwds] generates a valid example.
        * ``generate_observable``: from [system, *args, **kwds] generates whether the input is H or not.
        * ``generate_data``: from [size, time, ratio, memory, experiment] generates the output row for CSV
        * ``result_fle``: CSV writer to put the result
        * ``args``: arguments for the generate functions.
        * ``timeout``: timeout used during the lumping.
        * ``kwds``: named arguments for the generate functions.
         
        It stores the execution time of the lumping plus the reduction ratio. It stores the result on ``result_file``. It uses
        method ``generate_data`` to format the CSV output.
    '''
    print(datetime.now(), f"%%% [ddsim @ {name}] Computing DDSIM reduction for {name} and arguments {args} and {kwds}...", flush=True)
    experiment = generate_example(name, *args, **kwds)
    try: 
        true_size = experiment.correct_size()
    except:
        true_size = 2**experiment.size()

    print(datetime.now(), f"%%% [ddsim @ {name}] Creating the full circuit and job to simulate with DDSIM", flush=True)
    circuit, par = experiment.quantum()
    if par != None: circuit = circuit.bind_parameters({par: 1/(1000*true_size)})
    
    print(datetime.now(), f"%%% [ddsim @ {name}] Computing the reductio using DDSIM...", flush=True)
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            ctime = process_time()
            ## Executing the circuit one time
            L = CLUE_circuit(circuit, generate_observable(experiment, *args, **kwds), store = f"{name}_{args}")
            _ = CLUE_reduced(circuit, L)
            ctime = process_time()-ctime
    except TimeoutError:
        print(datetime.now(), f"%%% [ddsim @ {name}] Timeout reached for execution", flush=True)
        ctime = inf
    memory = tracemalloc.get_traced_memory()[1] / (2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(datetime.now(), f"%%% [ddsim @ {name}] Storing the data...", flush=True)
    result_file.writerow(generate_data(experiment, len(L)/2**experiment.size(), ctime, memory))

    return ctime

def direct_reduction(name: str, 
                   generate_example: Callable[[Any],Experiment], 
                   generate_observable: Callable[[Experiment, Any], tuple[SparseVector]], 
                   generate_data: Callable[[Experiment,Any], tuple], 
                   result_file, *args, timeout:float=0, **kwds) -> float: 
    r'''
        This method computes the CLUE lumping in a direct fashion

        This method compute the CLUE lumping of an example. The arguments are as follows:

        * ``generate_example``: from [name, *args, **kwds] generates a valid example.
        * ``generate_observable``: from [system, *args, **kwds] generates a valid observable.
        * ``generate_data``: from [size, time, ratio, memory, experiment] generates the output row for CSV
        * ``result_fle``: CSV writer to put the result
        * ``args``: arguments for the generate functions.
        * ``timeout``: timeout used during the lumping.
        * ``kwds``: named arguments for the generate functions.

        This method generates an instance of a problem, create a valid observable for the lumping and compute the lumping w.r.t.
        that observable. Then it stores the reduction ratio, time used, and memory spent in the execution.

        It stores the execution time of the lumping plus the reduction ratio. It stores the result on ``result_file``. It uses
        method ``generate_data`` to format the CSV output.
    '''
    print(datetime.now(), f"%%% [direct @ {name}] Computing Direct CLUE reduction for {name} and arguments {args} and {kwds}...", flush=True)
    experiment = generate_example(name, *args, **kwds)
    
    print(datetime.now(), f"%%% [direct @ {name}] Computing the lumped system...", flush=True)
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            ctime = process_time()
            L, _ = experiment.direct()
            ctime = process_time()-ctime
    except TimeoutError:
        print(datetime.now(), f"%%% [direct @ {name}] Timeout reached for execution", flush=True)
        ctime = inf
    memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(datetime.now(), f"%%% [direct @ {name}] Storing the data...", flush=True)
    result_file.writerow(generate_data(experiment, L.nrows/L.ncols, ctime, memory))

    return ctime

def clue_iteration(name: str, 
                   generate_example: Callable[[Any],Experiment], 
                   generate_observable: Callable[[Experiment, Any], tuple[SparseVector]], 
                   generate_data: Callable[[Experiment,Any], tuple],
                   result_file, iterations: int, *args, timeout:float=0, **kwds) -> float: 
    r'''
        This method computes the CLUE iteration.

        This method compute the CLUE lumping of an example. The arguments are as follows:

        * ``generate_example``: from [name, *args, **kwds] generates a valid example.
        * ``generate_observable``: from [system, *args, **kwds] generates a valid observable.
        * ``generate_data``: from [size, time, iters, time, memory, experiment] generates the output row for CSV
        * ``iterations``: number of iterations to compute.
        * ``result_fle``: CSV writer to put the result
        * ``args``: arguments for the generate functions.
        * ``timeout``: timeout used during the lumping.
        * ``kwds``: named arguments for the generate functions.

        This method generates an instance of a problem, create a valid observable for the lumping and compute the lumping w.r.t.
        that observable.

        It stores the result on ``result_file``. It uses method ``generate_data`` to format the CSV output.
    '''
    print(datetime.now(), f"%%% [full-clue @ {name}] Computing CLUE iterations ({iterations}) for {name} and arguments {args} and {kwds}...", flush=True)
    experiment = generate_example(name, *args, **kwds)
    
    try: 
        true_size = experiment.correct_size()
    except:
        true_size = None

    print(datetime.now(), f"%%% [full-clue @ {name}] Creating the full system to apply CLUE", flush=True)
    system = FODESystem.LinearSystem(experiment.matrix(), lumping_subspace=NumericalSubspace)
    obs = generate_observable(experiment, *args, **kwds)
    
    print(datetime.now(), f"%%% [full-clue @ {name}] Computing the lumped system...", flush=True)
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            lump_time = process_time()
            ## Executing the circuit one time
            lumped = system.lumping(obs, print_reduction=False, print_system=False)
            lump_time = process_time()-lump_time

            print(datetime.now(), f"%%% [full-clue @ {name}] Checking correct size (if possible)", flush=True)
            if true_size != None:
                if true_size != lumped.size:
                    print(datetime.now(), f"%%% [clue @ {name}] ERROR!! Found weird dimension in lumping -- \n%%% \t* Expected: {true_size}\n%%% \t* Got: {lumped.size}\n%%% \t* Experiment: {experiment}", flush=True)

            print(datetime.now(), f"%%% [full-clue @ {name}] Getting the reduced U_P", flush=True)
            U_P = lumped.construct_matrices("polynomial")[0].to_numpy(dtype=cdouble)
            print(datetime.now(), f"%%% [full-clue @ {name}] Getting the reduced U_B", flush=True)
            try:
                U_B = experiment.matrix_B(U_P)
            except:
                print(datetime.now(), f"%%% [full-clue @ {name}] No reduced U_B: going to identity", flush=True)
                U_B = eye(U_P.shape[0])
            
            print(datetime.now(), f"%%% [full-clue @ {name}] Computing the iteration (U_P*U_B)^iterations", flush=True)
            U = matmul(U_P, U_B)
            it_time = process_time()
            _ = matrix_power(U, iterations)
            it_time = process_time() - it_time
    except TimeoutError:
        print(datetime.now(), f"%%% [full-clue @ {name}] Timeout reached for execution", flush=True)
        lump_time = inf
        it_time = inf
    finally:
        memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
        tracemalloc.stop()
        result_file.writerow(generate_data(experiment, lump_time, iterations, it_time, memory))
        
    return lump_time + it_time

def ddsim_iteration(name: str, 
                   generate_example: Callable[[Any],Experiment], 
                   generate_observable: Callable[[Experiment, Any], bool|QuantumCircuit], 
                   generate_data: Callable[[Experiment,Any], tuple],
                   result_file, iterations: int, *args, timeout:float=0, **kwds) -> float: 
    r'''
        This method computes the DDSIM iteration.

        This method compute the DDSIM iteration of an example. The arguments are as follows:

        * ``generate_example``: from [name, *args, **kwds] generates a valid example.
        * ``generate_observable``: from [system, *args, **kwds] generates whether the input is H or not.
        * ``generate_data``: from [size, iters, time, memory, experiment] generates the output row for CSV
        * ``iterations``: number of iterations to compute.
        * ``result_fle``: CSV writer to put the result
        * ``args``: arguments for the generate functions.
        * ``timeout``: timeout used during the lumping.
        * ``kwds``: named arguments for the generate functions.

        This method generates an instance of a problem, create the associated quantum circuits and execute it ``iteration`` times.

        It stores the result on ``result_file``. It uses method ``generate_data`` to format the CSV output.
    '''
    print(datetime.now(), f"%%% [full-ddsim @ {name}] Computing DDSIM iterations ({iterations}) for {name} and arguments {args} and {kwds}...", flush=True)
    experiment = generate_example(name, *args, **kwds)

    print(datetime.now(), f"%%% [full-ddsim @ {name}] Creating the full circuit and job to simulate with DDSIM", flush = True)
    U_P, par = experiment.quantum()
    if par != None: U_P = U_P.bind_parameters({par: 1/(2**experiment.size()*10*iterations)})
    try:
        U_B, par = experiment.quantum_B()
        if par != None: U_B = U_B.bind_parameters({par: 1/(2**experiment.size()*10*iterations)})
        U_B.append(U_P, U_P.qregs[0]) # Now U_B is the alternate circuit U_B * U_P
    except NotImplementedError: # U_B do not exist
        U_B = U_P
    circuit = loop(U_B, iterations, generate_observable(experiment, *args, **kwds), True)
    backend = DDSIMProvider().get_backend("qasm_simulator")
    
    print(datetime.now(), f"%%% [full-ddsim @ {name}] Computing the simulation of the circuit...", flush = True)
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            ctime = process_time()
            ## Executing the circuit one time
            job = execute(circuit, backend, shots=1)
            job.result()
            ctime = process_time()-ctime
    except TimeoutError:
        print(datetime.now(), f"%%% [full-ddsim @ {name}] Timeout reached for execution", flush = True)
        ctime = inf
    memory = tracemalloc.get_traced_memory()[1] / (2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(datetime.now(), f"%%% [full-ddsim @ {name}] Storing the data...", flush = True)
    result_file.writerow(generate_data(experiment, iterations, ctime, memory))

    return ctime

def direct_iteration(name: str, 
                   generate_example: Callable[[Any],Experiment], 
                   generate_observable: Callable[[Experiment, Any], tuple[SparseVector]], 
                   generate_data: Callable[[Experiment,Any], tuple],
                   result_file, iterations: int, *args, timeout:float=0, **kwds) -> float: 
    r'''
        This method computes the CLUE iteration.

        This method compute the CLUE lumping of an example. The arguments are as follows:

        * ``generate_example``: from [name, *args, **kwds] generates a valid example.
        * ``generate_observable``: from [system, *args, **kwds] generates a valid observable.
        * ``generate_data``: from [size, time, iters, time, memory, experiment] generates the output row for CSV
        * ``iterations``: number of iterations to compute.
        * ``result_fle``: CSV writer to put the result
        * ``args``: arguments for the generate functions.
        * ``timeout``: timeout used during the lumping.
        * ``kwds``: named arguments for the generate functions.

        This method generates an instance of a problem, create a valid observable for the lumping and compute the lumping w.r.t.
        that observable.

        It stores the result on ``result_file``. It uses method ``generate_data`` to format the CSV output.
    '''
    print(datetime.now(), f"%%% [full-direct @ {name}] Computing Direct CLUE iterations ({iterations}) for {name} and arguments {args} and {kwds}...", flush=True)
    experiment = generate_example(name, *args, **kwds)
    
    print(datetime.now(), f"%%% [full-direct @ {name}] Computing the lumped system...", flush=True)
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            lump_time = process_time()
            _, U = experiment.direct()
            lump_time = process_time()-lump_time

            print(datetime.now(), f"%%% [full-direct @ {name}] Getting the reduced U_P", flush=True)
            U_P = U.to_numpy(dtype=cdouble)
            print(datetime.now(), f"%%% [full-direct @ {name}] Getting the reduced U_B", flush=True)
            try:
                U_B = experiment.matrix_B(U_P)
            except:
                print(datetime.now(), f"%%% [full-direct @ {name}] No reduced U_B: going to identity", flush=True)
                U_B = eye(U_P.shape[0])
            
            print(datetime.now(), f"%%% [full-direct @ {name}] Computing the iteration (U_P*U_B)^iterations", flush=True)
            U = matmul(U_P, U_B)
            it_time = process_time()
            _ = matrix_power(U, iterations)
            it_time = process_time() - it_time
    except TimeoutError:
        print(datetime.now(), f"%%% [full-direct @ {name}] Timeout reached for execution", flush=True)
        lump_time = inf 
        it_time = inf
    finally:
        memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
        tracemalloc.stop()
        result_file.writerow(generate_data(experiment, lump_time, iterations, it_time, memory))
        
    return lump_time + it_time

def generate_data(example: Experiment, *args) -> list:
    return [example.size()] + example.data() + [*args] + [repr(example)]

## Reading-argument methods
def get_size_bounds(*argv) -> tuple[int,int]:
    ## Getting the minimum
    if "-m" in argv:
        ind = argv.index("-m") + 1
        try:
            m = int(argv[ind])
            assert m > 0
        except:
            raise TypeError(f"Script argument [-m] not well used: we required a follow-up positive integer")
    else:
        m = 5
        
    ## Getting the maximum
    if "-M" in argv:
        ind = argv.index("-M") + 1
        try:
            M = int(argv[ind])
            assert M >= m
        except:
            raise TypeError(f"Script argument [-M] not well used: we required a follow-up integer bigger than the minimum (i.e., {m})")
    else:
        M = m
    return (m,M)

def get_method(*argv) -> tuple[str,Callable]:
    if "-t" in argv:
        ind = argv.index("-t") + 1
        if len(argv) > ind:
            if argv[ind] == "clue": return "clue", clue_reduction
            if argv[ind] == "ddsim": return "ddsim", ddsim_reduction
            if argv[ind] == "direct": return "direct", direct_reduction
            if argv[ind] == "full_clue": return "full_clue", clue_iteration
            if argv[ind] == "full_ddsim": return "full_ddsim", ddsim_iteration
            if argv[ind] == "full_direct": return "full_direct", direct_iteration

        raise TypeError("Script argument [-t] not well used: we required a follow-up name in ('clue','ddsim','direct','full_clue','full_ddsim','full_clue')")
    else:
        return "clue", clue_reduction
            
def get_timeout(*argv) -> int:
    if "-to" in argv:
        ind = argv.index("-to") + 1
        try:
            t = int(argv[ind])
            assert t > 0
        except:
            raise TypeError(f"Script argument [-to] not well used: we required a follow-up positive integer")
    else:
        t = None
    return t

def get_repeats(*argv) -> int:
    if "-r" in argv:
        ind = argv.index("-r") + 1
        try:
            r = int(argv[ind])
            assert r > 0
        except:
            raise TypeError(f"Script argument [-r] not well used: we required a follow-up positive integer")
    else:
        r = 1
    return r

def get_rem_timeout(rem_timeout, used_time):
    if rem_timeout != None:
        rem_timeout -= used_time
        if rem_timeout < 0:
            raise TimeoutError
        else:
            return ceil(rem_timeout)
    return None

def main_script(dir: str, filename: str, name: str,         # directory and filename of the script; name of the experiment
                methods: list[Callable] | tuple[Callable],  # list of methods with 5 methods to generate the appropriate script
                ttype: str, script: Callable,               # the type and the script to be run from misc.py
                m: int, M: int,                             # validated bounds for size to be executed
                timeout: int | None, repeats: int,          # validated timeout and number of repetitions
                observables : list[Any] | None,             # list of observables to be used
                **kwds                                      # other arguments for the script
    ):
    ## Processing arguments
    generate_header, generate_example, generate_observable_clue, generate_observable_ddsim = methods
    ## Running the script
    existed = os.path.exists(os.path.join(dir, "results", f"[result]{filename}_{ttype}.csv"))
    with open(os.path.join(dir, "results", f"[result]{filename}_{ttype}.csv"), "at" if existed else "wt") as result_file:
        csv_writer = writer(result_file)
        if not existed:
            generate_header(csv_writer, ttype)
        print(f"##################################################################################")
        print(f"### EXECUTION ON {name.upper()} [{m=}, {M=}, {repeats=}, method={ttype}]")
        print(f"##################################################################################")
        for size in range(m, M+1):
            obs_to_use = ([0] + ["H"] + list(range(1, 2**(generate_example(name, size, 0).quantum()[0].num_qubits)))) if (observables != None and len(observables) == 0) else ["def"] if observables == None else observables
            script_args = [name, generate_example, generate_observable_clue if not ("ddsim" in ttype) else generate_observable_ddsim, generate_data, csv_writer]
            for execution in range(1,repeats+1):
                rem_timeout = timeout
                try:
                    for i,observable in enumerate(obs_to_use):
                        print(datetime.now(), f"### ++ Starting execution {execution}/{repeats} ({size=}, observable={i+1}/{len(obs_to_use)})")
                        if not ("full" in ttype): # No iterations are required
                            used_time = script(*(script_args + [size] + ([observable] if observables != None else [])), timeout=rem_timeout, **kwds)
                            rem_timeout = get_rem_timeout(rem_timeout, used_time)
                            result_file.flush()
                        else:
                            used_time = 0
                            for it in (1,ceil(sqrt(2**size))):#,1000):#,10000)
                                print(datetime.now(), f"    ++++ Case with {it} iterations.")
                                it_used_time = script(*(script_args + [it, size] + ([observable] if observables != None else [])), timeout=rem_timeout, **kwds)
                                rem_timeout = get_rem_timeout(rem_timeout, it_used_time)
                                used_time += it_used_time
                            result_file.flush()
                        print(datetime.now(), f"### -- Finished execution {execution}/{repeats} ({size=}, observable={i+1}/{len(obs_to_use)}): took {used_time} s.")
                        if rem_timeout != None:
                            rem_timeout -= used_time
                            if rem_timeout < 0:
                                break
                            else:
                                rem_timeout = ceil(rem_timeout)
                except TimeoutError:
                    print(datetime.now(), f"### -- Finished execution {execution}/{repeats} ({size=}): reached Timeout.")
