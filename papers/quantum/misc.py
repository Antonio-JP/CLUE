r'''
    Some auxiliary methods
'''
from __future__ import annotations

from clue import FODESystem, NumericalSubspace, SparseRowMatrix, SparseVector
from csv import writer
from math import ceil,inf,sqrt
from mqt import ddsim
from numpy import cdouble, eye, matmul, ndarray
from numpy.linalg import matrix_power
from qiskit import transpile
from qiskit.circuit import Parameter, QuantumCircuit
from time import process_time
from typing import Any, Callable
import os, signal, tracemalloc

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
    elif state_preparation:
        for qreg in final.qregs:
            final.h(qreg)
    
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

## General Experiment interface
class Experiment:
    r'''Interface for experiments to run other methods of this module'''
    def size(self) -> int: raise NotImplementedError(f"Method for getting 'qbits size' not implemented")
    def correct_size(self) -> int: raise NotImplementedError(f"Method for getting 'correct_size' not implemented")
    def direct(self) -> tuple[SparseRowMatrix, SparseRowMatrix]: raise NotImplementedError(f"Method for getting 'direct lumping' not implemented")
    def matrix(self) -> SparseRowMatrix: raise NotImplementedError(f"Method for getting 'matrix' not implemented")
    def matrix_B(self, red_U: ndarray) -> ndarray: raise NotImplementedError(f"Method for getting 'matrix begin' not implemented")
    def build_quantum(self, *, iterations: int = 1, state_preparation: bool | QuantumCircuit = True, measure: bool = False, append_B: bool = False, par_value: float = None) -> QuantumCircuit:
        r'''
            Builds the quantum circuit for the experiment.

            INPUT:

            * ``iterations``: number of iterations to apply the circuit.
            * ``state_preparation``: if ``True``, applies H to all qubits before the circuit.
              If a circuit, it is used as a state preparation circuit.
            * ``measure``: if ``True``, measures all qubits at the end of the circuit.
            * ``append_B``: if ``True``, appends the begin Hamiltonian at the end of the circuit.
        '''
        circuit = QuantumCircuit(self.size())
        par = Parameter("t")
        if state_preparation is True:
            circuit.h(range(circuit.num_qubits))
        elif isinstance(state_preparation, QuantumCircuit):
            circuit.append(state_preparation, range(circuit.num_qubits))

        for _ in range(iterations):
            circuit, par = self.quantum(circuit, par)
            if append_B:
                try:
                    circuit, par = self.quantum_B(circuit, par)
                except NotImplementedError:
                    print(f"WARNING: quantum_B not implemented for {self.__class__.__name__}, skipping append_B")
                    pass

        if measure:
            circuit.measure_all()

        if par != None: circuit = circuit.assign_parameters({par: 1/(2**self.size()*10*iterations)})      

        return circuit
    def quantum(self, circuit: QuantumCircuit, dt: float | Parameter) -> tuple[QuantumCircuit, Parameter]: raise NotImplementedError(f"Method for getting 'quantum circuit' not implemented")
    def quantum_B(self, circuit: QuantumCircuit, dt: float | Parameter) -> tuple[QuantumCircuit, Parameter]: raise NotImplementedError(f"Method for getting 'quantum begin' not implemented")
    def data(self): pass

    @staticmethod
    def generate_example(name: str, size: int) -> Experiment:
        raise NotImplementedError

    @staticmethod
    def generate_header(csv_writer, ttype):
        raise NotImplementedError

    ## METHODS TO GENERATE OBSERVABLES
    @staticmethod
    def generate_observable_clue(graph: Experiment, size: int, obs: int) -> tuple[SparseVector]:
        raise NotImplementedError

    @staticmethod
    def generate_observable_ddsim(graph: Experiment, size: int, obs: int) -> bool:
        raise NotImplementedError

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
    print(f"%%% [clue @ {name}] Computing CLUE reduction for {name} and arguments {args} and {kwds}...", flush=True)
    experiment = generate_example(name, *args, **kwds)
    try: 
        true_size = experiment.correct_size()
    except:
        true_size = None

    print(f"%%% [clue @ {name}] Creating the full system to apply CLUE", flush=True)
    system = FODESystem.LinearSystem(experiment.matrix(), lumping_subspace=NumericalSubspace)
    obs = generate_observable(experiment, *args, **kwds)
    
    print(f"%%% [clue @ {name}] Computing the lumped system...", flush=True)
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            ctime = process_time()
            lumped = system.lumping(obs, print_reduction=False, print_system=False)
            ctime = process_time()-ctime
    except TimeoutError:
        print(f"%%% [clue @ {name}] Timeout reached for execution", flush=True)
        ctime = inf
    memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(f"%%% [clue @ {name}] Checking correct size (if possible)", flush=True)
    if true_size != None and ctime < inf:
        if true_size != lumped.size:
            print(f"%%% [clue @ {name}] ERROR!! Found weird dimension in lumping -- \n%%% \t* Expected: {true_size}\n%%% \t* Got: {lumped.size}\n%%% \t* Experiment: {experiment}", flush=True)
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
    print(f"%%% [ddsim @ {name}] Computing DDSIM reduction for {name} and arguments {args} and {kwds}...", flush=True)
    experiment = generate_example(name, *args, **kwds)
    try: 
        true_size = experiment.correct_size()
    except:
        true_size = 2**experiment.size()

    print(f"%%% [ddsim @ {name}] Creating the full circuit and job to simulate with DDSIM", flush=True)
    circuit = experiment.build_quantum(
        iterations=true_size, 
        state_preparation=generate_observable(experiment, *args, **kwds),
        measure=True,
        par_value=1/(1000*true_size))
    
    backend = ddsim.DDSIMProvider().get_backend("qasm_simulator")
    
    print(f"%%% [ddsim @ {name}] Computing the simulation of the circuit...", flush=True)
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            ctime = process_time()
            ## Executing the circuit one time
            circuit_tr = transpile(circuit, backend)
            job = backend.run(circuit_tr, shots=1)
            job.result()
            ctime = process_time()-ctime
    except TimeoutError:
        print(f"%%% [ddsim @ {name}] Timeout reached for execution", flush=True)
        ctime = inf
    memory = tracemalloc.get_traced_memory()[1] / (2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(f"%%% [ddsim @ {name}] Storing the data...", flush=True)
    result_file.writerow(generate_data(experiment, true_size/2**experiment.size(), ctime, memory))

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
    print(f"%%% [direct @ {name}] Computing Direct CLUE reduction for {name} and arguments {args} and {kwds}...", flush=True)
    experiment = generate_example(name, *args, **kwds)
    
    print(f"%%% [direct @ {name}] Computing the lumped system...", flush=True)
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            ctime = process_time()
            L, _ = experiment.direct()
            ctime = process_time()-ctime
    except TimeoutError:
        print(f"%%% [direct @ {name}] Timeout reached for execution", flush=True)
        ctime = inf
    memory = tracemalloc.get_traced_memory()[1]/(2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(f"%%% [direct @ {name}] Storing the data...", flush=True)
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
    print(f"%%% [full-clue @ {name}] Computing CLUE iterations ({iterations}) for {name} and arguments {args} and {kwds}...", flush=True)
    experiment = generate_example(name, *args, **kwds)
    
    try: 
        true_size = experiment.correct_size()
    except:
        true_size = None

    print(f"%%% [full-clue @ {name}] Creating the full system to apply CLUE", flush=True)
    system = FODESystem.LinearSystem(experiment.matrix(), lumping_subspace=NumericalSubspace)
    obs = generate_observable(experiment, *args, **kwds)
    
    print(f"%%% [full-clue @ {name}] Computing the lumped system...", flush=True)
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            lump_time = process_time()
            ## Executing the circuit one time
            lumped = system.lumping(obs, print_reduction=False, print_system=False)
            lump_time = process_time()-lump_time

            print(f"%%% [full-clue @ {name}] Checking correct size (if possible)", flush=True)
            if true_size != None:
                if true_size != lumped.size:
                    print(f"%%% [clue @ {name}] ERROR!! Found weird dimension in lumping -- \n%%% \t* Expected: {true_size}\n%%% \t* Got: {lumped.size}\n%%% \t* Experiment: {experiment}", flush=True)

            print(f"%%% [full-clue @ {name}] Getting the reduced U_P", flush=True)
            U_P = lumped.construct_matrices("polynomial")[0].to_numpy(dtype=cdouble)
            print(f"%%% [full-clue @ {name}] Getting the reduced U_B", flush=True)
            try:
                U_B = experiment.matrix_B(U_P)
            except:
                print(f"%%% [full-clue @ {name}] No reduced U_B: going to identity", flush=True)
                U_B = eye(U_P.shape[0])
            
            print(f"%%% [full-clue @ {name}] Computing the iteration (U_P*U_B)^iterations", flush=True)
            U = matmul(U_P, U_B)
            it_time = process_time()
            _ = matrix_power(U, iterations)
            it_time = process_time() - it_time
    except TimeoutError:
        print(f"%%% [full-clue @ {name}] Timeout reached for execution", flush=True)
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
    print(f"%%% [full-ddsim @ {name}] Computing DDSIM iterations ({iterations}) for {name} and arguments {args} and {kwds}...", flush=True)
    experiment = generate_example(name, *args, **kwds)

    print(f"%%% [full-ddsim @ {name}] Creating the full circuit and job to simulate with DDSIM", flush = True)
    circuit = experiment.build_quantum(
        iterations=iterations,
        state_preparation=generate_observable(experiment, *args, **kwds),
        measure=True,
        par_value=1/(2**experiment.size()*10*iterations)
    )

    backend = ddsim.DDSIMProvider().get_backend("qasm_simulator")
    
    print(f"%%% [full-ddsim @ {name}] Computing the simulation of the circuit...", flush = True)
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            ctime = process_time()
            ## Executing the circuit one time
            circuit_tr = transpile(circuit, backend)
            job = backend.run(circuit_tr, shots=1)
            job.result()
            ctime = process_time()-ctime
    except TimeoutError:
        print(f"%%% [full-ddsim @ {name}] Timeout reached for execution", flush = True)
        ctime = inf
    memory = tracemalloc.get_traced_memory()[1] / (2**20) # maximum memory usage in MB
    tracemalloc.stop()

    print(f"%%% [full-ddsim @ {name}] Storing the data...", flush = True)
    result_file.writerow(generate_data(experiment, iterations, ctime, memory))

    return ctime

def quokka_iteration(name: str, 
                   generate_example: Callable[[Any],Experiment], 
                   generate_observable: Callable[[Experiment, Any], bool|QuantumCircuit], 
                   generate_data: Callable[[Experiment,Any], tuple],
                   result_file, iterations: int, *args, timeout:float=0, **kwds) -> float: 
    r'''
        This method computes the Quokka# iteration.

        This method compute the Quokka# iteration of an example. The arguments are as follows:

        * ``generate_example``: from [name, *args, **kwds] generates a valid example.
        * ``generate_observable``: from [system, *args, **kwds] generates whether the input is H or not.
        * ``generate_data``: from [size, iters, time, memory, experiment] generates the output row for CSV
        * ``iterations``: number of iterations to compute.
        * ``result_fle``: CSV writer to put the result
        * ``args``: arguments for the generate functions.
        * ``timeout``: timeout used during the lumping.
        * ``kwds``: named arguments for the generate functions.

        This method generates an instance of a problem, create the associated quantum circuits and execute it ``iteration`` times using Quokka#

        It stores the result on ``result_file``. It uses method ``generate_data`` to format the CSV output.
    '''
    import quokka_sharp as qk
    import tempfile
    from qiskit import qasm2

    print(f"%%% [quokka# @ {name}] Computing DDSIM iterations ({iterations}) for {name} and arguments {args} and {kwds}...", flush=True)
    experiment = generate_example(name, *args, **kwds)
    print(f"%%% [quokka# @ {name}] {experiment.size()} - {experiment.data()}\n\t{experiment}", flush = True)

    print(f"%%% [quokka# @ {name}] Creating the full circuit and job to simulate with DDSIM", flush = True)
    circuit = experiment.build_quantum(
        iterations=iterations,
        state_preparation=generate_observable(experiment, *args, **kwds),
        measure=False,
        par_value=1/(2**experiment.size()*10*iterations)
    )

    with tempfile.NamedTemporaryFile(suffix=".qasm", delete=False) as f:
        file_name = f.name
        f.close()  # Close the file so that it can be used by qasm2
        print(f"%%% [quokka# @ {name}] Writing the circuit to a temporary QASM file...", flush = True)
        qasm2.dump(circuit, file_name)
        
        print(f"%%% [quokka# @ {name}] Computing the simulation of the circuit...", flush = True)
        tracemalloc.start()
        try:
            with(Timeout(timeout)):
                ## Encoding into CNF
                enc_time = process_time()
                print(f"%%% [quokka# @ {name}] Reading QASM file...", flush = True)
                circuit_cnf = qk.encoding.QASMparser(file_name, translate_ccx = True)
                print(f"%%% [quokka# @ {name}] Encoding the circuit into CNF...", flush = True)
                cnf = qk.encoding.QASM2CNF(circuit_cnf, computational_basis = False)
                print(f"%%% [quokka# @ {name}] Setting initial value to |0...0>...", flush = True)
                cnf.leftProjectAllZero()
                print(f"%%% [quokka# @ {name}] Adding measurement to circuit...", flush = True)
                cnf.add_measurement({0:0})
                enc_time = process_time() - enc_time
                ## Executing the circuit one time
                print(f"%%% [quokka# @ {name}] Simulating the circuit...", flush = True)
                ctime = process_time()
                qk.Simulate(cnf)
                ctime = process_time()-ctime
        except TimeoutError:
            print(f"%%% [quokka# @ {name}] Timeout reached for execution", flush = True)
            enc_time, ctime = inf, inf
        memory = tracemalloc.get_traced_memory()[1] / (2**20) # maximum memory usage in MB
        tracemalloc.stop()

    print(f"%%% [quokka# @ {name}] Storing the data...", flush = True)
    result_file.writerow(generate_data(experiment, iterations, enc_time, ctime, memory))

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
    print(f"%%% [full-direct @ {name}] Computing Direct CLUE iterations ({iterations}) for {name} and arguments {args} and {kwds}...", flush=True)
    experiment = generate_example(name, *args, **kwds)
    
    print(f"%%% [full-direct @ {name}] Computing the lumped system...", flush=True)
    tracemalloc.start()
    try:
        with(Timeout(timeout)):
            lump_time = process_time()
            _, U = experiment.direct()
            lump_time = process_time()-lump_time

            print(f"%%% [full-direct @ {name}] Getting the reduced U_P", flush=True)
            U_P = U.to_numpy(dtype=cdouble)
            print(f"%%% [full-direct @ {name}] Getting the reduced U_B", flush=True)
            try:
                U_B = experiment.matrix_B(U_P)
            except:
                print(f"%%% [full-direct @ {name}] No reduced U_B: going to identity", flush=True)
                U_B = eye(U_P.shape[0])
            
            print(f"%%% [full-direct @ {name}] Computing the iteration (U_P*U_B)^iterations", flush=True)
            U = matmul(U_P, U_B)
            it_time = process_time()
            _ = matrix_power(U, iterations)
            it_time = process_time() - it_time
    except TimeoutError:
        print(f"%%% [full-direct @ {name}] Timeout reached for execution", flush=True)
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
            if argv[ind] == "full_quokka#": return "full_quokka#", quokka_iteration

        raise TypeError("Script argument [-t] not well used: we required a follow-up name in ('clue','ddsim','direct','full_clue','full_ddsim','full_direct','full_quokka#')")
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
                methods: list[Callable] | tuple[Callable],  # list of methods with 4 methods to generate the appropriate script
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
                        print(f"### ++ Starting execution {execution}/{repeats} ({size=}, observable={i+1}/{len(obs_to_use)})")
                        if not ("full" in ttype): # No iterations are required
                            used_time = script(*(script_args + [size] + ([observable] if observables != None else [])), timeout=rem_timeout, **kwds)
                            rem_timeout = get_rem_timeout(rem_timeout, used_time)
                            result_file.flush()
                        else:
                            used_time = 0
                            for it in (1,ceil(sqrt(2**size))):#,1000):#,10000)
                                print(f"    ++++ Case with {it} iterations.")
                                it_used_time = script(*(script_args + [it, size] + ([observable] if observables != None else [])), timeout=rem_timeout, **kwds)
                                rem_timeout = get_rem_timeout(rem_timeout, it_used_time)
                                used_time += it_used_time
                            result_file.flush()
                        print(f"### -- Finished execution {execution}/{repeats} ({size=}, observable={i+1}/{len(obs_to_use)}): took {used_time} s.")
                        if rem_timeout != None:
                            rem_timeout -= used_time
                            if rem_timeout < 0:
                                break
                            else:
                                rem_timeout = ceil(rem_timeout)
                except TimeoutError:
                    print(f"### -- Finished execution {execution}/{repeats} ({size=}): reached Timeout.")
