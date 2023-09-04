import sys, os

SCRIPT_DIR = os.path.dirname(__file__) if __name__ != "__main__" else "./"
sys.path.insert(0, os.path.join(SCRIPT_DIR, "..", "..")) # clue is here

import logging, csv, re

from itertools import chain, combinations
from clue.clue import LDESystem
from clue.qiskit import *
from math import ceil, gcd, sqrt
from mqt import ddsim
from numpy import cdouble, load, save
from numpy.linalg import eig
from os import listdir
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, execute, Aer, qpy
from time import time
from random import choice, randint

logger = logging.getLogger("clue")

########################################################################
### AUXILIARY METHODS
########################################################################
def add_first(first, gen):
    yield first
    yield from gen

def powerset(iterable):
    "powerset([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1,len(s)+1))

def list_circuits(*argv):
    r'''List examples in the folder. It allows several arguments.'''
    def get_size(example):
        try:
            return int(example.removesuffix(".qasm").split("_")[-1])
        except:
            return -1

    examples = [file.removesuffix(".qasm") for file in listdir(os.path.join(SCRIPT_DIR, "circuits")) if file.endswith(".qasm")]
    examples.sort(); examples.sort(key=get_size)
    executed_examples = [file.removeprefix("[output]").removesuffix(".example.txt") for file in listdir(os.path.join(SCRIPT_DIR, "results")) if file.endswith(".example.txt")]

    full = False
    allowed_names = []; forbidden_names = []; executed = None; sizes = None
    i = 0
    ## Reading the arguments
    while i < len(argv): 
        if argv[i] == "-n":
            allowed_names.append(argv[i+1]); i += 2
        elif argv[i] == "-wn":
            forbidden_names.append(argv[i+1]); i += 2
        elif argv[i] == "-e":
            if executed != None and executed != True: raise TypeError("The command for executed tests were already given.")
            executed = True; i+= 1
        elif argv[i] == "-ne":
            if executed != None and executed != False: raise TypeError("The command for executed tests were already given.")
            executed = False; i+= 1
        elif argv[i] == "-a":
            full = True; i+= 1
        elif argv[i] == "-s":
            if sizes is None: sizes = []
            sizes.append(int(argv[i+1])); i+= 2
        else:
            raise TypeError(f"Option {argv[i]} not recognized. Check 'help' command for further information")


    ## Creating the filtering function
    filter = lambda example: (
                                (len(allowed_names) == 0 or any(name in example for name in allowed_names)) and 
                                (len(forbidden_names) == 0 or all(not (name in example) for name in forbidden_names)) and
                                (executed == None or ((example in executed_examples) == executed)) and
                                (sizes == None or (get_size(example) in sizes))
                             )

    ## Creating the string to be printed
    if full:
        lines = [["Example name", "# q-bits", "Executed"]]
        get_str = lambda example : (example, str(get_size(example)), "Yes" if example in executed_examples else "No")

        lines.extend([get_str(name) for name in examples if filter(name)])
        lines.append(["N.models", f"{len(lines)-1}", ""])
        n_elem = len(lines[0])
        max_length = [max(len(line[i]) if line[i] != None else 4 for line in lines) for i in range(n_elem)]

        lines.insert(1, [max_length[i]*"-" for i in range(n_elem)])
        lines.insert(len(lines)-1, [max_length[i]*"-" for i in range(n_elem)])

        for line in lines:
            print(" | ".join([(line[i] if line[i] != None else "None").ljust(max_length[i]) for i in range(n_elem)]))
    else:
        print(" ".join([name for name in examples if filter(name)]))

########################################################################
### METHODS FOR COMPRESSING
########################################################################
def reduce_circuit(circuit: str, observable: list[str] = None, **kwds):
    with open(os.path.join(SCRIPT_DIR, "results", f"[compress]{circuit}.example.txt"), "w") as out_file:
        out_file.write("########################################################################\n")
        out_file.write(f"### Compressing {circuit} with observable {observable}\n")
        out_file.write("########################################################################\n")
        out_file.write("*** Reading the Dynamical System...\n"); ctime = time()
        ## Reading the ".qasm" circuit
        system = DS_QuantumCircuit.from_qasm_file(os.path.join(SCRIPT_DIR, "circuits", f"{circuit}.qasm"), delta=kwds.pop("delta", 1e-10))
        read_time=time()-ctime; 
        ## Computing the lumping for the first state
        out_file.write("*** Computing lumping...\n"); ctime = time()
        observable = [system.variables[0]] if observable is None else observable
        lumped = system.lumping(observable, file=out_file)
        lumping_time=time()-ctime
        ## Writing the reduced circuit
        out_file.write("*** Building reduces .qpy...\n"); ctime = time()
        Ur = lumped.construct_matrices("polynomial")[0].to_numpy(dtype=cdouble); Ur
        Ur = extend_to_power(Ur)
        nqbits = int(log2(Ur.shape[0]))
        q = QuantumRegister(nqbits,'q')
        c = ClassicalRegister(nqbits,'c')
        red_circuit = QuantumCircuit(q,c)
        red_circuit.unitary(Ur, q)
        red_circuit.measure(q, c)
        with open(os.path.join(SCRIPT_DIR, "circuits", "reduced", f"{circuit}.qpy"), "wb") as f:
            qpy.dump(red_circuit, f)
        reduced_time = time()-ctime
        ## Writing the lumping matrix
        out_file.write("*** Saving the lumping data...\n"); ctime = time()
        save(os.path.join(SCRIPT_DIR, "circuits", "reduced", f"{circuit}.npy"), lumped.lumping_matrix.to_numpy(dtype=cdouble))
        matrix_time = time()-ctime

        out_file.write("########################################################################\n")
        out_file.write("*** Execution data:\n")
        out_file.write(f"*** Time reading: {read_time}\n")
        out_file.write(f"*** Time lumping: {lumping_time}\n")
        out_file.write(f"*** Time saving the reduced circuit: {reduced_time}\n")
        out_file.write(f"*** Time saving the lumping: {matrix_time}\n")
        out_file.write(f"*** Total time: {sum((read_time,lumping_time, reduced_time, matrix_time))}\n")
        out_file.write(f"*** Size of reduction (in qbits): {nqbits} (from {int(log2(system.size))})\n")
        out_file.write("########################################################################\n")
        
        return {"read": read_time, "lumping": lumping_time, "reduced": reduced_time, "matrix": matrix_time}

def simulate_qasm(circuit, backend="aer", shots=8192, output="counts"):
    ctime = time()
    ## Reading the quantum circuit
    circuit = QuantumCircuit.from_qasm_file(os.path.join(SCRIPT_DIR, "circuits",f"{circuit}.qasm"))
    reading_time = time()-ctime; print("Finished reading the circuit"); ctime=time()
    ## Simulating the quantum circuit
    backend = Aer.get_backend('aer_simulator') if backend == "aer" else ddsim.DDSIMProvider().get_backend("qasm_simulator") if backend == "ddsim" else backend
    job = execute(circuit, backend, shots=shots)
    simulation_time = time()-ctime; print("Finished executing the job"); ctime=time()
    ## Collecting the data to have proper output
    output = job.result().get_counts() if output == "counts" else job if output == "job" else job.result()
    output_time = time()-ctime; print("Finished obtaining results")
    
    return output, {"read": reading_time, "simulation": simulation_time, "output": output_time}, reading_time + simulation_time + output_time

def simulate_reduced(circuit_name, shots=8192):
    ctime = time()
    ## Reading the quantum circuit
    with open(os.path.join(SCRIPT_DIR, "circuits","reduced", f"{circuit_name}.qpy"), "rb") as f:
        circuit = qpy.load(f)[0]
    reading_time = time()-ctime; ctime=time()    
    ## Simulating the quantum circuit
    backend = Aer.get_backend('aer_simulator')
    job = execute(circuit, backend, shots=shots)
    simulation_time = time()-ctime; ctime=time()
    ## Collecting the data to have proper output
    #### Reading the lumping matrix
    Ur = load(os.path.join(SCRIPT_DIR, "circuits","reduced", f"{circuit_name}.npy"))
    nqbits = int(log2(Ur.shape[1]))
    matrix_time = time()-ctime; ctime=time()
    #### Measuring the data from qasm to the original states
    output = dict()
    for (out, times) in job.result().get_counts().items():
        for m,t in repeated_measure(Ur[int(out, 2)], times, out=dict).items():
            output[m] = output.get(m,0) + t
    output = {format(k, f"0{nqbits}b"): v for k,v in output.items()}
    measuring_time = time()-ctime
    output_time = matrix_time + measuring_time
    
    return output, {"read": reading_time, "simulation": simulation_time, "output": {"total": output_time, "matrix": matrix_time, "measure": measuring_time}}

def run_example(circuit: str, **kwds):
    system = DS_QuantumCircuit.from_qasm_file(os.path.join(SCRIPT_DIR, "circuits", f"{circuit}.qasm"), delta=kwds.pop("delta", 1e-10))
        
    with open(os.path.join(SCRIPT_DIR, "results", f"[output]{circuit}.example.txt"), "w") as out_file:
        if "observable" in kwds:
            first_obs = None
            full_observables = list()
            for obs in kwds["observable"]:
                observable = list()
                for part in obs.split("+"):
                    if part.startswith("|"):
                        observable.append(system.variables[int(part.removeprefix("|").removesuffix(">"))])
                    elif part.startswith("Q"):
                        observable.append(part)
                    elif part.startswith("H"):
                        observable.append("+".join(system.variables))
                    else:
                        raise NotImplementedError(f"Processing observable part {part} not implemented.")
                observable = "+".join(observable)
                full_observables.append(observable)
                generator = [tuple(full_observables)]
                all_obs = 1
        else:
            first_obs = ("+".join(system.variables),)
            all_obs = 2**len(system.variables); generator = powerset(system.variables)
            if all_obs > 10000: 
                all_obs = len(system.variables)+1
                generator = (tuple([el]) for el in system.variables)
            generator = add_first(first_obs, generator)
        
        try:
            for i,obs in enumerate(generator):
                out_file.write("################################################################\n")
                out_file.write(f"### Observable: {('H-State',) if obs == first_obs else obs}\n")
                out_file.write("################################################################\n")
                logger.info(f"[circuit_example ({i+1}/{all_obs})] Lumping {circuit} with observable {obs} ")
                start = time()
                lumped = system.lumping(obs,print_reduction=True,file=out_file)
                total = time()-start
                out_file.write(f"**  Reduction found: {system.size} --> {lumped.size}\n")
                out_file.write(f"**  Time spent: {total} s.\n")
                out_file.write("################################################################\n")
                out_file.flush()
        except KeyboardInterrupt:
            logger.info(f"[circuit_example] Interrupted {circuit} with Ctr-C")

    return

def __is_prime(n):
    return all(n%i != 0 for i in primes(int(ceil(sqrt(n)+1))))

__CACHED_PRIMES = [2,3,5,7,11,13]
@lru_cache(maxsize=None)
def primes(bound):
    if bound > __CACHED_PRIMES[-1]:
        for n in range(__CACHED_PRIMES[-1]+1, bound):
            if __is_prime(n): __CACHED_PRIMES.append(n)
            
    li = 0; ui = len(__CACHED_PRIMES) - 1
    if __CACHED_PRIMES[-1] <= bound: return __CACHED_PRIMES[:len(__CACHED_PRIMES)]
    while (ui-li) > 1:
        ci = (ui+li)//2
        if __CACHED_PRIMES[ci] < bound: li = ci
        elif __CACHED_PRIMES[ci] > bound: ui = ci
        else: ui = ci; break # found the element
    return __CACHED_PRIMES[:ui]

def __get_matrix_case_study(case: str, qbits: int):
    r'''Generate a unitary matrix for the case study and the data to reproduce the experiment'''
    if case == "search":
        # we generate the search function
        success = [randint(0, 2**qbits-1) for _ in range(qbits)] # sparse success search, but random
        f = lambda p : 1 if p in success else 0
        # we build the matrix from Grover's algorithm
        matrix = G(f, qbits)
        # we build the entangled state as observable
        observables = [("H-state", [[1. for _ in range(2**qbits)]])]
        data = {"success": success}
    elif case in ("order", "phase"):
        ## We generate a product of two primes that are smaller than 2**qbits
        logger.info(f"[get_matrix @ {case} - {qbits}] Trying to get the product of two different primes")
        candidates = [(p,q) for p,q in product(primes(ceil(2**qbits / 3))[1:], repeat=2) if (p < q and p*q >= 2**(qbits-1) and p*q < 2**qbits)]
        if len(candidates) == 0:
            raise ValueError(f"Impossible thing happened: no composed number between {2**(qbits-1)} and {2**qbits}")
        p,q = choice(candidates)
        N = p*q
                    
        logger.info(f"[get_matrix @ {case} - {qbits}] {N} = {p} * {q}")
        ## We generate a coprime number ´x´
        x = randint(2, N-1)
        while gcd(N, x) != 1: x = randint(2, N-1)

        ## We create the circuit for order finding of multiplying by x module N
        matrix = U(x, N)

        ## Now for each case we distinguish the data
        if case == "order":
            observables = [("|1>", [[0,1] + (2**qbits - 2)*[0]])]
            data = {"N": N, "x": x, "p": p, "q": q}
        else: # case is "phase"
            ## We compute a non-trivial eigenvector for the matrix
            vals, vects = eig(matrix)
            # we filter eigenvalues 1 and -1 (boring)
            valid = [i for i in range(len(vals)) if vals[i].round(5) not in (0., 1., -1.)]
            if len(valid) == 0: valid = [i for i in range(len(vals)) if vals[i].round(5) not in (0., 1.)] # if not possible, we allow -1 as eigenvalue
            i = choice(valid)
            eigenvalue = vals[i]
            u = vects[:,i]

            ## We now change to Kitaev's gate:
            matrix = K(matrix)
            observables = [("|0>|u>", [kron([1,0], u)]), ("|1>|u>", [kron([0,1], u)]), ("H|u>", [kron([1/sqrt(2), -1/sqrt(2)], u)])]
            data = {"N": N, "x": x, "u": u, "lambda": eigenvalue}
    else:
        raise NotImplementedError(f"The case {case} is not recognized")
    return matrix, observables, data

def __write_extra_data(case: str, qbits: int, observable, data, out_file):
    r'''Write in the result file the data corresponding to the given case study'''
    out_file.write(f"### Size: {qbits}\n")
    out_file.write(f"### Type: {case}\n")
    if case == "search":
        out_file.write(f"### Success search: {data['success']}\n")
    elif case == "order":
        out_file.write(f"### Multiplication by: {data['x']}\n")
        out_file.write(f"### Modulo: {data['N']}\n")
        out_file.write(f"### Factors: {data['p']}, {data['q']}\n")
    elif case == "phase":
        out_file.write(f"### Multiplication by: {data['x']}\n")
        out_file.write(f"### Modulo: {data['N']}\n")
        out_file.write(f"### Eigenvalue: {data['lambda']}\n")
    else:
        raise NotImplementedError

def __post_lumping_study(case: str, system: DS_QuantumCircuit, lumped: LDESystem, observable: Sequence[SparsePolynomial], data, out_file):
    r'''Method to do extra computations and print extra results depending on the case study'''
    if case == "search":
        logger.info(f"[case_study ({case})] Checking coherence of result")
        if lumped.size != 2:
            out_file.write(f"** Size error: {lumped.size}\n")
        else:
            logger.info(f"[case_study ({case})] Checking recovering information")
            L = lumped.lumping_matrix.to_numpy(dtype=cdouble)
            new_L = [L[0], L[1]]
            first_nonzero = min([i for i in range(len(new_L[0])) if new_L[0][i] != 0])
            new_L[1] = (new_L[1] - new_L[1][first_nonzero]/new_L[0][first_nonzero]*new_L[0]).round(10)
            first_nonzero = min([i for i in range(len(new_L[1])) if new_L[1][i] != 0])
            new_L[0] = new_L[0] - new_L[0][first_nonzero]/new_L[1][first_nonzero]*new_L[1]

            states = [set([i for i in range(len(row)) if row[i].round(10) != 0]) for row in new_L]
            success = set(data["success"])
            if success in states:
                out_file.write(f"** Success?: {success in states}\n")
            else:
                out_file.write(f"** Success?: {success in states} --> {states}\n")
    elif case == "order":
        N = data["N"]; x = data["x"]
        if x**lumped.size % N != 1:
            out_file.write(f"** Size error: {lumped.size} --> {x**lumped.size % N}\n")
        else:
            out_file.write(f"** Success?: True\n")
            if lumped.size % 2 == 0:
                guess = gcd(N, x**(lumped.size//2)+1)
                p, q = (guess, N//guess) if guess != 1 else (1,1)
                
                out_file.write(f"** Factorized?: {p*q == N}\n")
    elif case == "phase":
        N = data["N"]; x = data["x"]; u = data["u"]; eigenvalue = data["lambda"]
        ## getting the observable we are using
        if observable[0].linear_part_as_vec().nonzero_count() == len([i for i in range(len(u)) if u[i] != 0]): # it is |j> |u> for j in {0,1}
            if lumped.size != 2:
                out_file.write(f"** Size error: {lumped.size}\n")
            else:
                out_file.write(f"** Success?: True\n")
                ## Add recovering of phase by measuring the reduced system
        else: # it is the special case with lumping of size 1
            if lumped.size != 1:
                out_file.write(f"** Size error: {lumped.size}\n")
            else:
                out_file.write(f"** Success?: True\n")
    else:
        raise NotImplementedError

def run_case_study(case: str, qbits: int, repeats: int = 1):
    r'''
        Run the case studies with a given number of q-bits and repeat the experiment a given number of times
    '''
    ## Checking the arguments
    if not case in ("search", "order", "phase"):
        raise ValueError("We only allow 3 case studies: 'search' for Grover's circuit, 'order' for order finding circuit and 'phase' for Kitaev's circuit.")
    if not isinstance(qbits, int) or qbits <= 0:
        raise TypeError(f"The number of q-bits must be a positive integer. Got {qbits}")
    if not isinstance(repeats, int) or repeats <= 0:
        raise TypeError(f"The number of repetitions must be a positive integer. Got {qbits}")

    with open(os.path.join(SCRIPT_DIR, "results", f"[output]case_{case}[{qbits=}].example.txt"), "w") as out_file:
        for _ in range(repeats):
            matrix, observables, data = __get_matrix_case_study(case, qbits)

            system = DS_QuantumCircuit(matrix)

            try:
                for obs in observables:
                    if isinstance(obs, tuple) and isinstance(obs[0], str):
                        name, obs = obs[0], obs[1]
                    else:
                        name = str(obs)
                    obs = [SparsePolynomial.from_vector(v, system.variables, system.field) for v in obs]
                    out_file.write("################################################################\n")
                    out_file.write(f"### Observable: {name}\n")
                    __write_extra_data(case, qbits, obs, data, out_file)
                    out_file.write("################################################################\n")
                    logger.info(f"[case_study ({case}:{qbits})] Lumping with observable {name}")
                    start = time()
                    lumped = system.lumping(obs,print_reduction=True,file=out_file)
                    total = time()-start
                    out_file.write(f"**  Reduction found: {system.size} --> {lumped.size}\n")
                    out_file.write(f"**  Time spent: {total} s.\n")
                    __post_lumping_study(case, system, lumped, obs, data, out_file)
                    out_file.write("################################################################\n")
                    out_file.flush()
            except KeyboardInterrupt:
                logger.info(f"[case_study] Interrupted {circuit} with Ctr-C")

def compile_results():
    data = list()
    for file_name in listdir(os.path.join(SCRIPT_DIR, "results")):
        if file_name.startswith("[output]"):
            with open(os.path.join(SCRIPT_DIR, "results", file_name), "r") as file:
                circuit = file_name.removeprefix("[output]").removesuffix(".example.txt")
                ## We have two options for the circuit: a case study or a benchmark example
                if "_qiskit_" in circuit: # these are benchmark examples
                    qbits = int(circuit.split("_")[-1]); circuit = circuit.removesuffix(f"_{qbits}")
                else: # case study
                    qbits = int(circuit.split("[qbits=")[1].removesuffix("]"))
                    circuit = circuit.split("[")[0]
                logger.log(60, f"[compile] Starting reading result for {circuit} ({qbits})")
                line = file.readline().strip() # first line of file
                while line != "":
                    if line.startswith("################################################################"):
                        # We start an example
                        observable = file.readline().strip().removeprefix("### Observable: ").removeprefix("(").removesuffix(")").removeprefix("[").removesuffix("]").removesuffix(",")
                        observable = re.sub("Q_[01]*", lambda match : f"|{int(match.group(0).removeprefix('Q_'), base=2)}>", observable)
                        or_size, red_size, time_used = None, None, None
                        line = file.readline().strip()
                        while line != "################################################################": ## Cleaning until the end of header for example
                            line = file.readline().strip()
                        line = file.readline().strip() # first line of example
                        while line != "################################################################":
                            if line == "": 
                                logger.error(f"[compile - {circuit} ({qbits})] The format of the result file is not correct: EOF unexpected")
                                break
                            if line.startswith("**  Reduction found:"):
                                or_size, red_size = [int(el.strip()) for el in line.removeprefix("**  Reduction found:").split("-->")]
                            if line.startswith("**  Time spent: "):
                                time_used = float(line.removeprefix("**  Time spent: ").removesuffix("s."))

                            line = file.readline().strip() # next line
                        else:
                            data.append([circuit, qbits, observable.replace("', '", " & "), or_size, red_size, red_size/or_size, time_used])
                            line = file.readline().strip() # next line of file
                            continue
                        break
                    line = file.readline().strip() # next line of file
                
                logger.log(60, f"[compile] Finished results for {circuit} ({qbits})")
    
    ## We have ompiled all the data: we create a CSV for it
    
    logger.log(60, f"[compile] Dumping to CSV file")
    with open(os.path.join(SCRIPT_DIR, "compilation.csv"), "w") as file:
        writer = csv.writer(file, delimiter=",")
        writer.writerow(["circuit", "qbits", "observable", "size", "reduced", "ratio", "time"]) # Header
        for result in data:
            writer.writerow(result)

    return

if __name__ == "__main__":
    ## Reading the arguments
    if sys.argv[1] == "list":
        list_circuits(*sys.argv[2:])
    elif sys.argv[1] == "compile":
        compile_results()
    elif sys.argv[1] == "case":
        n = 2; nargs = len(sys.argv); cases = set(); sizes = set(); repeats = None
        while n < nargs:
            if sys.argv[n].startswith("-"):
                if sys.argv[n].endswith("s"):
                    sizes.add(int(sys.argv[n+1])); n += 2
                elif sys.argv[n].endswith("r"):
                    repeats = int(sys.argv[n+1]); n += 2
                else: n += 1
            else:
                cases.add(sys.argv[n]); n += 1

        if repeats == None: repeats = 10
        if len(sizes) == 0: sizes = set([5,6])
        if len(cases) == 0: cases = set(["search", "order", "phase"])

        for case, size in product(cases, sizes):
            run_case_study(case, size, repeats)
    elif sys.argv[1] == "compress":
        for circuit in sys.argv[2:]:
            logger.info(f"[compress] Compressing circuit {circuit}...")
            try:
                reduce_circuit(circuit)
            except Exception as e:
                logger.error(f"[compress] Error in circuit {circuit}: {e}")
            logger.info(f"[compress] Compressing circuit {circuit}... Done")
    else:
        n = 1; nargs = len(sys.argv); kwds = dict(); circuits = list()
        while n < nargs:
            if sys.argv[n].startswith("-"):
                if sys.argv[n].endswith("d"):
                    kwds["delta"] = float(sys.argv[n+1])
                    n += 2
                elif sys.argv[n].endswith("o"):
                    if "observable" not in kwds: kwds["observable"] = list()
                    kwds["observable"].append(sys.argv[n+1])
                    n += 2
                else:
                    n+=1
            else:
                circuits.append(sys.argv[n].removeprefix("circuits/").removesuffix(".qasm")); n+= 1
            
        for circuit in circuits:
            logger.info(f"[circuit_example] Running example from circuit {circuit}")
            run_example(circuit, **kwds)
            logger.info(f"[circuit_example] Finished example from circuit {circuit}")