import sys, os

SCRIPT_DIR = os.path.dirname(__file__) if __name__ != "__main__" else "./"
sys.path.insert(0, os.path.join(SCRIPT_DIR, "..", "..")) # clue is here

import logging, csv

from itertools import chain, combinations
from clue.qiskit import *
from os import listdir
from time import time

logger = logging.getLogger("clue")

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

def run_example(circuit: str, **kwds):
    system = DS_QuantumCircuit(os.path.join(SCRIPT_DIR, "circuits", f"{circuit}.qasm"), delta=kwds.pop("delta", 1e-10))

    with open(os.path.join(SCRIPT_DIR, "results", f"[output]{circuit}.example.txt"), "w") as out_file:
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
            logger.info(f"[circuit_example] Interrupted {circuit} with Ctr-C ")

    return

def compile_results():
    data = list()
    for file_name in listdir(os.path.join(SCRIPT_DIR, "results")):
        if file_name.startswith("[output]"):
            with open(os.path.join(SCRIPT_DIR, "results", file_name), "r") as file:
                circuit = file_name.removeprefix("[output]").removesuffix(".example.txt")
                logger.log(60, f"[compile] Starting reading result for {circuit}")
                line = file.readline().strip() # first line of file
                while line != "":
                    if line.startswith("################################################################"):
                        # We start an example
                        observable = file.readline().strip().removeprefix("### Observable: (").removesuffix(")").removesuffix(",")
                        or_size, red_size, time_used = None, None, None
                        file.readline() # line of #####
                        line = file.readline().strip() # first line of example
                        while line != "################################################################":
                            if line == "": 
                                logger.error(f"[compile - {circuit}] The format of the result file is not correct: EOF unexpected")
                                break
                            if line.startswith("**  Reduction found:"):
                                or_size, red_size = [int(el.strip()) for el in line.removeprefix("**  Reduction found:").split("-->")]
                            if line.startswith("**  Time spent: "):
                                time_used = float(line.removeprefix("**  Time spent: ").removesuffix("s."))

                            line = file.readline().strip() # next line
                        else:
                            data.append([circuit, observable.replace("', '", " & "), or_size, red_size, red_size/or_size, time_used])
                            line = file.readline().strip() # next line of file
                            continue
                        break
                    line = file.readline().strip() # next line of file
                
                logger.log(60, f"[compile] Finished results for {circuit}")
    
    ## We have ompiled all the data: we create a CSV for it
    
    logger.log(60, f"[compile] Dumping to CSV file")
    with open(os.path.join(SCRIPT_DIR, "compilation.csv"), "w") as file:
        writer = csv.writer(file, delimiter=",")
        writer.writerow(["circuit", "observable", "size", "reduced", "ratio", "time"]) # Header
        for result in data:
            writer.writerow(result)

    return

if __name__ == "__main__":
    ## Reading the arguments
    if sys.argv[1] == "list":
        list_circuits(*sys.argv[2:])
    elif sys.argv[1] == "compile":
        compile_results()
    else:

        n = 1; nargs = len(sys.argv); kwds = dict(); circuits = list()
        while n < nargs:
            if sys.argv[n].startswith("-"):
                if sys.argv[n].endswith("d"):
                    kwds["delta"] = float(sys.argv[n+1])
                    n += 2
                else:
                    n+=1
            else:
                circuits.append(sys.argv[n].removeprefix("circuits/").removesuffix(".qasm")); n+= 1
            
        for circuit in circuits:
            logger.info(f"[circuit_example] Running example from circuit {circuit}")
            run_example(circuit, **kwds)
            logger.info(f"[circuit_example] Finished example from circuit {circuit}")