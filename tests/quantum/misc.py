r'''
    Some auxiliary methods
'''
from qiskit.circuit import QuantumCircuit

def loop(circuit: QuantumCircuit, size: int, iterations: int, prepend_H: bool = True, measure: bool = False):
    r'''
        Creates a quantum circuit as a loop of a circuit with fixed number of iterations
    '''
    final = QuantumCircuit(size)
    bits = list(range(size))

    if prepend_H:
        for b in bits:
            final.h(b)
    
    for _ in range(iterations):
        final.append(circuit, bits)

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
