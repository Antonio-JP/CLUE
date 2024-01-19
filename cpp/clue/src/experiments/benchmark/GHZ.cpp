#include "experiments/benchmark/GHZ.hpp"

namespace ghz {
    qc::QuantumComputation* create(luint qubits) {
        qc::QuantumComputation* circuit = new qc::QuantumComputation(qubits);
        qc::Qubit last = static_cast<qc::Qubit>(qubits-1);
        circuit->h(last);
        for (qc::Qubit i = 1; i < qubits; i++) {
            circuit->cx(last - i + 1, last - i);
        }

        return circuit;
    }
}