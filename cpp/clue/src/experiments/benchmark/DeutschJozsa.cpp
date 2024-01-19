#include "experiments/benchmark/DeutschJozsa.hpp"

namespace dj { // Completed
    void dj_oracle(qc::QuantumComputation& circuit) {
        luint qbits = circuit.getNqubits();
        circuit.x(static_cast<qc::Qubit>(qbits-1));
    }
    void dj_algorithm(qc::QuantumComputation& circuit) {
        qc::Qubit n = static_cast<qc::Qubit>(circuit.getNqubits());

        circuit.x(n-1);
        circuit.h(n-1);

        for (qc::Qubit i = 0; i < n-1; i++) {
            circuit.h(i);
        }

        dj_oracle(circuit);

        for (qc::Qubit i = 0; i < n-1; i++) {
            circuit.h(i);
        }
    }
    qc::QuantumComputation* create(luint qubits) {
        qc::QuantumComputation* result = new qc::QuantumComputation(qubits);
        dj_algorithm(*result);
        return result;
    }
}