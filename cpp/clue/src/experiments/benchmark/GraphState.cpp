#include "experiments/benchmark/GraphState.hpp"

#include "experiments/CUTExperiment.hpp"

namespace graphstate { // Completed
    void apply_graph(vector<CCSparseVector>& A, qc::QuantumComputation& circuit) {
        qc::Qubit qbits = static_cast<qc::Qubit>(A.size());

        for (qc::Qubit i = 0; i < qbits; i++) {
            circuit.h(i);
        }
        for (qc::Qubit i = 0; i < qbits; i++) {
            for (qc::Qubit j = i+1; i < qbits; j++) {
                if (A[i][j] != CC(0)) {
                    circuit.cz(i, j);
                }
            }
        }
    }

    qc::QuantumComputation* create(luint qbits) {
        UndirectedGraph* G = UndirectedGraph::random(qbits, 0.3);
        vector<CCSparseVector> A = G->adjacency_matrix();

        qc::QuantumComputation* circuit = new qc::QuantumComputation(qbits);
        apply_graph(A, *circuit);

        return circuit;
    }
}