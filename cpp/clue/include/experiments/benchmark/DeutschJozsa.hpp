#ifndef EXP_BEN_DJ
#define EXP_BEN_DJ

#include "Types.hpp"
#include "QuantumComputation.hpp"

using namespace clue;

namespace dj { // Completed
    void dj_oracle(qc::QuantumComputation& circuit);
    void dj_algorithm(qc::QuantumComputation& circuit);
    qc::QuantumComputation* create(luint qubits);
}

#endif