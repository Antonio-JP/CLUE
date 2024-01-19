#ifndef EXP_BEN_GHZ
#define EXP_BEN_GHZ

#include "Types.hpp"
#include "QuantumComputation.hpp"

using namespace clue;

namespace ghz {
    qc::QuantumComputation* create(luint qubits);
}

#endif