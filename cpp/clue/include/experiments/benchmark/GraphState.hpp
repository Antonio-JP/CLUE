#ifndef EXP_BEN_GSTATE
#define EXP_BEN_GSTATE

#include "Types.hpp"
#include "Linalg.hpp"
#include "QuantumComputation.hpp"

using namespace std;
using namespace clue;

namespace graphstate { // Completed
    void apply_graph(vector<CCSparseVector>& A, qc::QuantumComputation& circuit);

    qc::QuantumComputation* create(luint qbits);
}

#endif