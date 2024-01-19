#ifndef EXP_BEN_QNN
#define EXP_BEN_QNN

#include "Types.hpp"
#include "QuantumComputation.hpp"

using namespace clue;

namespace qnn {
    qc::QuantumComputation* create(luint);
}

#endif