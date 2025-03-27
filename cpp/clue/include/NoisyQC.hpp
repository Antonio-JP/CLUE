#pragma once

#include "QuantumComputation.hpp"

typedef long unsigned int luint;
/*************************************************************************/

/* Class for noisy quantum circuits
 */
class NoisyQuantumComputation
{
protected:
    luint nQubits;
    std::vector<std::tuple<qc::QuantumComputation, double>> layers;

public:
    /* Constructor*/
    NoisyQuantumComputation(luint);

    /* Helper functions*/
    luint size() { return this->nQubits; }
    qc::QuantumComputation build_noisy_qc();          // Build a noisy quantum circuit based on the
    qc::QuantumComputation build_non_noisy_circuit(); // Get the quantum circuit as if no noise is present.

    void push_back(const qc::QuantumComputation qc, const double epsilon)
    {
        if (qc.size() != this->size())
            std::cerr << "The same if the circuit does not match the expected size." << std::endl;

        layers.push_back({qc, epsilon});
    }
};