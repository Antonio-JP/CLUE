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
    int size() { return this->layers.size(); }
    luint getNqubits() { return this->nQubits; } // Number of qubits in the circuit
    qc::QuantumComputation *build_noisy_qc();    // Build a noisy quantum circuit based on the epsilon values. Have to return a pointer to match with experiment class.
    qc::QuantumComputation build_non_noisy_qc(); // Get the quantum circuit as if no noise is present.

    /* Add a layer to the Noisy Quantum Computation
       Primary use of the noisy quantum computation. With this implementation we can each gate a specific epsilon value.

       This is a bit more cumbersome of an implementation, but it allows us to build it in similar fashion to the python implementation.
    */
    void push_back(const qc::QuantumComputation qc, const double epsilon)
    {
        if (qc.size() == 0)
            throw std::runtime_error("The quantum circuit is empty.");
        if (qc.getNqubits() != this->getNqubits())
            throw std::runtime_error("The qubits of the circuits does not match the expected number.");
        if (epsilon > 1.0 || epsilon < 0.0)
            throw std::runtime_error("Noise should be between 0 and 1.");

        layers.push_back({qc, epsilon});
    }
};