#include "NoisyQC.hpp"

NoisyQuantumComputation::NoisyQuantumComputation(luint _nQubits)
{
    this->nQubits = _nQubits;
}

qc::QuantumComputation NoisyQuantumComputation::build_noisy_qc()
{
    auto noisy_qc = qc::QuantumComputation(this->nQubits);

    // Generate random number between [0,1)
    std::random_device rd;
    std::mt19937 gen(rd());

    for (const auto &[layer, epsilon] : this->layers)
    {
        if (std::generate_canonical<double, 10>(gen) > epsilon)
        {
            noisy_qc.emplace_back(layer.front()->clone());
        }
    }

    return noisy_qc;
}

qc::QuantumComputation NoisyQuantumComputation::build_non_noisy_circuit()
{
    auto qc = qc::QuantumComputation(this->size());

    for (const auto &[layer, _] : this->layers)
    {
        qc.emplace_back(layer.front()->clone());
    }

    return qc;
}