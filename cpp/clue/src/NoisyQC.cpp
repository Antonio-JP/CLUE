#include "NoisyQC.hpp"

NoisyQuantumComputation::NoisyQuantumComputation(luint _nQubits)
{
    this->nQubits = _nQubits;
}

/*  When we build the noisy qc, we simply add the operation with probability 1-epsilon or else we do nothing.
    With this implementation, we are also requiring the qc to have only one operation in the layer.
    This is a bit more cumbersome of an implementation, but it allows us to build it in similar fashion to the python implementation.
    I thought about instead of doing nothing, we apply the identify gate to all the targets in the operation. Thoughts?*/
qc::QuantumComputation *NoisyQuantumComputation::build_noisy_qc()
{
    auto noisy_qc = new qc::QuantumComputation(this->nQubits);

    // Generate random number between [0,1)
    std::random_device rd;
    std::mt19937 gen(rd());

    for (const auto &[layer, epsilon] : this->layers)
    {
        if (std::generate_canonical<double, 10>(gen) > epsilon) // add the layer with probability 1-epsilon. Else do nothing.
        {
            noisy_qc->emplace_back(layer.front()->clone());
        }
    }

    return noisy_qc;
}

qc::QuantumComputation NoisyQuantumComputation::build_non_noisy_qc()
{
    auto qc = qc::QuantumComputation(this->size());

    for (const auto &[layer, _] : this->layers)
    {
        qc.emplace_back(layer.front()->clone());
    }

    return qc;
}