#include <iostream>
#include <string>

#include "NoisyQC.hpp"

// Test if a layer is added to the noisy_qc currectly with the push_back function
void add_layer_succes(int n, luint qbits)
{

    auto nqc = NoisyQuantumComputation(qbits);

    // Add layers to the Noisy QuantumComputation
    for (int i = 0; i < n; i++)
    {
        auto layer = qc::QuantumComputation(qbits);
        layer.h(0);
        double epsilon = 0.001;
        nqc.push_back(layer, epsilon);
    }

    // Throw error if we do not have the expected number of layers
    if (nqc.size() != n)
    {
        throw std::runtime_error("The noisy quantum computation layers did not get added currectly. Expected " + std::to_string(n) + " got " + std::to_string(nqc.size()));
    }
}

void add_layer_failure_incorrect_size()
{
    luint nqc_size = 3;
    luint qc_size = 4;
    double epsilon = 0.01;

    auto nqc = NoisyQuantumComputation(nqc_size);
    auto invalid_qc = qc::QuantumComputation(qc_size);

    try
    {
        nqc.push_back(invalid_qc, epsilon);
    }
    catch (std::runtime_error)
    {
        return;
    }
    throw std::runtime_error("Expected runtime_error was not thrown for mismatched sizes.");
}

void add_layer_failure_invalid_epsilon()
{
    luint nqc_size = 3;
    double invalid_epsilon = 1.5; // Invalid epsilon value (greater than 1.0)   throw std::runtime_error("Incorrect size of Noisy QC and QC didn't throw error");

    auto nqc = NoisyQuantumComputation(nqc_size);
    auto qc = qc::QuantumComputation(nqc_size);

    try
    {
        nqc.push_back(qc, invalid_epsilon);
    }
    catch (std::runtime_error)
    {
        return;
    }
    throw std::runtime_error("Expected runtime_error was not thrown for invalid epsilon.");
}

void add_layer_failure_empty_qc()
{
    luint nqc_size = 3;
    double epsilon = 0.01;

    auto nqc = NoisyQuantumComputation(nqc_size);
    auto empty_qc = qc::QuantumComputation(nqc_size);

    try
    {
        nqc.push_back(empty_qc, epsilon);
    }
    catch (std::runtime_error)
    {
        return;
    }
    throw std::runtime_error("Expected runtime_error was not thrown for empty quantum circuit.");
}

void build_non_noisy_circuit(luint qubits)
{
    double epsilon = 0.01;

    auto nqc = NoisyQuantumComputation(qubits);

    // Add layers to the Noisy QuantumComputation
    for (int i = 0; i < qubits; i++)
    {
        auto qc = qc::QuantumComputation(qubits);
        qc.h(i);
        nqc.push_back(qc, epsilon);
    }

    // Build the non-noisy circuit
    auto non_noisy_qc = nqc.build_non_noisy_qc();

    // Check if the size of the non-noisy circuit is equal to the number of layers
    if (non_noisy_qc.size() != nqc.size())
    {
        throw std::runtime_error("The size of the non-noisy circuit does not match the number of layers.");
    }
}

int main()
{
    add_layer_succes(3, 3);
    add_layer_failure_incorrect_size();
    add_layer_failure_invalid_epsilon();
    add_layer_failure_empty_qc();
    build_non_noisy_circuit(3);

    return 0;
}