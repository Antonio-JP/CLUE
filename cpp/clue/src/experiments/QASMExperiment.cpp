#include "experiments/QASMExperiment.hpp"

#include <boost/algorithm/string.hpp>
#include <cmath>
#include "CircuitOptimizer.hpp"
#include "dd/Operations.hpp"
#include "dd/FunctionalityConstruction.hpp"

/*** CODE FOR CLASS QASM_EXAMPLE ***/
luint QASMExperiment::correct_size() {
    return 0; // Unknown
}
luint QASMExperiment::bound_size() {
    return this->nstates; // No preknown bound
}

array<dd::CMat, 2U> QASMExperiment::direct() {
    throw logic_error("Direct type not implemented for QASM examples");
}

vector<CCSparseVector> QASMExperiment::matrix() {
    this->quantum(0.0);
    dd::mEdge circuit_dd = buildFunctionality(this->circuit, *this->package);
    dd::CMat unitary = circuit_dd.getMatrix(this->size()); 
    luint N = unitary.size();
    vector<CCSparseVector> U = vector<CCSparseVector>(N, N); // Initializing the final matrix
    for (luint i = 0; i < N; i++) {
        for (luint j = 0; j < N; j++) {
            U[i].set_value(j, unitary[i][j]);
        }
    }

    return U;
}
dd::CMat QASMExperiment::matrix_B(dd::CMat& U) {
    return identity_matrix(U.size()); // There is no begin hamiltonian: we use the identity
}

qc::QuantumComputation* QASMExperiment::quantum(double) {
    if (this->circuit == nullptr) {
        this->circuit = new qc::QuantumComputation(this->path);
        qc::CircuitOptimizer::removeFinalMeasurements(*this->circuit);
    }

    if (this->circuit == nullptr) {
        throw logic_error("Error loading file " + this->path);
    }
    return this->circuit;
}

qc::QuantumComputation* QASMExperiment::quantum_B(double) { // A nothing circuit
    qc::QuantumComputation* empty = new qc::QuantumComputation(this->size());

    return empty;
}

QASMExperiment* QASMExperiment::change_exec_type(ExperimentType new_type) {
    return new QASMExperiment(this->size(), this->name, this->path, this->observable, new_type, this->package);
}

/*** BUILDERS FOR BENCHMARK_EXAMPLE ***/
QASMExperiment::QASMExperiment(luint bQbits, string eName, string ePath, string eObservable, ExperimentType eType, dd::Package<>* ePackage) : 
    Experiment(eName, eObservable, 1UL, eType, ePackage) {
    this->path = ePath;
    this->qbits = bQbits;
    this->circuit = nullptr;
    this->nstates = static_cast<luint>(pow(2UL, bQbits));
}

string QASMExperiment::to_string() {
    stringstream output;
    output << "QASM Experiment (" << this->name << ") from file " << this->path << " with " << this->size() << " q-bits.";

    return output.str();
}