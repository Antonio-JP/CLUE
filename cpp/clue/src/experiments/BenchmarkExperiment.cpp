#include "experiments/BenchmarkExperiment.hpp"

#include <boost/algorithm/string.hpp>
#include "dd/Operations.hpp"

/*** CODE FOR CIRCUIT_TYPE ENUM ***/
CircuitType CircuitType_fromString(string name) {
    string name_upper = boost::to_upper_copy<std::string>(name);
    if (name_upper == "AE") {
        return CircuitType::AE;
    } else if (name_upper == "DJ") {
        return CircuitType::DJ;
    } else if (name_upper == "GHZ") {
        return CircuitType::GHZ;
    } else if (name_upper == "GRAPHSTATE") {
        return CircuitType::GRAPHSTATE;
    } else if (name_upper == "HHL") {
        return CircuitType::HHL;
    } else if (name_upper == "PORTFOLIOQAOA") {
        return CircuitType::PORTFOLIOQAOA;
    } else if (name_upper == "PORTFOLIOVQE") {
        return CircuitType::PORTFOLIOVQE;
    } else if (name_upper == "PRICINGCALL") {
        return CircuitType::PRICINGCALL;
    } else if (name_upper == "PRICINGPUT") {
        return CircuitType::PRICINGPUT;
    } else if (name_upper == "QFT") {
        return CircuitType::QFT;
    } else if (name_upper == "QNN") {
        return CircuitType::QNN;
    } else if (name_upper == "QPEEXACT") {
        return CircuitType::QPEEXACT;
    } else if (name_upper == "QPEINEXACT") {
        return CircuitType::QPEINEXACT;
    } else if (name_upper == "QWALK") {
        return CircuitType::QWALK;
    } else if (name_upper == "TSP") {
        return CircuitType::TSP;
    } else if (name_upper == "VQE") {
        return CircuitType::VQE;
    } else if (name_upper == "WSTATE") {
        return CircuitType::WSTATE;
    } else {
        throw domain_error("Circuit type " + name_upper + " not recognized.");
    }
}

string CircuitType_toString(CircuitType type) {
    switch (type)
    {
        case CircuitType::AE:
            return "AE";
        case CircuitType::DJ:
            return "DJ";
        case CircuitType::GHZ:
            return "GHZ";
        case CircuitType::GRAPHSTATE:
            return "GRAPHSTATE";
        case CircuitType::HHL:
            return "HHL";
        case CircuitType::PORTFOLIOQAOA:
            return "PORTFOLIOQAOA";
        case CircuitType::PORTFOLIOVQE:
            return "PORTFOLIOVQE";
        case CircuitType::PRICINGCALL:
            return "PRICINGCALL";
        case CircuitType::PRICINGPUT:
            return "PRICINGPUT";
        case CircuitType::QFT:
            return "QFT";
        case CircuitType::QNN:
            return "QNN";
        case CircuitType::QPEEXACT:
            return "QPEEXACT";
        case CircuitType::QPEINEXACT:
            return "QPEINEXACT";
        case CircuitType::QWALK:
            return "QWALK";
        case CircuitType::TSP:
            return "TSP";
        case CircuitType::VQE:
            return "VQE";
        case CircuitType::WSTATE:
            return "WSTATE";    
        default:
            throw domain_error("Circuit type not recognized.");
    }
}

string CircuitType_toName(CircuitType type) {
    switch (type)
    {
        case CircuitType::AE:
            return "Amplitude Estimation";
        case CircuitType::DJ:
            return "Deutsch-Jozsa";
        case CircuitType::GHZ:
            return "Greenberger–Horne–Zeilinger ";
        case CircuitType::GRAPHSTATE:
            return "Graph State";
        case CircuitType::HHL:
            return "Harrow–Hassidim–Lloyd Algorithm";
        case CircuitType::PORTFOLIOQAOA:
            return "Portfolio QAOA";
        case CircuitType::PORTFOLIOVQE:
            return "Portfolio QVE";
        case CircuitType::PRICINGCALL:
            return "Pricing-Call";
        case CircuitType::PRICINGPUT:
            return "Pricing-Put";
        case CircuitType::QFT:
            return "Quantum Fourier Transform";
        case CircuitType::QNN:
            return "Quantum Neural Network";
        case CircuitType::QPEEXACT:
            return "Quantum Phase Estimation (exact)";
        case CircuitType::QPEINEXACT:
            return "Quantum Phase Estimation (inexact)";
        case CircuitType::QWALK:
            return "Quantum Walk";
        case CircuitType::TSP:
            return "Travelling-Salesman Problem";
        case CircuitType::VQE:
            return "Variational Quantum Eigensolver";
        case CircuitType::WSTATE:
            return "W-State";    
        default:
            throw domain_error("Circuit type not recognized.");
    }
}

// Methods to create the circuits using the description in the Python implementation of mqt.bench
/** AE **/
qc::QuantumComputation* create_ae(luint) {
    return nullptr;
}
/** DJ **/
void dj_oracle(qc::QuantumComputation& circuit) {
    luint qbits = circuit.getNqubits();
    circuit.x(static_cast<qc::Qubit>(qbits-1));
}
void dj_algorithm(qc::QuantumComputation& circuit) {
    qc::Qubit n = static_cast<qc::Qubit>(circuit.getNqubits());

    circuit.x(n-1);
    circuit.h(n-1);

    for (qc::Qubit i = 0; i < n-1; i++) {
        circuit.h(i);
    }

    dj_oracle(circuit);

    for (qc::Qubit i = 0; i < n-1; i++) {
        circuit.h(i);
    }
}
qc::QuantumComputation* create_dj(luint qubits) {
    qc::QuantumComputation* result = new qc::QuantumComputation(qubits);
    dj_algorithm(*result);
    return result;
}

/** GHZ **/
qc::QuantumComputation* create_ghz(luint) {
    return nullptr;
}
/** GRAPHSTATE **/
qc::QuantumComputation* create_graphstate(luint) {
    return nullptr;
}
/** HHL **/
qc::QuantumComputation* create_hhl(luint) {
    return nullptr;
}
/** PRICINGPUT **/
qc::QuantumComputation* create_pricingput(luint) {
    return nullptr;
}
/** PRICINGCALL **/
qc::QuantumComputation* create_pricingcall(luint) {
    return nullptr;
}
/** PORTFOLIOQAOA **/
qc::QuantumComputation* create_portfolioqaoa(luint) {
    return nullptr;
}
/** PORTFOLIOVQE **/
qc::QuantumComputation* create_portfoliovqe(luint) {
    return nullptr;
}
/** QFT **/
qc::QuantumComputation* create_qft(luint) {
    return nullptr;
}
/** QPEEXACT **/
qc::QuantumComputation* create_qpeexact(luint) {
    return nullptr;
}
/** QPEINEXACT **/
qc::QuantumComputation* create_qpeinexact(luint) {
    return nullptr;
}
/** QWALK **/
qc::QuantumComputation* create_qwalk(luint) {
    return nullptr;
}
/** TSP **/
qc::QuantumComputation* create_tsp(luint) {
    return nullptr;
}
/** QNN **/
qc::QuantumComputation* create_qnn(luint) {
    return nullptr;
}
/** VQE **/
qc::QuantumComputation* create_vqe(luint) {
    return nullptr;
}
/** WSTATE **/
qc::QuantumComputation* create_wstate(luint) {
    return nullptr;
}
/*** CODE FOR CLASS BENCHMARK_EXAMPLE ***/
luint BenchmarkExperiment::correct_size() {
    return -1UL; // Unknown
}
luint BenchmarkExperiment::bound_size() {
    return this->nstates; // No preknown bound
}

array<dd::CMat, 2U> BenchmarkExperiment::direct() {
    throw logic_error("Direct type not implemented for Benchmark examples");
}

vector<CCSparseVector> BenchmarkExperiment::matrix() {
    throw logic_error("Matrix from circuit not yet implemented");
}
dd::CMat BenchmarkExperiment::matrix_B(dd::CMat& U) {
    return identity_matrix(U.size()); // There is no begin hamiltonian: we use the identity
}

qc::QuantumComputation* BenchmarkExperiment::quantum(double) {
    if (this->circuit == nullptr) {
        CircuitType circuit_type = CircuitType_fromString(this->name);
        switch (circuit_type)
        {
            case CircuitType::AE:
                this->circuit = create_ae(this->size());
                break;
            case CircuitType::DJ:
                this->circuit = create_dj(this->size());
                break;
            case CircuitType::GHZ:
                this->circuit = create_ghz(this->size());
                break;
            case CircuitType::GRAPHSTATE:
                this->circuit = create_graphstate(this->size());
                break;
            case CircuitType::HHL:
                this->circuit = create_hhl(this->size());
                break;
            case CircuitType::PORTFOLIOQAOA:
                this->circuit = create_portfolioqaoa(this->size());
                break;
            case CircuitType::PORTFOLIOVQE:
                this->circuit = create_portfoliovqe(this->size());
                break;
            case CircuitType::PRICINGCALL:
                this->circuit = create_pricingcall(this->size());
                break;
            case CircuitType::PRICINGPUT:
                this->circuit = create_pricingput(this->size());
                break;
            case CircuitType::QFT:
                this->circuit = create_qft(this->size());
                break;
            case CircuitType::QNN:
                this->circuit = create_qnn(this->size());
                break;
            case CircuitType::QPEEXACT:
                this->circuit = create_qpeexact(this->size());
                break;
            case CircuitType::QPEINEXACT:
                this->circuit = create_qpeinexact(this->size());
                break;
            case CircuitType::QWALK:
                this->circuit = create_qwalk(this->size());
                break;
            case CircuitType::TSP:
                this->circuit = create_tsp(this->size());
                break;
            case CircuitType::VQE:
                this->circuit = create_vqe(this->size());
                break;
            case CircuitType::WSTATE:
                this->circuit = create_wstate(this->size());
                break;
            default:
                throw domain_error("Circuit type not recognized.");
        }
    }

    if (this->circuit == nullptr) {
        throw logic_error("The case " + this->name + " is not yet implemented.");
    }
    return this->circuit;
}

qc::QuantumComputation* BenchmarkExperiment::quantum_B(double) { // A nothing circuit
    qc::QuantumComputation* empty = new qc::QuantumComputation(this->size());

    return empty;
}

BenchmarkExperiment* BenchmarkExperiment::change_exec_type(ExperimentType new_type) {
    return new BenchmarkExperiment(this->size(), this->name, this->observable, new_type, this->package);
}

/*** BUILDERS FOR BENCHMARK_EXAMPLE ***/
BenchmarkExperiment::BenchmarkExperiment(luint bQbits, string eName, string eObservable, ExperimentType eType, std::unique_ptr<dd::Package<>>& bPackage) : 
    Experiment(eName, eObservable, 1UL, eType) {
    this->qbits = bQbits;
    this->package = std::move(bPackage);
    this->circuit = nullptr;
    this->nstates = static_cast<luint>(pow(2UL, bQbits));
}

string BenchmarkExperiment::to_string() {
    stringstream output;
    output << CircuitType_toName(CircuitType_fromString(this->name)) << " with " << this->size() << " q-bits.";

    return output.str();
}