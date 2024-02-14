#include "experiments/BenchmarkExperiment.hpp"

#include <boost/algorithm/string.hpp>
#include <cmath>
#include "CircuitOptimizer.hpp"
#include "dd/Operations.hpp"
#include "dd/FunctionalityConstruction.hpp"
#include "experiments/CUTExperiment.hpp"
#include "experiments/benchmark/AmplitudeEstimation.hpp"
#include "experiments/benchmark/DeutschJozsa.hpp"
#include "experiments/benchmark/GHZ.hpp"
#include "experiments/benchmark/GraphState.hpp"
#include "experiments/benchmark/HHL.hpp"
#include "experiments/benchmark/PortfolioQAOA.hpp"
#include "experiments/benchmark/PortfolioVQE.hpp"
#include "experiments/benchmark/PricingCall.hpp"
#include "experiments/benchmark/PricingPut.hpp"
#include "experiments/benchmark/QPEExact.hpp"
#include "experiments/benchmark/QPEInexact.hpp"
#include "experiments/benchmark/QuantumFourierTransform.hpp"
#include "experiments/benchmark/QuantumNeuralNetwork.hpp"
#include "experiments/benchmark/QuantumWalk.hpp"
#include "experiments/benchmark/TravellingSalesman.hpp"
#include "experiments/benchmark/VariationalQuantumEigensolver.hpp"
#include "experiments/benchmark/WState.hpp"

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
/*** CODE FOR CLASS BENCHMARK_EXAMPLE ***/
luint BenchmarkExperiment::correct_size() {
    return 0; // Unknown
}
luint BenchmarkExperiment::bound_size() {
    return this->nstates; // No preknown bound
}

array<dd::CMat, 2U> BenchmarkExperiment::direct() {
    throw logic_error("Direct type not implemented for Benchmark examples");
}

vector<CCSparseVector> BenchmarkExperiment::matrix() {
    this->quantum(0.0);
    dd::mEdge circuit_dd = buildFunctionality(this->circuit, *this->package);
    dd::CMat unitary = circuit_dd.getMatrix(); 
    luint N = unitary.size();
    vector<CCSparseVector> U = vector<CCSparseVector>(N, N); // Initializing the final matrix
    for (luint i = 0; i < N; i++) {
        for (luint j = 0; j < N; j++) {
            U[i].set_value(j, unitary[i][j]);
        }
    }

    return U;
}
dd::CMat BenchmarkExperiment::matrix_B(dd::CMat& U) {
    return identity_matrix(U.size()); // There is no begin hamiltonian: we use the identity
}

qc::QuantumComputation* BenchmarkExperiment::quantum(double) {
    if (this->circuit == nullptr) {
        CircuitType circuit_type = CircuitType_fromString(this->name);
        string file_name = boost::to_lower_copy(CircuitType_toString(circuit_type)) + "_indep_qiskit_" + std::to_string(this->qbits) + ".qasm";
        this->circuit = new qc::QuantumComputation("../../../../tests/quantum/circuits/" + file_name);
        qc::CircuitOptimizer::removeFinalMeasurements(*this->circuit);
        // switch (circuit_type)
        // {
        //     case CircuitType::AE:
        //         this->circuit = ae::create(this->size());
        //         break;
        //     case CircuitType::DJ:
        //         this->circuit = dj::create(this->size());
        //         break;
        //     case CircuitType::GHZ:
        //         this->circuit = ghz::create(this->size());
        //         break;
        //     case CircuitType::GRAPHSTATE:
        //         this->circuit = graphstate::create(this->size());
        //         break;
        //     case CircuitType::HHL:
        //         this->circuit = hhl::create(this->size());
        //         break;
        //     case CircuitType::PORTFOLIOQAOA:
        //         this->circuit = portfolioqaoa::create(this->size());
        //         break;
        //     case CircuitType::PORTFOLIOVQE:
        //         this->circuit = portfoliovqe::create(this->size());
        //         break;
        //     case CircuitType::PRICINGCALL:
        //         this->circuit = pricingcall::create(this->size());
        //         break;
        //     case CircuitType::PRICINGPUT:
        //         this->circuit = pricingput::create(this->size());
        //         break;
        //     case CircuitType::QFT:
        //         this->circuit = qft::create(this->size());
        //         break;
        //     case CircuitType::QNN:
        //         this->circuit = qnn::create(this->size());
        //         break;
        //     case CircuitType::QPEEXACT:
        //         this->circuit = qpeexact::create(this->size());
        //         break;
        //     case CircuitType::QPEINEXACT:
        //         this->circuit = qpeinexact::create(this->size());
        //         break;
        //     case CircuitType::QWALK:
        //         this->circuit = qwalk::create(this->size());
        //         break;
        //     case CircuitType::TSP:
        //         this->circuit = tsp::create(this->size());
        //         break;
        //     case CircuitType::VQE:
        //         this->circuit = vqe::create(this->size());
        //         break;
        //     case CircuitType::WSTATE:
        //         this->circuit = wstate::create(this->size());
        //         break;
        //     default:
        //         throw domain_error("Circuit type not recognized.");
        // }
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
BenchmarkExperiment::BenchmarkExperiment(luint bQbits, string eName, string eObservable, ExperimentType eType, dd::Package<>* ePackage) : 
    Experiment(eName, eObservable, 1UL, eType, ePackage) {
    this->qbits = bQbits;
    this->circuit = nullptr;
    this->nstates = static_cast<luint>(pow(2UL, bQbits));
}

string BenchmarkExperiment::to_string() {
    stringstream output;
    output << CircuitType_toName(CircuitType_fromString(this->name)) << " with " << this->size() << " q-bits.";

    return output.str();
}