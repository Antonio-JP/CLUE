#ifndef CLUE_EX_BENCHMARK
#define CLUE_EX_BENCHMARK

#include "Experiment.hpp"
#include "boost/dynamic_bitset.hpp"

using namespace std;

enum CircuitType { // 5 per row
    AE, DJ, GHZ, GRAPHSTATE, HHL, 
    PRICINGPUT, PRICINGCALL, PORTFOLIOQAOA, PORTFOLIOVQE, QFT, 
    QPEEXACT, QPEINEXACT, QWALK, TSP, QNN,
    VQE, WSTATE
};
CircuitType CircuitType_fromString(string);
string CircuitType_toString(CircuitType);
string CircuitType_toName(CircuitType);



/**
 * Class for clauses in a formula. 
 * 
 * A clause is a a formula with the following shape: "x1 v x2 v x3...", where `vi` is a variable (possibly negated).
 * The variables can only appear once (since appearing together with different value implies the clause is always true).
*/
class BenchmarkExperiment : public Experiment{
    private:
        luint nstates;
        qc::QuantumComputation* circuit;
        std::unique_ptr<dd::Package<>> package;

    protected:
        luint qbits;

        /* Virtual methods from Experiment */
        luint size() {return this->qbits; }
        luint correct_size();
        luint bound_size();
        array<dd::CMat, 2U> direct();
        vector<CCSparseVector> matrix();
        dd::CMat matrix_B(dd::CMat&);
        qc::QuantumComputation* quantum() { return this->quantum(0.); }
        qc::QuantumComputation* quantum(double);
        qc::QuantumComputation* quantum_B() { return this->quantum_B(0.); }
        qc::QuantumComputation* quantum_B(double);
        BenchmarkExperiment* change_exec_type(ExperimentType);
    public:
        BenchmarkExperiment(luint, string, string, ExperimentType, std::unique_ptr<dd::Package<>>&);

        /* Method to get the string out of an experiment */
        string to_string();
};

#endif