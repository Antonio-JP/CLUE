#ifndef CLUE_EX_BENCHMARK
#define CLUE_EX_BENCHMARK

#include "Experiment.hpp"
#include "boost/dynamic_bitset.hpp"

using namespace std;

/**
 * Class for clauses in a formula. 
 * 
 * A clause is a a formula with the following shape: "x1 v x2 v x3...", where `vi` is a variable (possibly negated).
 * The variables can only appear once (since appearing together with different value implies the clause is always true).
*/
class QASMExperiment : public Experiment{
    private:
        luint nstates;
        qc::QuantumComputation* circuit;

    protected:
        luint qbits;
        string path;

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
        QASMExperiment* change_exec_type(ExperimentType);
    public:
        QASMExperiment(luint, string, string, string, ExperimentType, dd::Package<>*);

        /* Method to get the string out of an experiment */
        string to_string();
};

#endif