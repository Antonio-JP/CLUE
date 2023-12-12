#ifndef CLUE_EX_SEARCH
#define CLUE_EX_SEARCH

#include "Experiment.hpp"
#include "boost/dynamic_bitset.hpp"

using namespace std;

/**
 * Class for clauses in a formula. 
 * 
 * A clause is a a formula with the following shape: "x1 v x2 v x3...", where `vi` is a variable (possibly negated).
 * The variables can only appear once (since appearing together with different value implies the clause is always true).
*/
class QuantumSearch : public Experiment{
    protected:
        luint qbits;
        unordered_set<luint> success_set;

        /* Method that serves as an oracle for the search function */
        bool oracle(boost::dynamic_bitset<>);
        bool oracle(luint);
        void quantum_oracle(qc::QuantumComputation&); // Adds the oracle part to the circuit given as input
        void quantum_diffusion(qc::QuantumComputation&); // Adds the diffusion operator to the circuit given as input

        /* Overriden methods from Experiment */
        CCSparseVector clue_observable();
        dd::vEdge dd_observable(std::unique_ptr<dd::Package<>>&);
        /* Virtual methods from Experiment */
        luint size() {return this->qbits; }
        luint correct_size() { return 2UL; }
        luint bound_size() { return 2UL; }
        array<dd::CMat, 2U> direct();
        vector<CCSparseVector> matrix();
        dd::CMat matrix_B(dd::CMat&);
        qc::QuantumComputation* quantum(double);
        qc::QuantumComputation* quantum_B(double);
        QuantumSearch* change_exec_type(ExperimentType);
    public:
        QuantumSearch(luint, vector<luint>, luint, ExperimentType);

        static QuantumSearch* random(luint, ExperimentType);

        /* Method to get the string out of an experiment */
        string to_string();
};

#endif