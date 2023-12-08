#include "Experiment.hpp"

using namespace std;

/**
 * Class for clauses in a formula. 
 * 
 * A clause is a a formula with the following shape: "x1 v x2 v x3...", where `vi` is a variable (possibly negated).
 * The variables can only appear once (since appearing together with different value implies the clause is always true).
*/
class Clause {
    private:
        unordered_map<luint, bool> elements_in;

    public:
        Clause(vector<pair<luint,bool>>);
        Clause(string);

        /* Method that creates a new Clause with random elements */
        static Clause* random(luint,luint);

        /* Method to evaluate the clause */
        bool is_trivial() { return this->elements_in.size() == 0; }
        /* Method to evaluate the clause */
        bool eval(vector<bool>);

        /* Return a list of variable sin the clause */
        vector<luint> variables();
        /* Return a list of variables that are not negated */
        vector<luint> pos_variables();
        /* Return a list of variables negated */
        vector<luint> neg_variables();

        /* Check equality over clauses */
        bool operator==(const Clause&);

        /* Transforms a clause into a string */
        string to_string();
};

/**
 * Class for a SAT formula.
 * 
 * A SAT formula is a list of non-trivial Clauses. It also knows how many variables the formula depend on
*/
class SATFormula : public Experiment {
    private:
        vector<Clause*> clauses;
        luint max_variables;

    public:
        SATFormula(luint nvars, luint eIterations, ExperimentType eType) : Experiment("SAT", "H", eIterations, eType) {this->max_variables = nvars; }
        SATFormula(string, luint, ExperimentType);
        ~SATFormula();

        luint add_clause(Clause*);
        static SATFormula* random(luint, luint, bool = true, luint = 3, ExperimentType = ExperimentType::DIRECT);
        string to_string();

        /* Virtual methods from Experiment */
        luint size();
        luint correct_size();
        array<dd::CMat, 2U> direct();
        vector<CCSparseVector> matrix();
        dd::CMat matrix_B(dd::CMat&);
        qc::QuantumComputation quantum();
        qc::QuantumComputation quantum_B();
        SATFormula& change_exec_type(ExperimentType);
};