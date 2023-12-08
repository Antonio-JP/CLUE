#include "SATExperiment.hpp"

#include <cstdlib>

/*** CODE FOR CLASS CLAUSE ***/
Clause::Clause(vector<pair<luint,bool>> variables) {
    for (pair<luint, bool> var : variables) {
        if (this->elements_in.contains(var.first)) {
            if (this->elements_in[var.first] != var.second) {
                this->elements_in.clear();
                return;
            }
        } else {
            this->elements_in[var.first] = var.second;
        }
    }
}
Clause::Clause(string) {
    throw logic_error("Constructor from String not yet implemented");
}
/*static*/ Clause* Clause::random(luint max_variables, luint num_variables) {
    vector<pair<luint,bool>> variables;
    for (luint i = 0; i < num_variables; i++) {
        variables.push_back({rand() % max_variables, static_cast<bool>(rand() % 2)});
    }
    return new Clause(variables);
}
/* Method to evaluate the clause */
bool Clause::eval(vector<bool> values) {
    if (this->is_trivial()) {
        return true;
    }
    for (pair<luint, bool> var : this->elements_in) {
        if (values.size() < var.first) { throw logic_error("Not enough elements given to evaluate a clause "); }
        if (var.second == values[var.first]) { // We found a true variable --> all disjunction is True
            return true;
        }
    }
    return false; // Nothing found --> the clause is false
}

/* Return a list of variable sin the clause */
vector<luint> Clause::variables() {
    vector<luint> elements = vector<luint>(this->elements_in.size());
    luint i = 0;
    for (pair<luint, bool> var : this->elements_in) {
        elements[i] = var.first;
        i++;
    }
    return elements;
}
/* Return a list of variables that are not negated */
vector<luint> Clause::pos_variables() {
    vector<luint> elements = vector<luint>();
    for (pair<luint, bool> var : this->elements_in) {
        if (var.second) { elements.push_back(var.first); }
    }
    return elements;
}
/* Return a list of variables negated */
vector<luint> Clause::neg_variables() {
    vector<luint> elements = vector<luint>();
    for (pair<luint, bool> var : this->elements_in) {
        if (!var.second) { elements.push_back(var.first); }
    }
    return elements;
}

/* Check equality over clauses */
bool Clause::operator==(const Clause& other) {
    return this->elements_in == other.elements_in;
}

/* Transforms a clause into a string */
string Clause::to_string() {
    if (this->is_trivial()) { return "()"; }
    stringstream stream;
    stream << "(";
    std::unordered_map<luint,bool>::iterator it = this->elements_in.begin();
    if (!it->second) { stream << "-"; }
    stream << "x_" << it->first;

    it++;
    while (it != this->elements_in.end()) {
        stream << " v ";
        if (!it->second) { stream << "-"; }
        stream << "x_" << it->first;
        it++;
    }
    stream << ")";
    
    return stream.str();
}

SATFormula::SATFormula(string, luint eIterations, ExperimentType eType) : Experiment("SAT", "H", eIterations, eType) {
    throw logic_error("Constructor from string not yet defined");
}
SATFormula::~SATFormula() {
    for (Clause* clause : this->clauses) {
        delete clause;
    }
}

luint SATFormula::add_clause(Clause* to_add) {
    if (!to_add->is_trivial()) { // We avoid trivial clauses
        bool found = false;
        for (Clause* clause : this->clauses) { // We avoid repeated clauses
            if (*to_add == *clause) {
                found = true; break;
            }
        }
        if (!found) { this->clauses.push_back(to_add); }
    }
}

SATFormula* SATFormula::random(luint max_variables, luint max_clauses, bool force_clauses = true, luint max_per_clause = 3, ExperimentType type = ExperimentType::DIRECT) {
    luint iterations = static_cast<luint>(ceil(pow(2., max_variables/2.)));
    SATFormula * result = new SATFormula(max_variables, iterations, type);
    for (luint i = 0; i < max_clauses; i++) {
        result->add_clause(Clause::random(max_variables, max_per_clause));
    }

    if (force_clauses) {
        while (result->clauses.size() < max_clauses) {
            result->add_clause(Clause::random(max_variables, max_per_clause));
        }
    }

    return result;
}
string SATFormula::to_string() {
    stringstream stream;
    stream << "[";
    if (this->clauses.size()) {
        std::vector<Clause*>::iterator it = this->clauses.begin();
        stream << (*it)->to_string();
        it++;

        while (it != this->clauses.end()) {
            stream << " ^ " << (*it)->to_string();
            it++;
        }
    }
    stream << "]";
    return stream.str();
}

/* Virtual methods from Experiment */
luint SATFormula::size() {
    return this->max_variables;
}
luint SATFormula::correct_size() {
    return 0; // TODO
}
array<dd::CMat, 2U> SATFormula::direct() {
    return; // TODO
}
vector<CCSparseVector> SATFormula::matrix() {
    return; // TODO
}
dd::CMat SATFormula::matrix_B(dd::CMat&) {
    return; // TODO
}
qc::QuantumComputation SATFormula::quantum() {
    return; // TODO
}
qc::QuantumComputation SATFormula::quantum_B() {
    return; // TODO
}
SATFormula& SATFormula::change_exec_type(ExperimentType) {
    return; // TODO
}