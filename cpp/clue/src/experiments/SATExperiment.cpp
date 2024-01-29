#include "experiments/SATExperiment.hpp"

#include "operations/Operation.hpp"
#include <cstdlib>
#include <string>

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
Clause::Clause(string clause) {
    if (clause[0] != '(' || clause[clause.length()-1] != ')') {
        throw logic_error("Malformed clause: it must begin with '(' and end with ')'");
    }

    luint start = 1, next_or = clause.find('v');
    while (start < clause.length()) {
        if (next_or == string::npos) { next_or = clause.length()-1; }
        string variable = clause.substr(start, next_or-start);
        luint var_number = stoul(variable.substr(variable.find("_")));
        bool var_value = variable.find("-") == string::npos;
        if (this->elements_in.contains(var_number)) {
            if (this->elements_in[var_number] != var_value) {
                this->elements_in.clear();
                return;
            }
        } else {
            this->elements_in[var_number] = var_value;
        }
        start = next_or+1; next_or = clause.find('v', start);
    }
}
/*static*/ Clause* Clause::random(luint max_variables, luint num_variables) {
    vector<pair<luint,bool>> variables;
    for (luint i = 0; i < num_variables; i++) {
        variables.push_back({static_cast<luint>(rand()) % max_variables, static_cast<bool>(rand() % 2)});
    }
    return new Clause(variables);
}
/*static*/ Clause* Clause::copy(Clause* clause) {
    vector<pair<luint,bool>> to_add = vector<pair<luint,bool>>(clause->variables().size());
    for (luint pos : clause->pos_variables()) {
        to_add.push_back({pos, true});
    }
    for (luint neg : clause->neg_variables()) {
        to_add.push_back({neg, false});
    }
    Clause * result = new Clause(to_add);
    return result;
}
/* Method to evaluate the clause */
bool Clause::eval(boost::dynamic_bitset<> values) {
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

/*** CODE FOR CLASS SATFORMULA ***/
SATFormula::SATFormula(string formula, luint eIterations, ExperimentType eType, dd::Package<>* ePackage) : Experiment("SAT", "H", eIterations, eType, ePackage) {
    if (formula[0] != '[') {
        throw logic_error("A formula must start with an opening bracket '['");
    }
    luint beg = 1, i = 1;
    bool opened = false;
    while(formula[i] != ']') {
        if (!opened && formula[i] == ')') { 
            beg = i; opened = true;
        } else if (formula[i] == ')') {
            throw range_error("Expression is not well formed: found two opening parenthesis together.");
        } else if (opened && formula[i] == ')') {
            this->add_clause(new Clause(formula.substr(beg, i-beg+1)));
            opened = false;
        } else if (formula[i] == ')') {
            throw range_error("Expression is not well formed: found closing parenthesis without opener.");
        }
        i++;
    }
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
            for (luint var : clause->variables()) {
                if (var >= this->max_variables) { 
                    throw logic_error("The clause exceed the variables valid for this formula"); 
                }
            }
        }
        if (!found) { this->clauses.push_back(to_add); }
    }

    return this->clauses.size();
}

/*static*/ SATFormula* SATFormula::random(luint max_variables, luint max_clauses, bool force_clauses, luint max_per_clause, ExperimentType type, dd::Package<>* ePackage) {
    luint iterations = static_cast<luint>(ceil(pow(2., static_cast<double>(max_variables)/2.)));
    SATFormula * result = new SATFormula(max_variables, iterations, type, ePackage);
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
bool SATFormula::eval(boost::dynamic_bitset<> values) {
    if (values.size() < this->max_variables) {
        throw logic_error("Insufficient number of values provided.");
    }
    
    for (Clause* clause : this->clauses) {
        if (! clause->eval(values)) { return false; }
    }
    return true;
}
luint SATFormula::count(boost::dynamic_bitset<> values) {
    if (values.size() < this->max_variables) {
        throw logic_error("Insufficient number of values provided.");
    }
    
    luint value = 0;
    for (Clause* clause : this->clauses) {
        if (clause->eval(values)) { value++; }
    }
    return value;
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

void SATFormula::compute_possible_values() {
    if (this->possible_values.size() == 0) { // Only do something if this was not called before
        for (luint i = 0; i < static_cast<luint>(pow(2, this->max_variables)); i++) {
            luint new_value = this->count(boost::dynamic_bitset<>(this->max_variables, i));
            if (!this->possible_values.contains(new_value)) {
                this->possible_values[new_value] = vector<luint>();
            }
            this->possible_values[new_value].push_back(i);
        }
    }
    return;
}

/* Virtual methods from Experiment */
luint SATFormula::size() {
    return this->max_variables;
}
luint SATFormula::correct_size() {
    this->compute_possible_values();
    return this->possible_values.size();
}
luint SATFormula::bound_size() {
    return this->clauses.size();
}
array<dd::CMat, 2U> SATFormula::direct() {
    luint d = this->correct_size(); // This computes possible_values
    luint full_size = static_cast<luint>(pow(2, this->max_variables));
    dd::CMat L = dd::CMat(d), U = dd::CMat(d); 
    luint i = 0;
    for (std::pair<luint,vector<luint>> pair : this->possible_values) {
        L[i] = dd::CVec(full_size);
        CC value = CC(1./sqrt(pair.second.size()));
        for (luint j : pair.second) { L[i][j] = value; }
        U[i] = dd::CVec(d);
        U[i][i] = exp(CC(0, -static_cast<double>(pair.first)));
        i++;
    }
    return {L,U}; // TODO
}
vector<CCSparseVector> SATFormula::matrix() {
    luint full_size = static_cast<luint>(pow(2, this->max_variables));
    vector<CCSparseVector> result = vector<CCSparseVector>(full_size, full_size);
    for (luint i = 0; i < full_size; i++) {
        luint value = this->count(boost::dynamic_bitset(this->max_variables, i));
        result[i].set_value(i, CC(static_cast<double>(value)));
    }
    return result;
}
dd::CMat SATFormula::matrix_B(dd::CMat& Uhat) {
    if (is_diagonal(Uhat)) {
        // Not yet implemented: we return an identity
        dd::CMat result = dd::CMat(Uhat.size());
        for (luint i = 0; i < Uhat.size(); i++) {
            result[i] = dd::CVec(Uhat.size());
            result[i][i] = CC(1.);
        }
        return result;
        // cerr << "Found a diagonal matrix: \n" << matrix_to_string(Uhat) << endl;
        // throw logic_error("Begin matrix not implemented for diagonal matrices");
    }
    luint num_clauses = this->clauses.size();
    dd::CMat result = dd::CMat(Uhat.size());
    for (luint i = 0; i < Uhat.size(); i++) {
        result[i] = dd::CVec(Uhat[i].size());
        result[i][i] = CC(1./static_cast<double>(num_clauses));
    }
    result[0][0] = CC(static_cast<double>(num_clauses));

    return result;
}
qc::QuantumComputation* SATFormula::quantum(double par_val) {
    qc::QuantumComputation* circuit = new qc::QuantumComputation(this->max_variables); 

    for (Clause* clause : this->clauses) {
        vector<qc::Qubit> pos, neg; // Convert the positive and negative variables into the Qubits
        for (luint p : clause->pos_variables()) { pos.push_back(static_cast<qc::Qubit>(p)); }
        for (luint n : clause->neg_variables()) { neg.push_back(static_cast<qc::Qubit>(n)); }

        if (pos.size() == 3) { // TTT case 
            qc::Controls controls{};
            controls.emplace(pos[0]); controls.emplace(pos[1]);

            circuit->x(pos[0]);
            circuit->cx(pos[0], pos[1]);  
            circuit->mcx(controls, pos[2]);
            circuit->mcp(-par_val, controls, pos[2]);
            circuit->mcx(controls, pos[2]);
            circuit->cx(pos[0], pos[1]);  
            circuit->x(pos[0]);
        } else if (pos.size() == 2 && neg.size() == 1) { // FTT case
            qc::Controls controls{};
            controls.emplace(neg[0]); controls.emplace(pos[0]);

            circuit->cx(neg[0],pos[0]);
            circuit->mcx(controls,pos[1]);
            circuit->mcp(-par_val, controls,pos[1]);
            circuit->mcx(controls,pos[1]);
            circuit->cx(neg[0],pos[0]);
        } else if (pos.size() == 2 && neg.size() == 0) { // TT case
            circuit->x(pos[0]);
            circuit->cx(pos[0], pos[1]);
            circuit->cp(-par_val, pos[0], pos[1]);
            circuit->cx(pos[0], pos[1]);
            circuit->x(pos[0]);
        } else if (pos.size() == 1 && neg.size() == 2) { // FFT case
            qc::Controls controls{};
            controls.emplace(neg[0]); controls.emplace(neg[1]);

            circuit->mcx(controls, pos[0]);
            circuit->mcp(-par_val, controls, pos[0]);
            circuit->mcx(controls, pos[0]);
        } else if (pos.size() == 1 && neg.size() == 1) { // FT case
            circuit->cx(neg[0], pos[0]);
            circuit->cp(-par_val, neg[0], pos[0]);
            circuit->cx(neg[0], pos[0]);
        } else if (pos.size() == 1 && neg.size() == 0) { // T case
            circuit->x(pos[0]);
            circuit->p(-par_val, pos[0]);
            circuit->x(pos[0]);
        } else if (pos.size() == 0 && neg.size() == 3) { // FFF case
            qc::Controls controls{};
            controls.emplace(neg[0]); controls.emplace(neg[1]);

            circuit->mcp(-par_val, controls, neg[2]);
        } else if (pos.size() == 0 && neg.size() == 2) { // FF case
            circuit->cp(-par_val, neg[0], neg[1]);
        } else if (pos.size() == 0 && neg.size() == 1) { // F case
            circuit->p(-par_val, neg[0]);
        } else {
            throw logic_error("Error in a clause: unrecognized set up of variables");
        }
    }
    return circuit;
}
qc::QuantumComputation* SATFormula::quantum_B(double par_val) {
    // Create a circuit with the appropriate number of qbits
    qc::QuantumComputation* circuit = new qc::QuantumComputation(this->max_variables); 

    // We create the new gates
    for (luint i = 0; i < this->max_variables; i++) {
        luint count = 0;
        for (Clause* clause : this->clauses) {
            for (luint var : clause->variables()) {
                if (var == i) { count++; break; }
            }
        }

        circuit->h(static_cast<qc::Qubit>(i));
        circuit->x(static_cast<qc::Qubit>(i));
        circuit->p(-par_val*static_cast<double>(count), static_cast<qc::Qubit>(i));
        circuit->x(static_cast<qc::Qubit>(i));
        circuit->h(static_cast<qc::Qubit>(i));
    }
    
    return circuit;
}

SATFormula* SATFormula::change_exec_type(ExperimentType new_type) {
    SATFormula* result = new SATFormula(this->max_variables, this->iterations, new_type, this->package);
    for (Clause* clause : this->clauses) {
        result->add_clause(Clause::copy(clause));
    }
    return result;
}