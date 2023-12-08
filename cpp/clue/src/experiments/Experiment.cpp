#include "Experiment.hpp"

#include "dd/Simulation.hpp"

// Auxiliar method to convert time clocks into double time
double time_to_double(clock_t& init, clock_t& end) {
    return (double(end - init) / double(CLOCKS_PER_SEC));
}

// IMPLEMENTATION OF GENERIC METHODS OF CLASS Experiment

// PROTECTED METHODS
/* Method to get the observable for use with CLUE */
CCSparseVector Experiment::clue_observable() {
    CCSparseVector result = CCSparseVector(static_cast<luint>(pow(2, this->size())));
    if (this->observable == "H") {
        CC c = CC(1/sqrt(result.dimension()));
        for (luint i = 0; i < result.dimension(); i++) {
            result.set_value(i, c);
        }
    } else {
        result.set_value(stoi(this->observable), CC(1));
    }
    return result;
}
/* Method to get the observable for use with DD */
dd::vEdge Experiment::dd_observable() {
    std::unique_ptr<dd::Package<>> package(dd_package(this->size()));
    dd::vEdge init;
    if (this->observable == "H") {
        init = package->makeBasisState(this->size(), vector<dd::BasisStates>(this->size(), dd::BasisStates::plus), 0);
    } else {
        int val = stoi(this->observable);
        vector<bool> binary = vector<bool>(this->size(), false);
        for (luint i = this->size(); i > 0 && val > 0; i--) {
            binary[i-1] = (val % 2 == 1);
            val /= 2;
        }
        init = package->makeBasisState(this->size(), binary, 0);
    }

    return init;
}

// PRIVATE METHODS
/* Method that runs the CLUE reduction (only used when this->type == CLUE) */
void Experiment::run_clue() {
    cerr << "+++ [clue @ " << this->name << "] Computing CLUE reduction for " << this->name << endl;
    clock_t begin = clock();
    cerr << "+++ [clue @ " << this->name << "] Setting up observable and system..."<< endl;
    CCSparseVector observable = this->clue_observable();
    vector<CCSparseVector> U = this->matrix();
    luint dimension = static_cast<luint>(pow(2, this->size()));

    cerr << "+++ [clue @ " << this->name << "] Computing lumping...";
    clock_t b_lumping = clock();
    CCSubspace lumping = CCSubspace(dimension);
    vector<vector<CCSparseVector>> matrices = {U};
    lumping.absorb_new_vector(&observable);
    lumping.minimal_invariant_space(matrices);
    cerr << "+++ [clue @ " << this->name << "] Getting the reduced U_P...";
    dd::CMat Uhat = lumping.reduced_matrix(U);
    clock_t a_lumping = clock();

    cerr << "+++ [clue @ " << this->name << "] Getting the reduced U_B...";
    dd::CMat UB = this->matrix_B(Uhat);
    dd::CMat U_full = matmul(Uhat, UB);

    cerr << "+++ [clue @ " << this->name << "] Computing the iteration (U_P*U_B)^iterations...";
    clock_t b_iteration = clock();
    matrix_power(U_full, this->iterations);
    clock_t a_iteration = clock();

    clock_t end = clock();

    // We store the data
    this->red_ratio = lumping.dimension() / static_cast<double>(dimension);
    this->red_time = time_to_double(b_lumping, a_lumping);
    this->it_time = time_to_double(b_iteration, a_iteration);
    this->tot_time = time_to_double(begin, end);
 
    return;
}
/* Method that runs the DDSIM reduction (only used when this->type == DDSIM) */
void Experiment::run_ddsim() {
    cerr << "+++ [ddsim @ " << this->name << "] Computing DDSIM reduction for " << this->name << endl;
    clock_t begin = clock();
    cerr << "+++ [ddsim @ " << this->name << "] Setting up observable and system..."<< endl;
    dd::vEdge observable = this->dd_observable();
    qc::QuantumComputation U = this->quantum();
    luint dimension = static_cast<luint>(pow(2, this->size()));

    cerr << "+++ [ddsim @ " << this->name << "] Computing lumping...";
    clock_t b_lumping = clock();
    FullDDSubspace lumping = FullDDSubspace(dimension);
    vector<qc::QuantumComputation> circuits = {U};
    lumping.absorb_new_vector(&observable);
    lumping.minimal_invariant_space(circuits);
    cerr << "+++ [ddsim @ " << this->name << "] Getting the reduced U_P...";
    vector<vector<dd::ComplexValue>> or_Uhat = lumping.reduced_matrix(U);
    dd::CMat Uhat = ComplexValue_to_complex(or_Uhat);
    clock_t a_lumping = clock();

    cerr << "+++ [ddsim @ " << this->name << "] Getting the reduced U_B...";
    dd::CMat UB = this->matrix_B(Uhat);
    dd::CMat U_full = matmul(Uhat, UB);

    cerr << "+++ [ddsim @ " << this->name << "] Computing the iteration (U_P*U_B)^iterations...";
    clock_t b_iteration = clock();
    matrix_power(U_full, this->iterations);
    clock_t a_iteration = clock();

    clock_t end = clock();

    // We store the data
    this->red_ratio = lumping.dimension() / static_cast<double>(dimension);
    this->red_time = time_to_double(b_lumping, a_lumping);
    this->it_time = time_to_double(b_iteration, a_iteration);
    this->tot_time = time_to_double(begin, end);
 
    return;
}
/* Method that runs the DIRECT reduction (only used when this->type == DIRECT) */
void Experiment::run_direct() {
    cerr << "+++ [direct @ " << this->name << "] Computing direct reduction for " << this->name << endl;
    clock_t begin = clock();
    luint dimension = static_cast<luint>(pow(2, this->size()));

    cerr << "+++ [direct @ " << this->name << "] Computing lumping...";
    clock_t b_lumping = clock();
    array<dd::CMat, 2U> reduction = this->direct();
    dd::CMat& lumping = reduction[0];
    cerr << "+++ [direct @ " << this->name << "] Getting the reduced U_P...";
    dd::CMat& Uhat = reduction[1];
    clock_t a_lumping = clock();

    cerr << "+++ [direct @ " << this->name << "] Getting the reduced U_B...";
    dd::CMat UB = this->matrix_B(Uhat);
    dd::CMat U_full = matmul(Uhat, UB);

    cerr << "+++ [direct @ " << this->name << "] Computing the iteration (U_P*U_B)^iterations...";
    clock_t b_iteration = clock();
    matrix_power(U_full, this->iterations);
    clock_t a_iteration = clock();

    clock_t end = clock();

    // We store the data
    this->red_ratio = lumping.size() / static_cast<double>(dimension);
    this->red_time = time_to_double(b_lumping, a_lumping);
    this->it_time = time_to_double(b_iteration, a_iteration);
    this->tot_time = time_to_double(begin, end);
 
    return;
}
/* Method that runs the CLUE reduction (only used when this->type == DDSIM_ALONE) */
void Experiment::run_ddsim_alone() {
    cerr << "+++ [ddsim-only @ " << this->name << "] Computing DDSIM ONLY execution for " << this->name << endl;
    clock_t begin = clock();
    cerr << "+++ [ddsim-only @ " << this->name << "] Setting up observable and system..."<< endl;
    dd::vEdge observable = this->dd_observable();
    qc::QuantumComputation U_P = this->quantum();
    qc::QuantumComputation U_B = this->quantum_B();



    cerr << "+++ [ddsim-only @ " << this->name << "] Computing the iteration (U_P*U_B)^iterations...";
    clock_t b_iteration = clock();
    std::unique_ptr<dd::Package<>> package(dd_package(this->size()));
    dd::vEdge current = observable;
    for (luint i = 0; i < this->iterations; i++) {
        current = dd::simulate<>(&U_P, current, package);
        current = dd::simulate<>(&U_B, current, package);
    }
    clock_t a_iteration = clock();
    clock_t end = clock();

    // We store the data
    this->red_time = 0.0;
    this->it_time = time_to_double(b_iteration, a_iteration);
    this->tot_time = time_to_double(begin, end);
 
    return;
}

// PUBLIC METHODS
/* Method that runs the experiment */
void Experiment::run() {
    if (!this->executed) {
        cerr << "++ Executing experiment " << this->name << " of type " << this->type << endl;
        switch (this->type)
        {
        case ExperimentType::CLUE:
            run_clue();
            break;
        case ExperimentType::DDSIM:
            run_ddsim();
            break;
        case ExperimentType::DIRECT:
            run_direct();
            break;
        case ExperimentType::DDSIM_ALONE:
            run_ddsim_alone();
            break;
        default:
            throw logic_error("Unexpected value for type of Experiment");
        }
        this->executed = true; // We avoid repeating computations
    }
    return;
}
/* Method that generate the CSV row for this experiment */
string Experiment::to_csv(char delimiter = ',') {
    throw logic_error("Method 'to_csv' not yet implemented");
    stringstream stream;
    stream << this->size()     << delimiter <<
              this->name       << delimiter << 
              this->observable << delimiter <<
              this->red_time   << delimiter <<
              this->red_ratio  << delimiter <<
              this->iterations << delimiter <<
              this->it_time    << delimiter <<
              this->tot_time   << delimiter;

    return stream.str();
}