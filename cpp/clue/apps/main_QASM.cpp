#include <iostream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <bits/stdc++.h>

#include "QuantumComputation.hpp"
#include "dd/Package.hpp"
#include "Linalg.hpp"

using namespace std;

double ddsim(string name, luint size, string observable) {
    // Loading the circuit
    cout << "### ++ -- Reading the example " << name << endl;
    clock_t b_read = clock();
    qc::QuantumComputation circuit = qc::QuantumComputation("../../../../tests/quantum/circuits/" + name + ".qasm");
    clock_t a_read = clock();
    double read_time = double(a_read - b_read) / double(CLOCKS_PER_SEC);
    std::unique_ptr<dd::Package<>> package = std::make_unique<dd::Package<>>(size);
    vector<qc::QuantumComputation> circuits = {circuit};
    cout << "### ++ -- -- Read " << name << endl;

    // Creating the initial state
    clock_t b_init = clock();
    dd::vEdge init;
    if (observable == "H") {
        init = package->makeBasisState(size, vector<dd::BasisStates>(size, dd::BasisStates::plus), 0);
    } else {
        int val = stoi(observable);
        vector<bool> binary = vector<bool>(size, false);
        for (luint i = size; i > 0 && val > 0; i--) {
            binary[i-1] = (val % 2 == 1);
            val /= 2;
        }
        init = package->makeBasisState(size, binary, 0);
    }
    clock_t a_init = clock();
    double init_time = double(a_init - b_init) / double(CLOCKS_PER_SEC);

    // Computing the lumping
    clock_t b_lumping = clock();
    DDSubspace lumping = DDSubspace(size, package);
    lumping.absorb_new_vector(&init);

    lumping.minimal_invariant_space(circuits);
    clock_t a_lumping = clock();
    double lumping_time = double(a_lumping - b_lumping) / double(CLOCKS_PER_SEC);

    // Computing the reduced model
    clock_t b_reduce = clock();
    vector<vector<dd::ComplexValue>> Uhat = lumping.reduced_matrix(circuit);
    clock_t a_reduce = clock();
    double reducing_time = double(a_reduce - b_reduce) / double(CLOCKS_PER_SEC);

    // Return time
    cout << "### ++ -- Execution times:" << endl;
    cout << "### ++ -- \tReading  : " << read_time << endl;
    cout << "### ++ -- \tInital   : " << init_time << endl;
    cout << "### ++ -- \tLumping  : " << lumping_time << endl;
    cout << "### ++ -- \tReducing : " << reducing_time << endl;
    cout << "### ++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "### ++ -- Lumping size: " << lumping.dimension() << endl;
    
    return read_time + init_time + lumping_time + reducing_time;
}

int main_script(string name, string type, luint m, luint M, luint repeats, vector<string> observables) {
    double total_time = 0., current_time = 0.;
    cout << "##################################################################################" << endl;
    cout << "### EXECUTION ON " << boost::to_upper_copy<std::string>(name) << "[m=" << m << ", M=" << M << ", repeats=" << repeats << ", method=" << type << "]" << endl;
    cout << "##################################################################################" << endl;
    
    for (luint size = m; size <= M; size++) { // We repeat for each size
        for (luint execution = 1; execution <= repeats; execution++) { // We repeat "repeats" times
            for (luint i = 0; i < observables.size(); i++) { // We execute each observable
                string observable = observables[i];
                cout << "### ++ Starting execution " << execution << "/" << repeats << "(size=" << size << ", observable=" << i+1 << "/" << observables.size() << ")" << endl;
                // We distinguish each case
                if (type == "ddsim") {
                    current_time = ddsim(name, size, observable);
                } else {
                    return -1;
                }
                
                cout << "### -- Finished execution " << execution << "/" << repeats << "(size=" << size << ", observable=" << i+1 << "/" << 
                        observables.size() << "): took " << current_time << "s." << endl;
                total_time += current_time;
            }
        }
    }
    cout << "### Average execution time: " << total_time/static_cast<double>((M-m+1)*repeats*observables.size()) << endl;
    cout << "##################################################################################" << endl;
    return 0;
}

enum ArgumentValues {
    type, min, max, repeats
};

std::map<std::string, ArgumentValues> create_argument_map() {
    std::map<std::string, ArgumentValues> m;
    m["-t"] = ArgumentValues::type;
    m["-m"] = ArgumentValues::min;
    m["-M"] = ArgumentValues::max;
    m["-repeats"] = ArgumentValues::repeats;
    return m;
}
static std::map<std::string, ArgumentValues> s_mapArgumentValues = create_argument_map();

int main(int argc, char** argv) {
    string file = "maxcut_3_2", type = "ddsim";
    luint m = 0, M = 0, repeats = 1;

    if (argc > 1) {
        file = argv[1];
        int i = 2;
        while (i < argc) {
            switch (s_mapArgumentValues[argv[i]]) {
                case ArgumentValues::type:
                    type = argv[i+1];
                    i+=2;
                    break;
                case ArgumentValues::min:
                    m = stoul(argv[i+1]);
                    i+=2;
                    break;
                case ArgumentValues::max:
                    M = stoul(argv[i+1]);
                    i+=2;
                    break;
                case ArgumentValues::repeats:
                    repeats = stoul(argv[i+1]);
                    i+=2;
                    break;        
                default:
                    cout << "Error in arguments: found " << argv[i];
                    return -1;
            }
        }
    }

    return main_script(file, type, m, M, repeats, vector<string>({"H"}));;
}