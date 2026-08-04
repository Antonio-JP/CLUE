#include <iostream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <bits/stdc++.h>
#include <time.h>
#include <cstdlib>
#include <string>

#include "experiments/Experiment.hpp"
#include "experiments/QASMExperiment.hpp"

using namespace std;

Experiment* generate_example(string name, string path, luint size, ExperimentType type, string observable, dd::Package<>* package) {
    string upper = boost::to_upper_copy<std::string>(name);
    return new QASMExperiment(size, name, path, observable, type, package);
}

vector<string> generate_observables(string observable, luint size) {
    vector<string> result = vector<string>();
    if (observable == "all") {
        result.push_back("H");
        for (luint i = 0; i < static_cast<luint>(pow(2, size)); i++) {
            result.push_back(std::to_string(i));
        }
    } else {
        result.push_back(observable);
    }
    return result;
}

int main_script(string name, string path, ExperimentType type, luint m, luint M, luint repeats, string observable) {
    double total_time = 0.;
    ofstream out; 
    // PROCESSING THE OUTPUT FILE
    string upper = boost::to_upper_copy<std::string>(name);
    string lower = boost::to_lower_copy<std::string>(name);
    string filename = "qasm";
    
    filesystem::path out_path = filesystem::path("../../../../tests/quantum/results/[result-cpp]q_" + filename + "_" + name + "_" + ExperimentType_toString(type) + ".csv");
    out.open(out_path, std::ios::app);
    cout << "##################################################################################" << endl;
    cout << "### EXECUTION ON " << boost::to_upper_copy<std::string>(name) << "[m=" << m << ", M=" << M << ", repeats=" << repeats << ", method=" << type << "]" << endl;
    cout << "##################################################################################" << endl;
    
    for (luint size = m; size <= M; size++) { // We repeat for each size
        for (string obs : generate_observables(observable, size)) {
            for (luint execution = 1; execution <= repeats; execution++) { // We repeat "repeats" times
                try {
                    dd::Package<>* package = new dd::Package<>(size);
                    Experiment * experiment = generate_example(name, path, size, type, obs, package);
                    cout << "Generated example\n\t" << experiment->to_string() << endl;
                    experiment->run();
                        
                    cout << "### -- Finished execution " << execution << "/" << repeats << "(size=" << size << "): took " << experiment->total_time() << "s." << endl;

                    total_time += experiment->total_time();
                    out << experiment->to_csv() << endl;
                    delete experiment;
                    delete package;
                } catch (qc::QFRException &e) {
                    cout << "### -- Error in execution " << execution << "/" << repeats << "(size=" << size << "): " << e.what() << endl;
                }
            }
        }
    }
    double average_time = total_time/static_cast<double>((M-m+1)*repeats);
    cout << "### Average execution time: " << average_time << endl;
    cout << "##################################################################################" << endl;
    return 0;
}

enum ArgumentValues {
    type, min, max, repeats, observable
};

std::map<std::string, ArgumentValues> create_argument_map() {
    std::map<std::string, ArgumentValues> m;
    m["-t"] = ArgumentValues::type;
    m["-m"] = ArgumentValues::min;
    m["-M"] = ArgumentValues::max;
    m["-repeats"] = ArgumentValues::repeats;
    m["-r"] = ArgumentValues::repeats;
    m["-obs"] = ArgumentValues::observable;
    return m;
}
static std::map<std::string, ArgumentValues> s_mapArgumentValues = create_argument_map();

int main(int argc, char** argv) {
    srand (static_cast<unsigned>(time(NULL)));
    string test;
    string path;
    ExperimentType type = ExperimentType::DDSIM;
    luint m = 9, M = 9, repeats = 1;
    string observable = "H";

    if (argc > 1) {
        test = argv[1];
        path = argv[2];
        int i = 3;
        while (i < argc) {
            switch (s_mapArgumentValues[argv[i]]) {
                case ArgumentValues::type:
                    type = ExperimentType_fromString(string(argv[i+1]));
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
                case ArgumentValues::observable:
                    observable = argv[i+1];
                    i+=2;
                    break;        
                default:
                    cout << "Error in arguments: found " << argv[i];
                    return -1;
            }
        }
    }

    return main_script(test, path, type, m, M, repeats, observable);
}