#include <iostream>
#include "Linalg.hpp"


using namespace std;

dd::CMat read_matrix(string filename) {
    std::ifstream myfile; myfile.open(filename);
    string line;
    dd::CMat output;
    while( std::getline(myfile, line) ) {
        // Processing the line
        dd::CVec row;
        std::istringstream line_stream (line);
        string element;
        while(std::getline(line_stream, element, ' ')) { // We split by space
            luint c = element.find_first_of(':');
            double real = std::stod(element.substr(0, c));
            double imag = std::stod(element.substr(c+1));
            row.push_back(CC(real,imag));
        }
        output.push_back(row);
    }

    return output;
}

dd::CVec starting_Grover(luint nQbits) {
    luint N = static_cast<luint>(pow(2, nQbits)), N2 = N/2;
    dd::CVec output = dd::CVec(N);
    for (luint i = 0; i < N2; i++) {
        output[i] = CC(1);
    }

    return output;
}
dd::CVec starting_phi(luint nQbits) {
    luint N = static_cast<luint>(pow(2, nQbits));
    dd::CVec output = dd::CVec(N);
    for (luint i = 0; i < N; i++) {
        output[i] = CC(1);
    }

    return output;
}

luint run_example(luint nQbits, dd::CVec& starting, dd::CMat& U) {
    cout << "Translating from dense to DD..." << endl;
    dd::Package<> * package = dd_package(nQbits);
    dd::vEdge starting_DD = package->makeStateFromVector(starting);
    dd::mEdge circuit = package->makeDDFromMatrix(U);
    vector<dd::mEdge> circuits = {circuit};
    luint N = static_cast<luint>(pow(static_cast<luint>(2), nQbits));

    cout << "Creating the starting subspace..." << endl;
    DDSubspace lumping = DDSubspace(N);
    DDVector starting_vector = DDVector(nQbits, starting_DD);
    lumping.absorb_new_vector(starting_vector);

    cout << "Computing the lumping..." << endl;
    lumping.minimal_invariant_space(circuits);

    cout << "Found a reduction of size: " << lumping.dimension() << endl;    
    cout << "Norms of the vectors in lumping: " << endl;
    cout << "[";
    vector<double> lumping_norms = lumping.norms();
    for (double fl : lumping_norms) {
        cout << fl << ", ";
    }
    cout << "]" << endl;

    CC **inner = NULL;
    inner = static_cast<CC**>(calloc(sizeof(CC *), lumping.dimension()));
    for (luint i = 0; i < lumping.dimension(); i++) {
        inner[i] = static_cast<CC *>(calloc(sizeof(CC),lumping.dimension()));
        for (luint j = 0; j < lumping.dimension(); j++) {
            inner[i][j] = lumping.basis[i].inner_product(lumping.basis[j]);
        }
    }
    
    cout << "Inner product matrix:" << endl;
    cout << "-----------------------------------------------------------------------------------" << endl;
    for (luint i = 0; i < lumping.dimension(); i++) {
        for (luint j = 0; j <lumping.dimension(); j++) {
            cout << CC_to_string(inner[i][j]) << ",\t";
        }
        cout << endl;
    }
    cout << "-----------------------------------------------------------------------------------" << endl;

    // Freeing memory
    for (luint i = 0; i < lumping.dimension(); i++) {
        free(inner[i]);
    }
    free(inner);

    return lumping.dimension();
}

luint example_Grover_1() { // Should be 2 or 3?
    //Toffoli of size 10
    cout << "Reading the data and creating the basic structures..." << endl;
    luint nQbits = 11;
    dd::CVec initial = starting_Grover(nQbits);
    dd::CMat U = read_matrix("../../../../tests/quantum/Grover-11.txt");

    cout << "Created a vector of size " << initial.size() << endl;
    cout << "Creating a matrix of size (" << U.size() << "," << U[0].size() << ")" << endl;

    return run_example(nQbits, initial, U);
}

int main(int, char**) {

    if(dd_package(5) != dd_package(5)) {
        return -1;
    }
    example_Grover_1(); // Got
    return 0;
}