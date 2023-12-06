#include <iostream>
#include <bits/stdc++.h>
#include "Linalg.hpp"


using namespace std;

bool is_squared(vector<vector<dd::ComplexValue>>& M) {
    if (M.size()) {
        return M.size() == M[0].size();
    }
    return true;
}

vector<vector<dd::ComplexValue>> identity(luint size) {
    vector<vector<dd::ComplexValue>> result = vector<vector<dd::ComplexValue>>(size);
    for (luint i = 0; i < size; i++) {
        result[i] = vector<dd::ComplexValue>(size);
        result[i][i] = dd::ComplexValue(1);
    }
    return result;
}

vector<vector<dd::ComplexValue>> matrix_prod(vector<vector<dd::ComplexValue>>& A, vector<vector<dd::ComplexValue>>& B) {
    if (A[0].size() != B.size()) {
        throw logic_error("The dimensions do not match!");
    }
    vector<vector<dd::ComplexValue>> result = vector<vector<dd::ComplexValue>>(A.size());
    for (luint i = 0; i < A.size(); i++) {
        result[i] = vector<dd::ComplexValue>(B[0].size()); // Number of columns of B
        for (luint j = 0; j < B[0].size(); j++) {
            for (luint k = 0; k < B.size(); k++) {
                result[i][j] += A[i][k]*B[k][i];
            }
        }
    }
    return result;
}

vector<vector<dd::ComplexValue>> Power(vector<vector<dd::ComplexValue>>& M, int t) {
    if (! is_squared(M)) {
        throw logic_error("The matrix must be square");
    }
    vector<vector<dd::ComplexValue>> R, I = identity(M.size()), B;
    int i;
    if(t & 1) {        //Handle odd values of t (this saves a multiplication later)
        R = M;
        t = t & ~1;    //Clear the least significant bit of t
    } else {
        R = I;
    }
    i=1;
    B=M;                //B will always be M^i, where i is a power of 2
    while (t!=0) {
        i = i*2;         //Advance i to the next power of 2
        B = matrix_prod(B,B);         //B was M^(i/2) and is now M^i

        if(t & i) {       //i is of the form 2^j. Is the j-th bit of t set?
            R = matrix_prod(R,B);      //Multiply the result with B=A^i
            t = t & ~i;   //Clear the j-th bit of t
        }
    }

    return R;
}

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
    for (luint i = N2; i < N; i++) {
        output[i] = CC(1/sqrt(N2));
    }

    return output;
}
dd::CVec starting_phi(luint nQbits) {
    luint N = static_cast<luint>(pow(2, nQbits));
    dd::CVec output = dd::CVec(N);
    for (luint i = 0; i < N; i++) {
        output[i] = CC(1/sqrt(N));
    }

    return output;
}

luint run_example(luint nQbits, dd::CVec& starting, dd::CMat& U) {
    clock_t start = clock();
    cout << "Translating from dense to DD..." << endl;
    dd::Package<> * package = dd_package(nQbits);
    dd::vEdge starting_DD = package->makeStateFromVector(starting);
    dd::mEdge circuit = package->makeDDFromMatrix(U);
    vector<dd::mEdge> circuits = {circuit};
    luint N = static_cast<luint>(pow(static_cast<luint>(2), nQbits));

    clock_t reading = clock();

    cout << "Creating the starting subspace..." << endl;
    DDSubspace lumping = DDSubspace(N);
    DDVector starting_vector = DDVector(nQbits, starting_DD);
    lumping.absorb_new_vector(&starting_vector);

    cout << "Computing the lumping..." << endl;
    lumping.minimal_invariant_space(circuits);
    vector<vector<dd::ComplexValue>> Uhat = lumping.reduced_matrix(circuit);

    clock_t time_lumping = clock();

    Power(Uhat, static_cast<int>(floor(pow(2., static_cast<double>(nQbits)/2.))));

    clock_t simulating = clock();

    dd::CMat inner = dd::CMat(lumping.dimension());
    for (luint i = 0; i < lumping.dimension(); i++) {
        inner[i] = dd::CVec(lumping.dimension());
        for (luint j = 0; j < lumping.dimension(); j++) {
            inner[i][j] = lumping.basis[i]->inner_product(*(lumping.basis[j]));
        }
    }
    
    cout << "Inner product matrix:" << endl;
    cout << "-----------------------------------------------------------------------------------" << 
            endl << 
            matrix_to_string(inner) << 
            "-----------------------------------------------------------------------------------" << 
            endl;
    
    cout << "Reduced matrix:" << endl;
    cout << "-----------------------------------------------------------------------------------" << 
            endl << 
            matrix_to_string(Uhat) << 
            "-----------------------------------------------------------------------------------" << 
            endl;
    
    clock_t end = clock();
    // Counting the time spent in the function
    cout << "Total time execution:\t" << (double(end - start) / double(CLOCKS_PER_SEC)) << endl;
    cout << "Total time reading:\t" << (double(reading - start) / double(CLOCKS_PER_SEC)) << endl;
    cout << "Total time lumping:\t" << (double(time_lumping - reading) / double(CLOCKS_PER_SEC)) << endl;
    cout << "Total time simulating:\t" << (double(simulating - time_lumping) / double(CLOCKS_PER_SEC)) << endl;

    return lumping.dimension();
}

luint example_Grover_1() { // Should be 2
    //Toffoli of size 10
    cout << "Reading the data and creating the basic structures..." << endl;
    luint nQbits = 11;
    dd::CVec initial = starting_Grover(nQbits);
    dd::CMat U = read_matrix("../../../../tests/quantum/Grover-11.txt");

    cout << "Created a vector of size " << initial.size() << endl;
    cout << "Creating a matrix of size (" << U.size() << "," << U[0].size() << ")" << endl;

    return run_example(nQbits, initial, U);
}

luint example_Grover_2() { // Should be 2
    //Toffoli of size 3
    cout << "Reading the data and creating the basic structures..." << endl;
    luint nQbits = 4;
    dd::CVec initial = starting_Grover(nQbits);
    dd::CMat U = read_matrix("../../../../tests/quantum/Grover-4.txt");

    cout << "Created a vector of size " << initial.size() << endl;
    cout << "Creating a matrix of size (" << U.size() << "," << U[0].size() << ")" << endl;

    return run_example(nQbits, initial, U);
}

int main(int, char**) {
    example_Grover_1(); // Got
    // example_Grover_2(); // Got 2
    return 0;
}