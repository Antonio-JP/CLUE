#include <iostream>
#include "RationalFunctions.h"

using namespace std;

string say_hello(){
    return "Hello, from clue!\n";
}

int main(int, char**) {
    vector<int> degrees = {5, 4, 3, 0, 1};
    Monomial mon = Monomial(degrees);
    vector<string> varnames = {"x", "y", "z", "a", "b"};
    cout << "My monomial: " << mon.to_string(varnames) << endl;;
}