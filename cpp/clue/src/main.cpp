#include <iostream>
#include "Linalg.hpp"

using namespace std;

int main(int, char**) {
    QQSparseVector dim1 = QQSparseVector(1);

    if (dim1.dimension() != 1) {
        cout << "Found an error on dimension?" << endl;
    }
    cout << "Got a dimension " << dim1.dimension() << " vector. Really?" << endl;
}