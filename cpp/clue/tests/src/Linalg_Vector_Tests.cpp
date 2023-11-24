
#include <iostream>

#include "Linalg.hpp"
#include <boost/rational.hpp>

void SV_Creation() {
    QQSparseVector dim1 = QQSparseVector(1);
    QQSparseVector dim2 = QQSparseVector(2);
    QQSparseVector dim3 = QQSparseVector(3);

    if (dim1.to_string() != "(0)") { throw std::runtime_error("Error creating a vector (dimension 1) : " + dim1.to_string()); }
    if (dim2.to_string() != "(0, 0)") { throw std::runtime_error("Error creating a vector (dimension 2) : " + dim2.to_string()); }
    if (dim3.to_string() != "(0, 0, 0)") { throw std::runtime_error("Error creating a vector (dimension 3) : " + dim3.to_string()); }
}


int main() {

    SV_Creation();

    return 0;
}