#include <iostream>
#include "Linalg.hpp"

using namespace std;

int main(int, char**) {
    QQSparseVector u = QQSparseVector({QQ(1),QQ(),QQ(),QQ(2)});

    cout << "My vector density: " << u.coeff_to_string(QQ(4,2)) << endl;
}