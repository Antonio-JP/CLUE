
#include <iostream>
#include <vector>

#include "Linalg.hpp"
#include <boost/rational.hpp>

void SV_QQ_Creation() {
    QQSparseVector dim1 = QQSparseVector(1);
    QQSparseVector dim2 = QQSparseVector(2);
    QQSparseVector dim3 = QQSparseVector(3);

    if (dim1.to_string() != "(0)") { throw std::runtime_error("Error creating a vector (dimension 1) : " + dim1.to_string()); }
    if (dim2.to_string() != "(0, 0)") { throw std::runtime_error("Error creating a vector (dimension 2) : " + dim2.to_string()); }
    if (dim3.to_string() != "(0, 0, 0)") { throw std::runtime_error("Error creating a vector (dimension 3) : " + dim3.to_string()); }

    vector<QQ> input = {QQ(1), QQ(2), QQ(3)};
    QQSparseVector from_vector = QQSparseVector(input);
    if (from_vector.to_string() != "(1, 2, 3)") { throw std::runtime_error("Error creating a vector (from vector) : " + from_vector.to_string()); }

    QQSparseVector from_SV = QQSparseVector(from_vector);
    if (from_SV.to_string() != "(1, 2, 3)") { throw std::runtime_error("Error creating a vector (from vector) : " + from_vector.to_string()); }
}

void SV_QQ_Arithmetic() {
    QQSparseVector u = QQSparseVector({QQ(1), QQ(2), QQ(3)});
    QQSparseVector v = QQSparseVector({QQ(1,2), QQ(3,2), QQ(5,2)});
    QQSparseVector result = u+v;

    if (result[0] != QQ(3,2)) { throw std::runtime_error("Error in the sum : " + result.to_string()); }
    if (result[1] != QQ(7,2)) { throw std::runtime_error("Error in the sum : " + result.to_string()); }
    if (result[2] != QQ(11,2)) { throw std::runtime_error("Error in the sum : " + result.to_string()); }

    result = u-v;
    if (result[0] != QQ(1,2)) { throw std::runtime_error("Error in the difference : " + result.to_string()); }
    if (result[1] != QQ(1,2)) { throw std::runtime_error("Error in the difference : " + result.to_string()); }
    if (result[2] != QQ(1,2)) { throw std::runtime_error("Error in the difference : " + result.to_string()); }

    QQSparseVector diff = u-v;
    result += diff;
    if (result[0] != QQ(1)) { throw std::runtime_error("Error in the in-place sum : " + result.to_string()); }
    if (result[1] != QQ(1)) { throw std::runtime_error("Error in the in-place sum : " + result.to_string()); }
    if (result[2] != QQ(1)) { throw std::runtime_error("Error in the in-place sum : " + result.to_string()); }

    if (diff*QQ(2) == result) { }// all good
    else { throw std::runtime_error("Error in equality"); }
    if (diff*QQ(2) != result) { throw std::runtime_error("Error in equality"); }

    if (diff == result) { throw std::runtime_error("Error in equality"); }
    QQSparseVector conj = result.conjugate();
    if (conj != result) { throw std::runtime_error("Error in conjugation"); }
}

void SV_QQ_TrueSparse() {
    QQSparseVector u = QQSparseVector(10), v = QQSparseVector(10);
    u.set_value(0, QQ(1)); u.set_value(5, QQ(2,5));

    
    if (u.density() != (float).2) { throw std::runtime_error("Error in the density : " + to_string(u.density())); }
    if (v.density() != (float) 0) { throw std::runtime_error("Error in the density : " + to_string(v.density())); }

    u.set_value(0, QQ()); // Removes the element 0
    if (u.density() != (float).1) { throw std::runtime_error("Error in the density : " + to_string(u.density())); }

    if (u[0] != QQ()) { throw std::runtime_error("Error getting a non-given value : " + u.to_string());  }
    if (u[2] != QQ()) { throw std::runtime_error("Error getting a non-given value : " + u.to_string());  }
    if (u[5] != QQ(2,5)) { throw std::runtime_error("Error getting a given value : " + u.to_string());  }
}

void SV_QQ_InnerProduct() {
    QQSparseVector u = QQSparseVector({QQ(1),QQ(),QQ(),QQ(2)}), v = QQSparseVector({QQ(1),QQ(),QQ(3),QQ()});
    QQSparseVector w = QQSparseVector({QQ(),QQ(1),QQ(5),QQ()}), x = QQSparseVector({QQ(),QQ(),QQ(1),QQ(7)});

    if (u.inner_product(v) != QQ(1) ) { throw std::runtime_error("Error in the inner_product : " + u.coeff_to_string(u.inner_product(v))); }
    if (u.inner_product(w) != QQ()  ) { throw std::runtime_error("Error in the inner_product : " + u.coeff_to_string(u.inner_product(w))); }
    if (u.inner_product(x) != QQ(14)) { throw std::runtime_error("Error in the inner_product : " + u.coeff_to_string(u.inner_product(x))); }
    if (v.inner_product(u) != QQ(1) ) { throw std::runtime_error("Error in the inner_product : " + u.coeff_to_string(v.inner_product(u))); }
    if (w.inner_product(u) != QQ()  ) { throw std::runtime_error("Error in the inner_product : " + u.coeff_to_string(w.inner_product(u))); }
    if (x.inner_product(u) != QQ(14)) { throw std::runtime_error("Error in the inner_product : " + u.coeff_to_string(x.inner_product(u))); }
    
    if (v.inner_product(w) != QQ(15)) { throw std::runtime_error("Error in the inner_product : " + u.coeff_to_string(u.inner_product(v))); }
    if (v.inner_product(x) != QQ(3) ) { throw std::runtime_error("Error in the inner_product : " + u.coeff_to_string(u.inner_product(v))); }

    if (w.inner_product(x) != QQ(5) ) { throw std::runtime_error("Error in the inner_product : " + u.coeff_to_string(u.inner_product(v))); }
}

void SV_QQ_Normalize() {
    QQSparseVector u = QQSparseVector({QQ(2,3),QQ(10,9),QQ(14,27)});
    QQSparseVector normalized = u.normalize();
    QQSparseVector true_normalized = QQSparseVector({QQ(1),QQ(5,3),QQ(7,9)});

    if (normalized != true_normalized) { throw std::runtime_error("Error in QQ normalization: " + normalized.to_string()); }

    u.normalize_in();
    if (u != normalized) { throw std::runtime_error("Error in inplace QQ normalization: " + normalized.to_string()); }

}

void SV_CC_Creation() {
    CCSparseVector dim1 = CCSparseVector(1);
    CCSparseVector dim2 = CCSparseVector(2);
    CCSparseVector dim3 = CCSparseVector(3);

    if (dim1.to_string() != "(0)") { throw std::runtime_error("Error creating a vector (dimension 1) : " + dim1.to_string()); }
    if (dim2.to_string() != "(0, 0)") { throw std::runtime_error("Error creating a vector (dimension 2) : " + dim2.to_string()); }
    if (dim3.to_string() != "(0, 0, 0)") { throw std::runtime_error("Error creating a vector (dimension 3) : " + dim3.to_string()); }

    vector<CC> input = {CC(1), CC(2), CC(3)};
    CCSparseVector from_vector = CCSparseVector(input);
    if (from_vector.to_string() != "(1.000000, 2.000000, 3.000000)") { throw std::runtime_error("Error creating a vector (from vector) : " + from_vector.to_string()); }

    CCSparseVector from_SV = CCSparseVector(from_vector);
    if (from_SV.to_string() != "(1.000000, 2.000000, 3.000000)") { throw std::runtime_error("Error creating a vector (from vector) : " + from_vector.to_string()); }

    CCSparseVector complex_vector = CCSparseVector({CC(1,1), CC(2,-1), CC(3), CC(5,-2)});
    if (complex_vector.to_string() != "(1.000000 + 1.000000*i, 2.000000 - 1.000000*i, 3.000000, 5.000000 - 2.000000*i)") { throw std::runtime_error("Error creating a vector (from vector) : " + complex_vector.to_string()); }
}

void SV_CC_Arithmetic() {
    CCSparseVector u = CCSparseVector({CC(1), CC(2), CC(3)});
    CCSparseVector v = CCSparseVector({CC(1,2), CC(3,2), CC(5,2)});
    CCSparseVector result = u+v;

    if (result[0] != CC(2,2)) { throw std::runtime_error("Error in the sum : " + result.to_string()); }
    if (result[1] != CC(5,2)) { throw std::runtime_error("Error in the sum : " + result.to_string()); }
    if (result[2] != CC(8,2)) { throw std::runtime_error("Error in the sum : " + result.to_string()); }

    result = u-v;
    if (result[0] != CC(0,-2)) { throw std::runtime_error("Error in the difference : " + result.to_string()); }
    if (result[1] != CC(-1,-2)) { throw std::runtime_error("Error in the difference : " + result.to_string()); }
    if (result[2] != CC(-2,-2)) { throw std::runtime_error("Error in the difference : " + result.to_string()); }

    CCSparseVector diff = u-v;
    result += diff;
    if (result[0] != CC(0,-4)) { throw std::runtime_error("Error in the in-place sum : " + result.to_string()); }
    if (result[1] != CC(-2,-4)) { throw std::runtime_error("Error in the in-place sum : " + result.to_string()); }
    if (result[2] != CC(-4,-4)) { throw std::runtime_error("Error in the in-place sum : " + result.to_string()); }

    if (diff*CC(2) == result) { }// all good
    else { throw std::runtime_error("Error in equality"); }
    if (diff*CC(2) != result) { throw std::runtime_error("Error in equality"); }

    CCSparseVector scaled = CCSparseVector({CC(-2,1), CC(-2,3), CC(-2,5)});
    CCSparseVector conjugated = CCSparseVector({CC(-2,-1), CC(-2,-3), CC(-2,-5)});
    result = v*CC(0,1);
    if (result != scaled) { throw std::runtime_error("Error in the scalar product : " + result.to_string()); }

    if (diff == result) { throw std::runtime_error("Error in equality"); }
    CCSparseVector conj = result.conjugate();
    if (conj != conjugated) { throw std::runtime_error("Error in conjugation"); }
}

void SV_CC_TrueSparse() {
    CCSparseVector u = CCSparseVector(10), v = CCSparseVector(10);
    u.set_value(0, CC(1)); u.set_value(5, CC(2,5));

    
    if (u.density() != (float).2) { throw std::runtime_error("Error in the density : " + to_string(u.density())); }
    if (v.density() != (float) 0) { throw std::runtime_error("Error in the density : " + to_string(v.density())); }

    u.set_value(0, CC()); // Removes the element 0
    if (u.density() != (float).1) { throw std::runtime_error("Error in the density : " + to_string(u.density())); }

    if (u[0] != CC()) { throw std::runtime_error("Error getting a non-given value : " + u.to_string());  }
    if (u[2] != CC()) { throw std::runtime_error("Error getting a non-given value : " + u.to_string());  }
    if (u[5] != CC(2,5)) { throw std::runtime_error("Error getting a given value : " + u.to_string());  }
}

void SV_CC_InnerProduct() {
    CCSparseVector u = CCSparseVector({CC(1),CC(),CC(),CC(2,1)}), v = CCSparseVector({CC(1),CC(),CC(3,1),CC()});
    CCSparseVector w = CCSparseVector({CC(),CC(1),CC(5,1),CC()}), x = CCSparseVector({CC(),CC(),CC(1),CC(7,1)});

    if (u.inner_product(v) != CC(1) ) { throw std::runtime_error("Error in the inner_product (1) : " + u.coeff_to_string(u.inner_product(v))); }
    if (u.inner_product(w) != CC()  ) { throw std::runtime_error("Error in the inner_product (2) : " + u.coeff_to_string(u.inner_product(w))); }
    if (u.inner_product(x) != CC(15,5)) { throw std::runtime_error("Error in the inner_product (3) : " + u.coeff_to_string(u.inner_product(x))); }
    if (v.inner_product(u) != CC(1) ) { throw std::runtime_error("Error in the inner_product (4) : " + u.coeff_to_string(v.inner_product(u))); }
    if (w.inner_product(u) != CC()  ) { throw std::runtime_error("Error in the inner_product (5) : " + u.coeff_to_string(w.inner_product(u))); }
    if (x.inner_product(u) != CC(15,-5)) { throw std::runtime_error("Error in the inner_product (6) : " + u.coeff_to_string(x.inner_product(u))); }
    
    if (v.inner_product(w) != CC(16,2)) { throw std::runtime_error("Error in the inner_product (7) : " + u.coeff_to_string(v.inner_product(w))); }
    if (v.inner_product(x) != CC(3,1) ) { throw std::runtime_error("Error in the inner_product (8) : " + u.coeff_to_string(v.inner_product(x))); }
    if (w.inner_product(v) != CC(16,-2)) { throw std::runtime_error("Error in the inner_product (9) : " + u.coeff_to_string(w.inner_product(v))); }
    if (x.inner_product(v) != CC(3,-1) ) { throw std::runtime_error("Error in the inner_product (10) : " + u.coeff_to_string(x.inner_product(v))); }

    if (w.inner_product(x) != CC(5,1) ) { throw std::runtime_error("Error in the inner_product (11) : " + u.coeff_to_string(w.inner_product(x))); }
    if (x.inner_product(w) != CC(5,-1) ) { throw std::runtime_error("Error in the inner_product (12) : " + u.coeff_to_string(x.inner_product(w))); }
}

void SV_CC_Normalize() {
    CCSparseVector u = CCSparseVector({CC(2,3),CC(10,9),CC(14,27)});
    CCSparseVector normalized = u.normalize();
    CCSparseVector true_normalized = u*(CC((double) 1/u.norm()));

    if ((normalized - true_normalized).norm() > 1e-8) { throw std::runtime_error("Error in CC normalization: " + (normalized-true_normalized).to_string()); }

    u.normalize_in();
    if (u != normalized) { throw std::runtime_error("Error in inplace CC normalization: " + normalized.to_string()); }

}

int main() {

    SV_QQ_Creation();
    SV_QQ_Arithmetic();
    SV_QQ_TrueSparse();
    SV_QQ_InnerProduct();
    SV_QQ_Normalize();

    SV_CC_Creation();
    SV_CC_Arithmetic();
    SV_CC_TrueSparse();
    SV_CC_InnerProduct();
    SV_CC_Normalize();

    return 0;
}