#include <vector>

#include "Linalg.hpp"

/*******************************************************************************************************************
 * 
 * CLASS FOR SUBSPACE
 * 
********************************************************************************************************************/
/*********************************************************************/
/* ATTRIBUTE/PROPERTIES */
vector<float> CCSubspace::densities() {
    vector<float> result = vector<float>(this->dimension());

    for (int i = 0; i < this->dimension(); i++) {
        result[i] = this->basis[i].density();
    }

    return result;
}

/*********************************************************************/
/* GETTING/SETTING DATA METHODS */
void CCSubspace::reduce_vector(CCSparseVector* vector) {
    /* Method that reduced a vector according to 'this' in-place. */
    for (int i = 0; i < this->dimension(); i++) {
        CCSparseVector to_rem = this->basis[i] * this->basis[i].inner_product(*vector);
        vector->operator-=(to_rem);
    }
}
bool CCSubspace::contains(CCSparseVector& vector) {
    /* Returns whether a vector is in the space or not */
    CCSparseVector copy = CCSparseVector(vector);
    this->reduce_vector(&copy);
    return copy.norm() < this->max_error;
}

CCSparseVector CCSubspace::find_in(CCSparseVector& vector) {
    /* Returns a vector representing the coordinates of vector in "this"*/
    if (! this->contains(vector)) { throw std::logic_error("The vector is not in the space"); }
    CCSparseVector result = CCSparseVector(this->dimension());
    for (int i = 0; i < this->dimension(); i++) {
        result.set_value(i, this->basis[i].inner_product(vector));
    }

    return result;
}

/**
 * Method that adds (if necessary) a new vector to "this"
 * 
 * This method guarantees that the new vector added is in "normal form",
 * which means that it is fully reduced w.r.t. "this" and with norm 1.
 * This guarantees some other properties of the basis.
*/
bool CCSubspace::absorb_new_vector(CCSparseVector& vector) {
    CCSparseVector copy = CCSparseVector(vector);

    this->reduce_vector(&copy); // We reduce the vector
    if (! this->contains(copy)) {
        this->basis[this->dimension()] = copy;
        return true;
    }
    return false;
}