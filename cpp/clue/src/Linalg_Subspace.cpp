#include "Linalg.hpp"

#include <queue>
#include <vector>
#include <iostream>

/*******************************************************************************************************************
 * 
 * CLASS FOR SUBSPACE
 * 
********************************************************************************************************************/
/*********************************************************************/
/* ATTRIBUTE/PROPERTIES */
vector<double> CCSubspace::densities() {
    vector<double> result = vector<double>(this->dimension());

    for (int i = 0; i < this->dimension(); i++) {
        result[i] = this->basis[i].density();
    }

    return result;
}
vector<double> CCSubspace::norms() {
    vector<double> result = vector<double>(this->dimension());

    for (int i = 0; i < this->dimension(); i++) {
        result[i] = this->basis[i].norm();
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
    this->reduce_vector(&vector); // We reduce the vector
    if (! this->contains(vector)) {
        cout << "Adding new element with norm " << vector.norm() << endl;
        vector.normalize_in();
        this->basis.push_back(vector);
        return true;
    } else {
        cout << "Found an element inside: " << vector.norm() << endl;
    }
    return false;
}

/** 
 * Method that computes the minimal invariant subspace under some matrices.
 * 
 * Given a vector space `V` (given by `this`) and a set of matrices `(M_i)_i`, this method computes the minimal 
 * extension of `V` (say `W`) that is invariant under every matrix `M_i`, i.e., for every `v \in W`, `M_i v \in W`.
 * 
 * A matrix is a vector of vectors. In this case, we will use a ``std::vector`` of ``CCSparseVector`` as a row 
 * representation of the matrix. Hence, the multiplication `M_i v` can be done via inner products of the rows of `M_i`
 * with the vector `v` (which will also be defined with a CCSparseVector).
 * 
 * This method extends in-place the current subspace we are working with. This method returns the new dimension of the 
 * vector space.
*/
int CCSubspace::minimal_invariant_space(vector<vector<CCSparseVector>>& matrices) {
    // We create a queue with the first round we need to check
    queue<CCSparseVector> to_process;

    for (CCSparseVector current : this->basis) {
        for (vector<CCSparseVector> matrix : matrices) {
            // We do multiplication matrix*current
            CCSparseVector result = CCSparseVector(matrix.size()); // Dimension is number of rows of the matrix
            for (int i = 0; i < result.dimension(); i++) {
                result.set_value(i, matrix[i].inner_product(current));
            }
            // We add this vector to the queue
            to_process.push(result);
        }
    }

    //We now iterate on the queue until this is empty
    while ((!to_process.empty()) && (this->dimension() < this->ambient_dimension())) {
        cout << "Starting computation of minimal invariant with " << to_process.size() << " vectors. " << endl;
        CCSparseVector current = to_process.front(); to_process.pop(); // We take the first element
        
        cout << "Current vector: " << current.norm() << endl;
        bool absorbed = this->absorb_new_vector(current);
        cout << "Was absorbed: " << absorbed << endl;
        if (absorbed) { // We have increased the dimension, we need to add new vectors
            for (vector<CCSparseVector> matrix : matrices) {
                // We do multiplication matrix*current
                CCSparseVector result = CCSparseVector(matrix.size()); // Dimension is number of rows of the matrix
                for (int i = 0; i < result.dimension(); i++) {
                    result.set_value(i, matrix[i].inner_product(current));
                }
                // We add this vector to the queue
                to_process.push(result);
            } 
        }
    }

    return this->dimension();
}