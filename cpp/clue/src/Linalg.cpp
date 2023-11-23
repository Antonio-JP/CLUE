#include <cmath>
#include "Linalg.hpp"

template <typename T>
SparseVector<T>::SparseVector(int dim) {
    if (dim <= 0) {
        throw std::invalid_argument("Vectors of non-positive dimension do not exist");
    }

    this->dim = dim;
}

template <typename T>
void SparseVector<T>::reduce(T coef, SparseVector<T>& other) { // method that computes inplace this - coef*other
    for (pair<int, T> ppair : other.nonzero) {
        T value = coef*ppair.second; // Coefficient coming from arguments
        // We check if the key is in "this"
        typename map<int,T>::const_iterator search = this->nonzero.find(ppair.first);
        if (search != this->nonzero.end()) { 
            // We compute the new value for this input
            value += search->second;

            if (value == (T)0) { // If the value becomes zero we remove it
                this->nonzero.erase(ppair.first);
                continue;
            }
        }
        // The variable value must be added to the nonzero map
        this->nonzero[ppair.first] = value;
    }
}

template <typename T>
T SparseVector<T>::inner_product(SparseVector<T>& other) {
    if (this->dim != other.dim) {
        throw std::invalid_argument("Inner product require two vectors of same dimension");
    }

    T total;
    // First we conjugate the other vector
    SparseVector<T> other_conj = other.conjugate();
    // Now we go through the intersection adding the products
    for (pair<int,T> this_pair : this->nonzero) {
        int key = this_pair.first;
        
        // We look for "key" in "other"
        typename map<int,T>::const_iterator search = other_conj.nonzero.find(key); 
        if (search != other_conj.nonzero.end()) {
            total += search->second;
        }
    }
    return total;
}

template <typename T>
T SparseVector<T>::get_value(int index) {
    // We check the index fit into the current dimension
    if (index < 0 || index >= this->dim) {
        throw std::invalid_argument("The index must be valid for the current vector dimension.");
    }
    // We look if "index" is in the "nonzero" map
    typename map<int,T>::const_iterator search = this->nonzero.find(index);
    if (search != this->nonzero.end()) { //Found
        return search->second;
    }
    // Else, we know it is zero
    return (T) 0;
}

template <typename T>
SparseVector<T> SparseVector<T>::operator+(SparseVector<T>& other) {
    SparseVector<T> new_vector = SparseVector<T>(*this);
    new_vector += other;
    return new_vector;
}
template <typename T>
void SparseVector<T>::operator+=(SparseVector<T>& other) {
    // We update the values of this using the data in 
    for (pair<int, T> ppair : other.nonzero) {
        int index = ppair.first;
        this->set_value(index, (*this)[index] + ppair.second);
    }
    return;
}
template <typename T>
SparseVector<T> SparseVector<T>::operator-() { /* unary minus: just negation of object */
    SparseVector<T> result = SparseVector<T>((*this) * ((T)-1));
    return result;
}
template <typename T>
SparseVector<T> SparseVector<T>::operator-(SparseVector<T>& other) {
    SparseVector<T> new_vector = SparseVector<T>(*this);
    new_vector -= other;
    return new_vector;
}
template <typename T>
void SparseVector<T>::operator-=(SparseVector<T>& other) {
    // We update the values of this using the data in 
    for (pair<int, T> ppair : other.nonzero) {
        int index = ppair.first;
        this->set_value(index, (*this)[index] - ppair.second);
    }
    return;
}
template <typename T>
SparseVector<T> SparseVector<T>::operator*(T other) { // Scalar product by a constant
    if (other == ((T) 0)) {
        return SparseVector<T>(this->dim);
    }
    // If the coefficient is not zero
    SparseVector<T> output = SparseVector<T>(*this); // We copy "this" and use in-place operations
    output *= other;
    return output;
}
template <typename T>
void SparseVector<T>::operator*=(T other) {
    if (other == ((T) 0)) {
        this->nonzero.clear();
    } else { // We know that other != 0, hence we only need to multiply things in the nonzero map
        for (pair<int,T> ppair : this->nonzero) {
            this->nonzero[ppair.first] = ppair.second * other;
        }
    }
}
template <typename T>
bool SparseVector<T>::operator==(SparseVector<T>& other) {
    if (this->nonzero_count() != other.nonzero_count()) { return false; } // The number of non-zeros is different
    for (pair<int, T> ppair : this->nonzero) {
        typename map<int,T>::const_iterator search = other.nonzero.find(ppair.first);
        if (search == other.nonzero.end()) { return false; } // Found something zero in other that is not zero in self

        if (ppair.second != search->second) { return false; } // Found a different element
    }
    return true; // All elements are equal
}

template <typename T>
SparseVector<T> SparseVector<T>::conjugate() { /* Conjugate the vector and return the new structure */
    SparseVector<T> output = SparseVector<T>(*this);
    output.conjugate_in();
    return output;
}
        
template class SparseVector<rational<int>>;
template class SparseVector<complex<double>>;