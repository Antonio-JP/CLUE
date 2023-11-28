#include <cmath>
#include "Linalg.hpp"

/*******************************************************************************************************************
 * 
 * ABSTRACT CLASS FOR SPARSE VECTOR
 * 
********************************************************************************************************************/
/*******************************************************************************************************************/
/* CONSTRUCTORS */
template <typename T>
SparseVector<T>::SparseVector(int dim) {
    if (dim <= 0) {
        throw std::invalid_argument("Vectors of non-positive dimension do not exist");
    }

    this->dim = dim;
}

/*******************************************************************************************************************/
/* ATTRIBUTE/PROPERTIES */
template <typename T>
void SparseVector<T>::reduce(T coeff, SparseVector<T>& other) { // method that computes inplace this - coeff*other
    for (pair<int, T> ppair : other.nonzero) {
        T value = coeff*ppair.second; // Coefficient coming from arguments
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
    // Now we go through the intersection adding the products
    for (pair<int,T> this_pair : this->nonzero) {
        int key = this_pair.first;
        
        // We look for "key" in "other"
        typename map<int,T>::const_iterator search = other.nonzero.find(key); 
        if (search != other.nonzero.end()) {
            total += this_pair.second * this->conjugate_coeff(search->second);
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
    if (this->nonzero_indices.find(index) != this->nonzero_indices.end()) { // This search is O(1)
        return this->nonzero[index]; // We do a new search of O(log(N))
    }

    return (T) 0;
}

template <typename T>
void SparseVector<T>::set_value(int index, T value)  {
    if (value == (T) 0) {
        if (this->nonzero_indices.find(index) != this->nonzero_indices.end()) { // This search is O(1)
            this->nonzero.erase(index);
            this->nonzero_indices.erase(index);
        }
    } else {
        this->nonzero[index] = value;
        this->nonzero_indices.insert(index);
    }
}

template <typename T>
vector<T> SparseVector<T>::to_list() { 
    vector<T> dense;
    for (int i = 0; i < this->dim; i++) {
        dense[i] = (*this)[i];
    }
    return dense;
}

/*******************************************************************************************************************/
/* ARITHMETIC/OPERATOR METHODS */
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
void SparseVector<T>::operator-=(SparseVector<T>& other) {
    // We update the values of this using the data in 
    for (pair<int, T> ppair : other.nonzero) {
        int index = ppair.first;
        this->set_value(index, (*this)[index] - ppair.second);
    }
    return;
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
void SparseVector<T>::conjugate_in()  {
    for (pair<int,T> ppair : this->nonzero) {
        this->set_value(ppair.first, this->conjugate_coeff(ppair.second));
    }
}


/*******************************************************************************************************************
 * 
 * CLASS QQ-SPARSE VECTOR
 * 
********************************************************************************************************************/
/*******************************************************************************************************************/
/* VIRTUAL METHODS */
double QQSparseVector::norm() {
    QQ squared_norm = this->inner_product(*this);
    return sqrt(squared_norm.numerator() / (double) squared_norm.denominator());
}
QQSparseVector& QQSparseVector::normalize() {
    static QQSparseVector normalized = QQSparseVector(*this);
    normalized.normalize_in();

    return normalized;
}
void QQSparseVector::normalize_in() {
    if (this->is_zero()) { return; } // Nothing to do

    int gcd_n, gcd_d;
    map<int,QQ>::iterator iterator = this->nonzero.begin();
    gcd_n = iterator->second.numerator(); gcd_d = iterator->second.denominator();
    iterator++;
    while (iterator != this->nonzero.end()) {
        gcd_n = gcd(gcd_n, iterator->second.numerator());
        gcd_d = gcd(gcd_d, iterator->second.denominator());
        iterator++;
    }

    this->operator*=(QQ(gcd_d, gcd_n));
}
string QQSparseVector::coeff_to_string(QQ element) {
    string output = std::to_string(element.numerator());
    if (element.denominator() != 1 && element.numerator() != 0) {
        output += "/" + std::to_string(element.denominator());
    }
    return output;
}
/*******************************************************************************************************************/
/* ARITHMETIC METHODS */
QQSparseVector QQSparseVector::operator+(QQSparseVector& other) {
    QQSparseVector new_vector = QQSparseVector(*this);
    new_vector += other;
    return new_vector;
}
QQSparseVector QQSparseVector::operator-() { /* unary minus: just negation of object */
    QQSparseVector result = QQSparseVector((*this) * QQ(-1));
    return result;
}
QQSparseVector QQSparseVector::operator-(QQSparseVector& other) {
    QQSparseVector new_vector = QQSparseVector(*this);
    new_vector -= other;
    return new_vector;
}
QQSparseVector QQSparseVector::operator*(QQ other) { // Scalar product by a constant
    if (other == QQ(0)) {
        return QQSparseVector(this->dimension());
    }
    // If the coefficient is not zero
    QQSparseVector output = QQSparseVector(*this); // We copy "this" and use in-place operations
    output *= other;
    return output;
}

/*******************************************************************************************************************
 * 
 * CLASS CC-SPARSE VECTOR
 * 
********************************************************************************************************************/
/*******************************************************************************************************************/
/* VIRTUAL METHODS */
double CCSparseVector::norm()  {
    CC squared_norm = this->inner_product(*this);
    return sqrt(squared_norm.real());
}
CCSparseVector& CCSparseVector::normalize() {
    static CCSparseVector normalized = CCSparseVector(*this);
    normalized.normalize_in();

    return normalized;
}
void CCSparseVector::normalize_in() {
    /* In CC, we divide the coefficients by the norm */
    double norm = this->norm();
    CC c_norm = CC((double)(1/norm));

    this->operator*=(c_norm);
}
string CCSparseVector::coeff_to_string(CC element) {
    if (element.real() == (double) 0 && element.imag() == (double) 0) {
        return "0";
    } else if (element.real() == (double) 0) {
        return std::to_string(element.imag()) + "*i";
    } else if (element.imag() == (double) 0) {
        return std::to_string(element.real());
    } else {
        if (element.imag() < (double) 0) {
            return std::to_string(element.real()) + " - " + std::to_string(-element.imag()) + "*i";
        }
        return std::to_string(element.real()) + " + " + std::to_string(element.imag()) + "*i";
    }
}
/*******************************************************************************************************************/
/* ARITHMETIC METHODS */
CCSparseVector CCSparseVector::operator+(CCSparseVector& other) {
    CCSparseVector new_vector = CCSparseVector(*this);
    new_vector += other;
    return new_vector;
}
CCSparseVector CCSparseVector::operator-() { /* unary minus: just negation of object */
    CCSparseVector result = CCSparseVector((*this) * CC(-1));
    return result;
}
CCSparseVector CCSparseVector::operator-(CCSparseVector& other) {
    CCSparseVector new_vector = CCSparseVector(*this);
    new_vector -= other;
    return new_vector;
}
CCSparseVector CCSparseVector::operator*(CC other) { // Scalar product by a constant
    if (other == CC(0)) {
        return CCSparseVector(this->dimension());
    }
    // If the coefficient is not zero
    CCSparseVector output = CCSparseVector(*this); // We copy "this" and use in-place operations
    output *= other;
    return output;
}
CCSparseVector CCSparseVector::conjugate() { /* Conjugate the vector and return the new structure */
    CCSparseVector output = CCSparseVector(*this);
    output.conjugate_in();
    return output;
}
        
template class SparseVector<QQ>;
template class SparseVector<CC>;

