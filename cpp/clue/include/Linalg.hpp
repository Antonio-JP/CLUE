#ifndef CLUE_LA_H
#define CLUE_LA_H

#include <map>
#include <vector>
#include <complex>
#include <sstream>
#include <string>
#include <boost/rational.hpp>


using namespace std;
using namespace boost;

/* Class for vector spaces and vectors */
template <typename T>
class SparseVector {
    private:
        int dim;
    protected:
        map<int, T> nonzero;
    public:
        /*********************************************************************/
        /* CONSTRUCTORS */
        SparseVector(int dim);
        SparseVector(vector<T> dense_vector) : SparseVector(dense_vector.size()) {
            for (int i = 0; i < this->dim; i++) {
                if (dense_vector[i] != (T)0) {
                    this->set_value(i, dense_vector[i]);
                }
            }
        }
        SparseVector(SparseVector<T>& sparse_vector) : SparseVector(sparse_vector.dim) {
            for (pair<int, T> ppair : sparse_vector.nonzero) {
                this->set_value(ppair.first, ppair.second);
            }
        }
        
        /*********************************************************************/
        /* ATTRIBUTE/PROPERTIES */
        int nonzero_count() { return this->nonzero.size(); }
        int first_nonzero() { return this->nonzero.begin()->first; }
        float density() { return (float)this->nonzero_count() / (float)this->dim; }

        bool is_zero() { return this->nonzero_count() == 0; }
        
        void reduce(T coef, SparseVector<T>& other); // method that computes inplace this - coef*other
        T inner_product(SparseVector<T>&);
        virtual float norm() { throw std::logic_error("Method 'norm' not implemented"); }

        /*********************************************************************/
        /* GETTING/SETTING DATA METHODS */
        T get_value(int index);
        void set_value(int index, T value) {
            if (value == (T) 0) {
                typename map<int,T>::const_iterator search = this->nonzero.find(index);
                if (search != this->nonzero.end()) {
                    this->nonzero.erase(index);
                }
            } else {
                this->nonzero[index] = value;
            }
        }

        vector<T> to_list() { 
            vector<T> dense;
            for (int i = 0; i < this->dim; i++) {
                dense[i] = (*this)[i];
            }
            return dense;
        }

        /*********************************************************************/
        /* ARITHMETIC/OPERATOR METHODS */
        SparseVector<T> operator+(SparseVector<T>&);
        void operator+=(SparseVector<T>&);
        SparseVector<T> operator-(); /* unary minus: just negation of object */
        SparseVector<T> operator-(SparseVector<T>&);
        void operator-=(SparseVector<T>&);
        SparseVector<T> operator*(T);
        void operator*=(T);
        T operator[](int index) { return this->get_value(index); }
        bool operator==(SparseVector<T>&);

        SparseVector<T> conjugate(); /* Conjugate the vector and return the new structure */
        virtual void conjugate_in() { return; } /* Conjugate the vector inplace */
        
        /*********************************************************************/
        /* REPRESENTATION METHODS */
        string to_string() {
            stringstream output = stringstream("(");
            if (this->dim > 0) {
                output << this->coeff_to_string((*this)[0]);
                for (int i = 1; i < this->dim; i++) {
                    output << ", " << this->coeff_to_string((*this)[i]);
                }
            }
            output << ")";

            return output.str();
        }

        virtual string coeff_to_string(T element) { throw std::logic_error("Method 'coeff_to_string' not implemented"); }
};

class QQSparseVector : public SparseVector<rational<int>> {
    public:
        QQSparseVector(int dim) : SparseVector(dim) { };

        string coeff_to_string(rational<int> element) {
            string output = std::to_string(element.numerator());
            if (element.denominator() != 1) {
                output += "/" + std::to_string(element.denominator());
            }
            return output;
        }

        float norm() {
            rational<int> squared_norm = this->inner_product(*this);
            return sqrt(squared_norm.numerator() / (float) squared_norm.denominator());
        }
};

class CCSparseVector : public SparseVector<complex<double>> {
    public:
        // using SparseVector<complex<double>>::SparseVector;

        void conjugate_in() {
            for (pair<int,complex<double>> ppair : this->nonzero) {
                this->set_value(ppair.first, complex<double>(ppair.second.real(), -ppair.second.imag()));
            }
        }

        string coeff_to_string(complex<double> element) {
            if (element.real() == (double) 0 && element.imag() == (double) 0) {
                return "0";
            } else if (element.real() == (double) 0) {
                return std::to_string(element.imag()) + "*i";
            } else if (element.imag() == (double) 0) {
                return std::to_string(element.real());
            } else {
                return std::to_string(element.real()) + " + " + std::to_string(element.imag()) + "*i";
            }
        }

        float norm() {
            complex<double> squared_norm = this->inner_product(*this);
            return (float) sqrt(squared_norm.real());
        }
};

// template <typename T>
// class SparseRowMatrix {
//     private:

//     public:

//     protected:
// };

#endif