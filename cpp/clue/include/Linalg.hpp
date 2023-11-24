#ifndef CLUE_LA_H
#define CLUE_LA_H

#include <map>
#include <unordered_set>
#include <vector>
#include <complex>
#include <sstream>
#include <string>
#include <boost/rational.hpp>


using namespace std;
using namespace boost;

/* Class for vector */
template <typename T>
class SparseVector {
    private:
        int dim;
        unordered_set<int> nonzero_indices;
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
        SparseVector(const SparseVector<T>& sparse_vector) : SparseVector(sparse_vector.dim) {
            for (pair<int, T> ppair : sparse_vector.nonzero) {
                this->set_value(ppair.first, ppair.second);
            }
        }
        
        /*********************************************************************/
        /* ATTRIBUTE/PROPERTIES */
        int dimension() { return this->dim; }
        int nonzero_count() { return this->nonzero.size(); }
        int first_nonzero() { return this->nonzero.begin()->first; }
        unordered_set<int>::const_iterator nonzero_iterator() { return this->nonzero_indices.begin(); }
        float density() { return this->nonzero_count() / (float)this->dim; }

        bool is_zero() { return this->nonzero_count() == 0; }
        
        void reduce(T coeff, SparseVector<T>& other); // method that computes inplace this - coeff*other
        T inner_product(SparseVector<T>&);

        /* Abstract methods */
        virtual float norm() = 0;

        /*********************************************************************/
        /* GETTING/SETTING DATA METHODS */
        T get_value(int index);
        void set_value(int index, T value);

        vector<T> to_list();

        /*********************************************************************/
        /* ARITHMETIC/OPERATOR METHODS */
        void operator+=(SparseVector<T>&);
        void operator-=(SparseVector<T>&);
        void operator*=(T);
        T operator[](int index) { return this->get_value(index); }
        bool operator==(SparseVector<T>&);
        void conjugate_in(); /* Conjugate the vector inplace */

        virtual T conjugate_coeff(T coeff) = 0;
        /* Abstract methods */
        
        /*********************************************************************/
        /* REPRESENTATION METHODS */
        string to_string() {
            stringstream output; output << "(";
            if (this->dim > 0) {
                output << this->coeff_to_string((*this)[0]);
                for (int i = 1; i < this->dim; i++) {
                    output << ", " << this->coeff_to_string((*this)[i]);
                }
            }
            output << ")";

            return output.str();
        }

        /* Abstract methods */
        virtual string coeff_to_string(T element) = 0;
};

class QQSparseVector : public SparseVector<rational<int>> {
    public:
        using SparseVector<rational<int>>::SparseVector;

        /* VIRTUAL METHODS */
        float norm();
        rational<int> conjugate_coeff(rational<int> coeff) { return coeff; }
        string coeff_to_string(rational<int> element);

        /* ARITHMETIC METHODS THAT RETURN AN OBJECT */
        QQSparseVector operator+(QQSparseVector&);
        QQSparseVector operator-(); /* unary minus: just negation of object */
        QQSparseVector operator-(QQSparseVector&); 
        QQSparseVector operator*(rational<int>); 
        QQSparseVector conjugate() { return (*this); } // No conjugation in a rational vector
};

class CCSparseVector : public SparseVector<complex<double>> {
    public:
        using SparseVector<complex<double>>::SparseVector;
        
        /* VIRTUAL METHODS */
        float norm();
        complex<double> conjugate_coeff(complex<double> coeff) { return complex<double>(coeff.real(), -coeff.imag()); }
        string coeff_to_string(complex<double> element);
        
        /* ARITHMETIC METHODS THAT RETURN AN OBJECT */
        CCSparseVector operator+(CCSparseVector&);
        CCSparseVector operator-(); /* unary minus: just negation of object */
        CCSparseVector operator-(CCSparseVector&); 
        CCSparseVector operator*(complex<double>); 
        CCSparseVector conjugate();
};

#endif