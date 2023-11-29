#ifndef CLUE_LA_H
#define CLUE_LA_H

#include <map>
#include <unordered_set>
#include <vector>
#include <sstream>
#include <string>
#include "Types.hpp"
#include "CircuitSimulator.hpp"

using namespace std;
using namespace clue;

typedef long unsigned int luint;

/*************************************************************************/
/* Class for vector */
template <typename T>
class SparseVector {
    private:
        luint dim;
        unordered_set<luint> nonzero_indices;
    protected:
        map<luint, T> nonzero;
    public:
        /*********************************************************************/
        /* CONSTRUCTORS */
        SparseVector(luint dim);
        SparseVector(vector<T> dense_vector) : SparseVector(dense_vector.size()) {
            for (luint i = 0; i < this->dim; i++) {
                if (dense_vector[i] != T(0)) {
                    this->set_value(i, dense_vector[i]);
                }
            }
        }
        SparseVector(const SparseVector<T>& sparse_vector) : SparseVector(sparse_vector.dim) {
            for (pair<luint, T> ppair : sparse_vector.nonzero) {
                this->set_value(ppair.first, ppair.second);
            }
        }

        virtual SparseVector<T>& operator=(const SparseVector<T>&) = default;
        virtual ~SparseVector() = default;
        
        /*********************************************************************/
        /* ATTRIBUTE/PROPERTIES */
        luint dimension() { return this->dim; }
        luint nonzero_count() { return this->nonzero.size(); }
        luint first_nonzero() { return this->nonzero.begin()->first; }
        unordered_set<luint>::iterator nonzero_iterator() { return this->nonzero_indices.begin(); }
        double density() { return static_cast<double>(this->nonzero_count()) / static_cast<double>(this->dim); }

        bool is_zero() { return this->nonzero_count() == 0; }
        
        void reduce(T coeff, SparseVector<T>& other); // method that computes inplace this - coeff*other
        T inner_product(SparseVector<T>&);

        /* Abstract methods */
        virtual double norm() = 0;
        virtual SparseVector<T>& normalize() = 0;
        virtual void normalize_in() = 0;

        /*********************************************************************/
        /* GETTING/SETTING DATA METHODS */
        T get_value(luint index);
        void set_value(luint index, T value);

        vector<T> to_list();

        /*********************************************************************/
        /* ARITHMETIC/OPERATOR METHODS */
        void operator+=(SparseVector<T>&);
        void operator-=(SparseVector<T>&);
        void operator*=(T);
        T operator[](luint index) { return this->get_value(index); }
        bool operator==(SparseVector<T>&);
        bool operator!=(SparseVector<T>& other) { return !((*this) == other); }
        void conjugate_in(); /* Conjugate the vector inplace */

        virtual T conjugate_coeff(T coeff) = 0;
        /* Abstract methods */
        
        /*********************************************************************/
        /* REPRESENTATION METHODS */
        std::string to_string() {
            stringstream output; output << "(";
            if (this->dim > 0) {
                output << this->coeff_to_string((*this)[0]);
                for (luint i = 1; i < this->dim; i++) {
                    output << ", " << this->coeff_to_string((*this)[i]);
                }
            }
            output << ")";

            return output.str();
        }

        /* Abstract methods */
        virtual std::string coeff_to_string(T element) = 0;
};


class QQSparseVector : public SparseVector<QQ> {
    public:
        using SparseVector<QQ>::SparseVector;
        using SparseVector<QQ>::operator=;

        /* VIRTUAL METHODS */
        double norm();
        QQSparseVector& normalize();
        void normalize_in();
        QQ conjugate_coeff(QQ coeff) { return coeff; }
        std::string coeff_to_string(QQ element);

        /* ARITHMETIC METHODS THAT RETURN AN OBJECT */
        QQSparseVector operator+(QQSparseVector&);
        QQSparseVector operator-(); /* unary minus: just negation of object */
        QQSparseVector operator-(QQSparseVector&); 
        QQSparseVector operator*(QQ); 
        QQSparseVector conjugate() { return (*this); } // No conjugation in a rational vector
        // QQSparseVector& operator=(const QQSparseVector&) = default;
};


class CCSparseVector : public SparseVector<CC> {
    public:
        using SparseVector<CC>::SparseVector;
        using SparseVector<CC>::operator=;
        
        /* VIRTUAL METHODS */
        double norm();
        CCSparseVector& normalize();
        void normalize_in();
        CC conjugate_coeff(CC coeff) { return CC(coeff.real(), -coeff.imag()); }
        std::string coeff_to_string(CC element);
        
        /* ARITHMETIC METHODS THAT RETURN AN OBJECT */
        CCSparseVector operator+(CCSparseVector&);
        CCSparseVector operator-(); /* unary minus: just negation of object */
        CCSparseVector operator-(CCSparseVector&); 
        CCSparseVector operator*(CC); 
        CCSparseVector conjugate();
        // CCSparseVector& operator=(const CCSparseVector&) = default;
};


/*************************************************************************/


/*************************************************************************/
/* Class for Subspace */
class CCSubspace {
    private:
        luint dim;
        double max_error;
    public:
        vector<CCSparseVector> basis; // Temporary: move to private
        CCSubspace(luint ambient_dimension, double error) { this->dim = ambient_dimension; this->max_error = error; }
        CCSubspace(luint ambient_dimension) : CCSubspace(ambient_dimension, 1e-6) { }

        /*********************************************************************/
        /* ATTRIBUTE/PROPERTIES */
        luint ambient_dimension() { return this->dim; }
        luint dimension() { return this->basis.size(); }
        vector<double> densities();
        vector<double> norms();

        /*********************************************************************/
        /* GETTING/SETTING DATA METHODS */
        void reduce_vector(CCSparseVector*);
        bool contains(CCSparseVector&);
        CCSparseVector find_in(CCSparseVector&);
        bool absorb_new_vector(CCSparseVector&);

        /*********************************************************************/
        /* COMPUTATIONAL METHODS */
        luint minimal_invariant_space(vector<vector<CCSparseVector>>& matrices);
};

class DDSubspace {
    private:
        luint dim;
        double max_error;
    public:
        vector<CCSparseVector> basis; // Temporary: move to private
        DDSubspace(luint ambient_dimension, double error) { this->dim = ambient_dimension; this->max_error = error; }
        DDSubspace(luint ambient_dimension) : DDSubspace(ambient_dimension, 1e-6) { }

        /*********************************************************************/
        /* ATTRIBUTE/PROPERTIES */
        luint ambient_dimension() { return this->dim; }
        luint dimension() { return this->basis.size(); }
        vector<double> norms();

        /*********************************************************************/
        /* GETTING/SETTING DATA METHODS */
        void reduce_vector(CCSparseVector*);
        bool contains(CCSparseVector&);
        CCSparseVector find_in(CCSparseVector&);
        bool absorb_new_vector(CCSparseVector&);

        /*********************************************************************/
        /* COMPUTATIONAL METHODS */
        luint minimal_invariant_space(vector<vector<CCSparseVector>>& matrices);
};

#endif