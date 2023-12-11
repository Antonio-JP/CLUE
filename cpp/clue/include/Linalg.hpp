#ifndef CLUE_LA_H
#define CLUE_LA_H

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <sstream>
#include <string>
#include "Types.hpp"
#include "QuantumComputation.hpp"
#include "dd/Node.hpp"
#include "dd/Package.hpp"

using namespace std;
using namespace clue;

typedef long unsigned int luint;

class CacheDDPackage {
    protected:
        /* Private constructor */
        CacheDDPackage() = default;
        static CacheDDPackage* singleton_;
        unordered_map<luint, dd::Package<>*> cache;
    public:
        /* Removing cloning and assignment methods */
        CacheDDPackage(CacheDDPackage &) = delete;
        CacheDDPackage& operator=(const CacheDDPackage &) = delete;
        ~CacheDDPackage() = default;

        /* Method for getting instance */
        static CacheDDPackage* GetInstance();

        /* Logic of the Cache */
        dd::Package<>* get_dd_package(luint nQbits);
};
dd::Package<> * dd_package(luint nQbits);

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
        unordered_set<luint>::iterator nonzero_iterator_end() { return this->nonzero_indices.end(); }
        double density() { return static_cast<double>(this->nonzero_count()) / static_cast<double>(this->dim); }

        bool is_zero() { return this->nonzero_count() == 0; }
        
        void reduce(T coeff, SparseVector<T>& other); // method that computes inplace this - coeff*other
        T inner_product(SparseVector<T>&);

        /* Abstract methods */
        virtual double norm() = 0;
        virtual SparseVector<T>* normalize() = 0;
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
        QQSparseVector* normalize();
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
        CCSparseVector* normalize();
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
class DDVector {
    private:
        luint qbits;
        unordered_map<dd::vEdge, CC> components;

    public:
        DDVector(luint nQbits) { this->qbits = nQbits; }
        DDVector(luint nQbits, unordered_map<dd::vEdge, CC> parts);
        DDVector(luint nQbits, vector<dd::vEdge> parts);
        DDVector(luint nQbits, dd::vEdge& part) : DDVector(nQbits, vector<dd::vEdge>({part})) { }
        DDVector(const DDVector &);
        ~DDVector() = default;

        CC inner_product(DDVector&);
        CC inner_product(dd::vEdge& vector);

        luint nQbits() {return this->qbits; }
        double norm();
        DDVector* conjugate();
        void conjugate_in();
        DDVector* normalize();
        void normalize_in();

        DDVector apply_circuit(const dd::mEdge&);

        DDVector operator+(const DDVector& other);
        DDVector operator-();
        DDVector operator-(const DDVector& other);
        DDVector operator*(CC& to_scale);
        void operator+=(const DDVector& other);
        void operator-=(const DDVector& other);
        void operator*=(CC& to_scale);

        string to_string();
        friend std::ostream& operator<<(std::ostream&, DDVector&);
};

/*************************************************************************/
/* Methods and functions for vectors/matrices */
string CC_to_string(CC& number);
string CC_to_string(dd::ComplexValue& number);
string CC_to_string(dd::Complex& number);

string vector_to_string(dd::CVec& vector);
string vector_to_string(CCSparseVector& vector);
string vector_to_string(vector<dd::ComplexValue>& vector);
string vector_to_string(vector<dd::Complex>& vector);
string vector_to_string(DDVector& vector);
string vector_to_string(dd::vEdge& vector);

string matrix_to_string(dd::CMat&);
string matrix_to_string(vector<CCSparseVector>&);
string matrix_to_string(vector<vector<dd::ComplexValue>>&);
string matrix_to_string(vector<vector<dd::Complex>>&);

bool is_diagonal(dd::CMat&);
bool is_diagonal(vector<CCSparseVector>&);

bool is_square(dd::CMat&);
bool is_square(vector<CCSparseVector>&);
dd::CMat sparse_to_dense(vector<CCSparseVector>&);
dd::CMat ComplexValue_to_complex(vector<vector<dd::ComplexValue>>&);

dd::CMat matmul(dd::CMat&, dd::CMat&);
dd::CMat matmul(vector<CCSparseVector>&, dd::CMat&);
dd::CMat matmul(dd::CMat&, vector<CCSparseVector>&);
dd::CMat matmul(vector<CCSparseVector>&, vector<CCSparseVector>&);

dd::CMat identity_matrix(luint);
dd::CMat matrix_power(dd::CMat&, luint);
dd::CMat matrix_power(vector<CCSparseVector>&, luint);


/*************************************************************************/
/* Class for Subspace */
template <typename V, typename M, typename C>
class Subspace {
    protected:
        luint dim; // Ambient dimension of the vector space
        double max_error; // Maximal error allow for vectors to be in the space (0 indicates exact algorithms)

        /* Operations required for public methods */
        virtual double norm(V*) = 0; // Compute the norm of a vector from its pointer
        virtual C coeff(double) = 0; // Compute the norm of a vector from its pointer
        virtual V* apply(V*, M&) = 0; // Compute the application of M to V (i.e., M*V)
        virtual V* scale(V*, C) = 0; // Scales a vector using a complex number
        virtual V* add(V*, V*) = 0; // Compute the addition of two vectors
        virtual C inner_product(V*, V*) = 0; // Compute the inner product of two vectors
        virtual V* conjugate(V*) = 0; // Conjugate a vector
        virtual void free_vector(V* vector) { delete vector; }
        virtual string print_vector(V* vector) = 0;
        virtual ~Subspace() { for (V* vector : this->basis) {delete vector; } }

    public:
        vector<V*> basis; // Basis of the subspace
        
        /*********************************************************************/
        /* Constructors for Subspaces */
        Subspace(luint ambient_dimension, double error) { this->dim = ambient_dimension; this->max_error = error; }
        Subspace(luint ambient_dimension) : Subspace<V,M,C>(ambient_dimension, 1e-6) { }
        

        /*********************************************************************/
        /* ATTRIBUTE/PROPERTIES */
        luint ambient_dimension() { return this->dim; }
        luint dimension() { return this->basis.size(); }
        vector<double> norms();
        
        /*********************************************************************/
        /* GETTING/SETTING DATA METHODS */
        V* reduce_vector(V*); // Takes the pointer to a vector and return the reduction of the vector w.r.t. this
        bool contains(V*); // Check whether a vector is in this or not
        bool absorb_new_vector(V*); // Absorbs a new vector to the space, leaving it still if the vector is in.
        
        /*********************************************************************/
        /* COMPUTATIONAL METHODS */
        luint minimal_invariant_space(vector<M>&);
        vector<V> lumping_matrix();
        vector<vector<C>> reduced_matrix(M&);
};

class CCSubspace : public Subspace<CCSparseVector, vector<CCSparseVector>, CC> {
    protected:
        double norm(CCSparseVector*); // Compute the norm of a vector from its pointer
        CC coeff(double); // Compute the norm of a vector from its pointer
        CCSparseVector* apply(CCSparseVector*, vector<CCSparseVector>&); // Compute the application of M to V (i.e., M*V)
        CCSparseVector* scale(CCSparseVector*, CC); // Scales a vector using a complex number
        CCSparseVector* add(CCSparseVector*, CCSparseVector*); // Compute the addition of two vectors
        CC inner_product(CCSparseVector*, CCSparseVector*); // Compute the inner product of two vectors
        CCSparseVector* conjugate(CCSparseVector*); // Conjugate a vector
        string print_vector(CCSparseVector* vector) { return vector->to_string(); }
    
    public:
        using Subspace<CCSparseVector, vector<CCSparseVector>, CC>::Subspace;
};

class DDSubspace : public Subspace<DDVector, dd::mEdge, dd::ComplexValue> {
    protected:
        double norm(DDVector*); // Compute the norm of a vector from its pointer
        dd::ComplexValue coeff(double); // Compute the norm of a vector from its pointer
        DDVector* apply(DDVector*, dd::mEdge&); // Compute the application of M to V (i.e., M*V)
        DDVector* scale(DDVector*, dd::ComplexValue); // Scales a vector using a complex number
        DDVector* add(DDVector*, DDVector*); // Compute the addition of two vectors
        dd::ComplexValue inner_product(DDVector*, DDVector*); // Compute the inner product of two vectors
        DDVector* conjugate(DDVector*); // Conjugate a vector
        string print_vector(DDVector* vector) { return vector->to_string(); }
    
    public:
        using Subspace<DDVector, dd::mEdge, dd::ComplexValue>::Subspace;
};

class FullDDSubspace : public Subspace<dd::vEdge, qc::QuantumComputation, dd::ComplexValue> {
    protected:
        luint nQbits;
        std::unique_ptr<dd::Package<>> package;

        double norm(dd::vEdge*); // Compute the norm of a vector from its pointer
        dd::ComplexValue coeff(double); // Compute the norm of a vector from its pointer
        dd::vEdge* apply(dd::vEdge*, qc::QuantumComputation&); // Compute the application of M to V (i.e., M*V)
        dd::vEdge* scale(dd::vEdge*, dd::ComplexValue); // Scales a vector using a complex number
        dd::vEdge* add(dd::vEdge*, dd::vEdge*); // Compute the addition of two vectors
        dd::ComplexValue inner_product(dd::vEdge*, dd::vEdge*); // Compute the inner product of two vectors
        dd::vEdge* conjugate(dd::vEdge*); // Conjugate a vector
        void free_vector(dd::vEdge*) { return; }
        string print_vector(dd::vEdge* vector) { dd::CVec v = vector->getVector(); return vector_to_string(v); }
    
    public:
        FullDDSubspace(luint qbits, double error) : 
            Subspace<dd::vEdge, qc::QuantumComputation, dd::ComplexValue>(static_cast<luint>(pow(2, qbits)), error) { 
                this->nQbits = qbits; 
                this->package = std::unique_ptr<dd::Package<>>(dd_package(qbits));
        }
        FullDDSubspace(luint qbits) : FullDDSubspace(qbits, 1e-6) { }
};

#endif