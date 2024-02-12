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
template <typename V, typename M, typename C>
vector<double> Subspace<V,M,C>::norms() {
    vector<double> result = vector<double>(this->dimension());

    for (luint i = 0; i < this->dimension(); i++) {
        result[i] = this->norm(this->basis[i]);
    }

    return result;
}

/*********************************************************************/
/* GETTING/SETTING DATA METHODS */
template <typename V, typename M, typename C>
V* Subspace<V,M,C>::reduce_vector(V* vector) {
    V* result = this->scale(vector, this->coeff(1.)); // We copy the input vector

    // cout <<"\t\tReducing vector using " << this->dimension() << " vectors:" << endl;
    
    /* This method does a MGS reduction of a vector, becoming numericaly stable*/
    for (luint i = 0; i < this->dimension(); i++) {
        // cout <<"\t\t\t[R] " << this->print_vector(this->basis[i]) << endl;
        C coeff = this->inner_product(this->basis[i], result);
        // cout <<"\t\t\t    (c): " << CC_to_string(coeff) << endl;
        V* to_rem = this->scale(this->basis[i], coeff*this->coeff(-1.));
        // cout <<"\t\t\t    ( ): " << this->print_vector(to_rem) << endl;
        V* aux = this->add(result, to_rem);
        this->free_vector(result); result = aux;
        // cout <<"\t\t\t    " << this->print_vector(result) << endl;
    }
    // cout <<"\t\t\t" << this->print_vector(vector) << "\n\t\t\t->\n\t\t\t" << this->print_vector(result) << endl;

    return result;
}

template <typename V, typename M, typename C>
bool Subspace<V,M,C>::contains(V* vector) {
    V* reduced = this->reduce_vector(vector);
    bool result = this->norm(reduced) <= this->max_error;
    this->free_vector(reduced); //Removing memory for reduced vector
    return result;
}

template <typename V, typename M, typename C>
bool Subspace<V,M,C>::absorb_new_vector(V* vector) {
    V* reduced = this->reduce_vector(vector);
    bool result = false;
    double norm_reduced = this->norm(reduced);
    // cerr << endl <<"\t\tNorm: " << this->norm(reduced) << " -- " << norm_reduced << "[" << this->max_error << "]";
    if (norm_reduced > this->max_error) {
        V* to_add = this->scale(reduced, this->coeff(1/this->norm(reduced)));
        this->basis.push_back(to_add);
        result = true;
    } 

    this->free_vector(reduced); //Removing memory for reduced vector
    return result;
}

/*********************************************************************/
/* COMPUTATIONAL METHODS */
template <typename V, typename M, typename C>
luint Subspace<V,M,C>::minimal_invariant_space(vector<M>& matrices) {
    // We create a queue with the first round we need to check
    queue<V*> to_process;

    for (V* current : this->basis) {
        // cout <<"Element of the basis: " << endl << this->print_vector(current) << endl;
        // cout <<"Generated:" << endl;
        for (M matrix : matrices) {
            // We do multiplication matrix*current
            V* result = this->apply(current, matrix);
            // cout <<"\t[-] " << this->print_vector(result) << endl; 
            // We add this vector to the queue
            to_process.push(result);
        }
    }

    //We now iterate on the queue until this is empty
    while ((!to_process.empty()) && (this->dimension() < this->ambient_dimension())) {
        V* current = to_process.front(); to_process.pop(); // We take the first element
        
        // cout <<this->print_vector(current) << endl;
        bool absorbed = this->absorb_new_vector(current);
        // cout <<"Was absorbed: " << absorbed << endl;
        // We release the memory for the processed vector
        this->free_vector(current);
        if (absorbed) { // We have increased the dimension, we need to add new vectors
            current = this->basis[this->dimension()-1];
            // cout <<"Added vector:\n\t[+] " << this->print_vector(current) << endl;
            // cout <<"Generated:" << endl;
            for (M matrix : matrices) {
                // We do multiplication matrix*current
                // cerr << endl << "\tStarting application";
                // cerr << endl << "\tSizes in basis: (";
                // for (luint i = 0; i < this->dimension(); i++) { cerr << this->vector_dim(this->basis[i]) << ", "; }
                // cerr << ")";
                V* result = this->apply(current, matrix);
                // cout <<"\t[-] " << this->print_vector(result) << endl; 
                // We add this vector to the queue
                to_process.push(result);
            } 
        }
    }
    // We release the memory if any vector remains on the queue
    while (!to_process.empty()) {
        V* current = to_process.front(); to_process.pop();
        this->free_vector(current);
    }

    // We return the new dimension of the subspace
    return this->dimension();
}

template <typename V, typename M, typename C>
vector<V> Subspace<V,M,C>::lumping_matrix() {
    vector<V> result = vector<V>();
    for (luint i = 0; i < this->dimension(); i++) {
        result.push_back(*this->basis[i]);
    }
    return result;
}

template <typename V, typename M, typename C>
vector<vector<C>> Subspace<V,M,C>::reduced_matrix(M& matrix) {
    // We create the strucutre for the matrix
    vector<vector<C>> result = vector<vector<C>>(this->dimension());
    vector<V*> conjugates = vector<V*>();
    for (luint i = 0; i < this->dimension(); i++) {
        V* conj = this->conjugate(this->basis[i]);
        conjugates.push_back(conj);
    }

    for (luint j = 0; j < this->dimension(); j++) {
        // V* conj = this->conjugate(this->basis[j]);
        V* Uld = this->apply(conjugates[j], matrix);
        // V* Uld_conj = this->conjugate(Uld);
        for (luint i = 0; i < this->dimension(); i++) {
            // C to_put = this->inner_product(this->basis[i], Uld_conj);
            C to_put = this->inner_product(Uld, conjugates[i]);
            result[i].push_back(to_put);
        }
        // Freeing memory if needed
        this->free_vector(Uld); 
        // this->free_vector(Uld_conj);
        // this->free_vector(conj);
    }
    for (V* v : conjugates) { this->free_vector(v); }

    return result;
}

template class Subspace<CCSparseVector, vector<CCSparseVector>, CC>;
template class Subspace<dd::vEdge, qc::QuantumComputation, dd::ComplexValue>;

/*******************************************************************************************************************
 * 
 * CLASS FOR CC-SUBSPACE
 * 
********************************************************************************************************************/

double CCSubspace::norm(CCSparseVector* vector) {
    return vector->norm();
}
CC CCSubspace::coeff(double coeff) {
    return CC(coeff);
}
CCSparseVector* CCSubspace::apply(CCSparseVector* u, vector<CCSparseVector>& M) {
    // Computes M*u so the result is nrows(M) x 1
    CCSparseVector * result = new CCSparseVector(M.size());
    CCSparseVector conj = u->conjugate(); // We conjugate so when the inner_product occurs, we do M*u
    for (luint i = 0; i < M.size(); i++) {
        result->set_value(i, M[i].inner_product(conj));
    }

    return result;
}
CCSparseVector* CCSubspace::scale(CCSparseVector* u, CC c) {
    CCSparseVector * result = new CCSparseVector(u->dimension());

    std::unordered_set<luint>::iterator it = u->nonzero_iterator();
    while (it != u->nonzero_iterator_end()) {
        result->set_value(*it, c*u->get_value(*it));
        it++;
    } 

    return result;
}
CCSparseVector* CCSubspace::add(CCSparseVector* u, CCSparseVector* v) {
    CCSparseVector *result = new CCSparseVector(u->dimension());
    
    *result = u->operator+(*v); // Using the copy operator on the sum of two vectors

    return result;
}
CC CCSubspace::inner_product(CCSparseVector* u, CCSparseVector*v) {
    return u->inner_product(*v);
}
CCSparseVector* CCSubspace::conjugate(CCSparseVector* u) {
    CCSparseVector* result = new CCSparseVector(u->dimension());
    *result = u->conjugate();

    return result;
}
