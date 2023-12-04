#include "Linalg.hpp"

#include "dd/ComplexValue.hpp"
#include "dd/Package.hpp"

CacheDDPackage* CacheDDPackage::singleton_ = nullptr;
CacheDDPackage* CacheDDPackage::GetInstance() {
    if(singleton_==nullptr) {
        singleton_ = new CacheDDPackage();
    }

    return singleton_;
}
dd::Package<>* CacheDDPackage::get_dd_package(luint nQbits) {
    if (this->cache.count(nQbits) == 0) {
        this->cache[nQbits] = new dd::Package<>(nQbits);
    }

    return this->cache[nQbits];
}
dd::Package<> * dd_package(luint nQbits) {
    return CacheDDPackage::GetInstance()->get_dd_package(nQbits);
}
/******************************************************************************************************************/
DDVector::DDVector(luint nQbits, unordered_map<dd::vEdge, CC> parts) : DDVector(nQbits) {
    for (std::pair<dd::vEdge,CC> part : parts) {
        this->components[part.first] = part.second;
    }
}

DDVector::DDVector(luint nQbits, vector<dd::vEdge> parts) : DDVector(nQbits) {
    for (dd::vEdge part : parts) {
        this->components[part] = CC(1);
    }
}

DDVector::DDVector(const DDVector& other) : DDVector(other.qbits) {
    for (std::pair<dd::vEdge, CC> pair : other.components) {
        this->components[pair.first] = pair.second;
    }
}

CC DDVector::inner_product(DDVector& vector) {
    CC result = CC(0);
    dd::Package<> * package = dd_package(this->qbits);
    for (std::pair<dd::vEdge,CC> this_part : this->components) {
        for (std::pair<dd::vEdge,CC> other_part : vector.components) {
            dd::ComplexValue value = package->innerProduct(this_part.first, other_part.first);
            result += (this_part.second*other_part.second)*CC(value.r, value.i);
        }
    }
    return result;
}

CC DDVector::inner_product(dd::vEdge& vector) {
    DDVector aux = DDVector(this->qbits, vector);
    return this->inner_product(aux);
}

/* Abstract methods */
double DDVector::norm() {
    CC sqrd_norm = this->inner_product(*this);

    return sqrt(sqrd_norm.real());
}
DDVector& DDVector::conjugate() {
    throw std::logic_error("Method 'conjugate' not yet implemented");
}
DDVector& DDVector::normalize() {
    static DDVector copy = (*this);
    copy.normalize_in();

    return copy;
}
void DDVector::normalize_in() {
    CC inv_norm = CC(1/this->norm());
    this->operator*=(inv_norm);
}

DDVector DDVector::apply_circuit(const dd::mEdge& circuit) {
    dd::Package<> * package = dd_package(this->qbits);

    unordered_map<dd::vEdge, CC> new_parts;
    for (std::pair<dd::vEdge, CC> pair : this->components) {
        dd::vEdge new_diagram = package->multiply<dd::mNode,dd::vNode>(circuit, pair.first, 0, false);
        CC new_value = pair.second;
        if (new_parts.count(new_diagram) > 0) {
            new_value += new_parts[new_diagram];
        }
        new_parts[new_diagram] = new_value;
    }

    return DDVector(this->qbits, new_parts);
}

DDVector DDVector::operator+(const DDVector& other) {
    DDVector result = (*this);
    result.operator+=(other);

    return result;
}
DDVector DDVector::operator-() {
    CC minus_one = CC(-1);
    return this->operator*(minus_one);
}
DDVector DDVector::operator-(const DDVector& other) {
    DDVector result = (*this);
    result.operator-=(other);

    return result;
}
DDVector DDVector::operator*(CC& to_scale) {
    DDVector result = (*this);
    result.operator*=(to_scale);

    return result;
}
void DDVector::operator+=(const DDVector& other) {
    for (std::pair<dd::vEdge, CC> other_part : other.components) {
        CC this_value = CC(0);
        if (this->components.count(other_part.first) != 0) {
            this_value = this->components.at(other_part.first);
        }
        this->components[other_part.first] = (this_value + other_part.second);
    }
}
void DDVector::operator-=(const DDVector& other) {
    for (std::pair<dd::vEdge, CC> other_part : other.components) {
        CC this_value = CC(0);
        if (this->components.count(other_part.first) != 0) {
            this_value = this->components.at(other_part.first);
        }
        this->components[other_part.first] = (this_value - other_part.second);
    }
}
void DDVector::operator*=(CC& to_scale) {
    for (std::pair<dd::vEdge, CC> this_part : this->components) {
        this->components[this_part.first] = this_part.second * to_scale;
    }
}

string DDVector::to_string() {
    stringstream stream;
    stream << "Vector of dimension " << static_cast<luint>(pow(2, this->qbits)) << " with " << this->components.size() << " components (total norm: " << this->norm() << ")" << endl;
    for (std::pair<dd::vEdge, CC> pair : this->components) {
        stream << "\t (" << CC_to_string(pair.second) << ") * [";
        for (CC coeff : pair.first.getVector()) {
            stream << CC_to_string(coeff) << ", ";
        }
        stream << "]";
    }
    return stream.str();
}

std::ostream& operator<<(ostream& stream, DDVector& vector) {
    stream << vector.to_string();
    return stream;
}

/******************************************************************************************************************/
/*********************************************************************/
/* ATTRIBUTE/PROPERTIES */
vector<double> DDSubspace::norms() {
    vector<double> result = vector<double>(this->dimension());

    for (luint i = 0; i < this->dimension(); i++) {
        result[i] = sqrt(this->basis[i].inner_product(this->basis[i]).real());
    }
    
    return result;
}

/*********************************************************************/
/* GETTING/SETTING DATA METHODS */
void DDSubspace::reduce_vector(DDVector* vector) {
    /* Method that reduced a vector according to 'this' in-place. */
    /* This method does a MGS reduction of a vector, becoming numericaly stable*/
    cout << "--------------------- Reducing a vector: " << vector->norm();
    for (luint i = 0; i < this->dimension(); i++) {
        DDVector bas = this->basis[i]; // Copy of the basis vector
        CC to_scale = bas.inner_product(*vector); // Inner product
        bas.operator*=(to_scale); // We scale
        vector->operator-=(bas); // We remove
    }
    cout << " --> " << vector->norm() << endl;
}

bool DDSubspace::contains(DDVector& vector) {
    /* Returns whether a vector is in the space or not */
    return vector.norm() < this->max_error;
}

CCSparseVector DDSubspace::find_in(DDVector& vector) {
    /* Returns a vector representing the coordinates of vector in "this"*/
    if (! this->contains(vector)) { throw std::logic_error("The vector is not in the space"); }
    CCSparseVector result = CCSparseVector(this->dimension());
    for (luint i = 0; i < this->dimension(); i++) {
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
bool DDSubspace::absorb_new_vector(DDVector& vector) {
    this->reduce_vector(&vector); // We reduce the vector
    if (! this->contains(vector)) {
        vector.normalize_in();
        this->basis.push_back(vector);
        return true;
    } else {
        cout << "Found an element inside: " << vector.norm() << endl;
    }
    return false;
}

/*********************************************************************/
/* COMPUTATIONAL METHODS */
/** 
 * Method that computes the minimal invariant subspace under some matrices.
 * 
 * Given a vector space `V` (given by `this`) and a set of matrices `(M_i)_i`, this method computes the minimal 
 * extension of `V` (say `W`) that is invariant under every matrix `M_i`, i.e., for every `v \in W`, `M_i v \in W`.
 * 
 * A matrix is a vector of vectors. In this case, we will use a ``std::vector`` of ``DDVector`` as a row 
 * representation of the matrix. Hence, the multiplication `M_i v` can be done via inner products of the rows of `M_i`
 * with the vector `v` (which will also be defined with a DDVector).
 * 
 * This method extends in-place the current subspace we are working with. This method returns the new dimension of the 
 * vector space.
*/
luint DDSubspace::minimal_invariant_space(vector<dd::mEdge>& circuits) {
    // We create a queue with the first round we need to check
    queue<DDVector> to_process;

    for (DDVector current : this->basis) {
        for (dd::mEdge circuit : circuits) {
            // We do multiplication matrix*current
            DDVector result = current.apply_circuit(circuit);
            // We add this vector to the queue
            to_process.push(result);
        }
    }

    //We now iterate on the queue until this is empty
    while ((!to_process.empty()) && (this->dimension() < this->ambient_dimension())) {
        cout << "Starting computation of minimal invariant with " << to_process.size() << " vectors. " << endl;
        DDVector current = to_process.front(); to_process.pop(); // We take the first element
        
        cout << "Current vector: " << current.norm() << endl;
        bool absorbed = this->absorb_new_vector(current);
        cout << "Was absorbed: " << absorbed << endl;
        if (absorbed) { // We have increased the dimension, we need to add new vectors
            for (dd::mEdge circuit : circuits) {
                // We do multiplication matrix*current
                DDVector result = current.apply_circuit(circuit);
                // We add this vector to the queue
                to_process.push(result);
            } 
        }
    }

    return this->dimension();
}

dd::CMat DDSubspace::reduced_matrix(dd::mEdge& circuit) {
    vector<DDVector> ULd;
    for (luint i = 0; i < this->dimension(); i++) {
        ULd.push_back(this->basis[i].apply_circuit(circuit));
    }

    dd::CMat result = dd::CMat(this->dimension());
    for (luint i = 0; i < this->dimension(); i++) {
        result[i] = dd::CVec(this->dimension());
        for (luint j = 0; j < this->dimension(); j++) {
            result[i][j] = this->basis[i].inner_product(ULd[j]); // TODO: Careful with the conjugation: is it there?
        }
    }

    return result;
}