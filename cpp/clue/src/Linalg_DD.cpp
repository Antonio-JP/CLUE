#include "Linalg.hpp"

#include <memory>

#include "dd/ComplexValue.hpp"
#include "dd/Package.hpp"
#include "dd/Simulation.hpp"

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

string vector_to_string(dd::CVec& vector) {
    stringstream stream;
    stream << "[";
    if (vector.size()) {
        stream << CC_to_string(vector[0]);
        for (luint i = 1; i < vector.size(); i++) {
            stream << ", " << CC_to_string(vector[i]);
        }
    }
    stream << "]";

    return stream.str();
}
string vector_to_string(vector<dd::ComplexValue>& vector) {
    stringstream stream;
    stream << "[";
    if (vector.size()) {
        stream << CC_to_string(vector[0]);
        for (luint i = 1; i < vector.size(); i++) {
            stream << ", " << CC_to_string(vector[i]);
        }
    }
    stream << "]";

    return stream.str();
}
string vector_to_string(vector<dd::Complex>& vector) {
    stringstream stream;
    stream << "[";
    if (vector.size()) {
        stream << CC_to_string(vector[0]);
        for (luint i = 1; i < vector.size(); i++) {
            stream << ", " << CC_to_string(vector[i]);
        }
    }
    stream << "]";

    return stream.str();
}

string matrix_to_string(dd::CMat& matrix) {
    stringstream stream;
    stream << "[" << endl;
    for (dd::CVec row : matrix) {
        stream << "\t" << vector_to_string(row) << endl;
    }
    stream << "]" << endl;
    return stream.str();
}
string matrix_to_string(vector<vector<dd::ComplexValue>>& matrix) {
    stringstream stream;
    stream << "[" << endl;
    for (vector<dd::ComplexValue> row : matrix) {
        stream << "\t" << vector_to_string(row) << endl;
    }
    stream << "]" << endl;
    return stream.str();
}
string matrix_to_string(vector<vector<dd::Complex>>& matrix) {
    stringstream stream;
    stream << "[" << endl;
    for (vector<dd::Complex> row : matrix) {
        stream << "\t" << vector_to_string(row) << endl;
    }
    stream << "]" << endl;
    return stream.str();
}

dd::vEdge& conjugate_edge(dd::vEdge& v, unique_ptr<dd::Package<>>& package) {
    return v;
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
    std::unique_ptr<dd::Package<>> package(dd_package(this->qbits));

    for (std::pair<dd::vEdge,CC> this_part : this->components) {
        for (std::pair<dd::vEdge,CC> other_part : vector.components) {
            dd::ComplexValue value = package->innerProduct(this_part.first, conjugate_edge(other_part.first, package));
            result += (this_part.second*CC(other_part.second.real(), -other_part.second.imag()))*CC(value.r, value.i);
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
    static DDVector copy = (*this);
    copy.conjugate_in();

    return copy;
}
void DDVector::conjugate_in() {
    std::unordered_map<dd::vEdge, CC> new_components;
    std::unique_ptr<dd::Package<>> package(dd_package(this->nQbits()));
    for (std::pair<dd::vEdge, CC> ppair : this->components) {
        dd::vEdge conj_edge = conjugate_edge(ppair.first, package);
        new_components[conj_edge] = CC(ppair.second.real(), -ppair.second.imag());
    }
    this->components.clear();
    for (std::pair<dd::vEdge, CC> ppair : new_components) {
        this->components[ppair.first] = ppair.second;
    }
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
// DDSubspace
// Virtual methods
double DDSubspace::norm(DDVector* vector) {
    return vector->norm();
}
dd::ComplexValue DDSubspace::coeff(double c) {
    return dd::ComplexValue(c);
}
DDVector* DDSubspace::apply(DDVector* v, dd::mEdge& M) {
    DDVector* result = new DDVector(v->apply_circuit(M));

    return result;
}
DDVector* DDSubspace::scale(DDVector* v, dd::ComplexValue c) {
    DDVector* result = new DDVector(*v);
    CC coeff = CC(c.r, c.i);
    result->operator*=(coeff); // We do in-place multiplication

    return result; // Return the new vector
}
DDVector* DDSubspace::add(DDVector* u, DDVector* v) {
    DDVector* result = new DDVector(*u);
    result->operator+=(*v); // We apply the in-place operation with the second vector
    return result;
}
dd::ComplexValue DDSubspace::inner_product(DDVector* u, DDVector* v) {
    CC value = u->inner_product(*v);
    return dd::ComplexValue(value.real(), value.imag());
}
DDVector* DDSubspace::conjugate(DDVector* v) {
    DDVector* result = new DDVector(*v); // We copy the vector
    result->conjugate_in();

    return result;
}

/******************************************************************************************************************/
// FullDDSubspace
double FullDDSubspace::norm(dd::vEdge* v) {
    return sqrt(this->package->innerProduct(*v, *v).r);
}
dd::ComplexValue FullDDSubspace::coeff(double c) {
    return dd::ComplexValue(c);
}
dd::vEdge* FullDDSubspace::apply(dd::vEdge* v, qc::QuantumComputation& M) {
    static dd::vEdge output = dd::simulate<>(&M, *v, this->package);
    return &output;
}
dd::vEdge* FullDDSubspace::scale(dd::vEdge* v, dd::ComplexValue c) {
    dd::Complex new_value = this->package->cn.lookup(c * dd::ComplexValue(v->w));
    dd::vEdge* output = new dd::vEdge{v->p, new_value}; // Creating the edge with the new weight
    return output;
}
dd::vEdge* FullDDSubspace::add(dd::vEdge* u, dd::vEdge* v) {
    static dd::vEdge result = this->package->add(*u, *v);
    return &result;
}
dd::ComplexValue FullDDSubspace::inner_product(dd::vEdge* u, dd::vEdge* v) {
    dd::vEdge v_conj = conjugate_edge(*v, this->package);
    return this->package->innerProduct(*u, v_conj);
}
dd::vEdge* FullDDSubspace::conjugate(dd::vEdge* v) {
    static dd::vEdge conj = conjugate_edge(*v, this->package);
    return &conj;
}