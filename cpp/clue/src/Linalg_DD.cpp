#include "Linalg.hpp"

#include <memory>

#include "dd/ComplexValue.hpp"
#include "dd/Package.hpp"
#include "dd/Simulation.hpp"

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
string vector_to_string(CCSparseVector& vector) {
    stringstream stream;
    stream << "[";
    if (vector.dimension()) {
        CC value = vector[0];
        stream << CC_to_string(value);
        for (luint i = 1; i < vector.dimension(); i++) {
            value = vector[i];
            stream << ", " << CC_to_string(value);
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
string vector_to_string(DDVector& vector) {
    return vector.to_string();
}
string vector_to_string(dd::vEdge& vector) {
    dd::CVec to_print = vector.getVector();
    return vector_to_string(to_print);
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
string matrix_to_string(vector<CCSparseVector>& matrix) {
    stringstream stream;
    stream << "[" << endl;
    for (CCSparseVector row : matrix) {
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

bool is_diagonal(dd::CMat& matrix) {
    for (luint i = 0; i < matrix.size(); i++) {
        for (luint j = 0; j < matrix[i].size(); j++) {
            if (matrix[i][j] != CC(0) && i != j) {
                return false;
            }
        }
    }
    return true;
}
bool is_diagonal(vector<CCSparseVector>& matrix) {
    for (luint i = 0; i < matrix.size(); i++) {
        CCSparseVector& row = matrix[i];
        if (row.nonzero_count() > 1 || (row.nonzero_count() == 1 && row.first_nonzero() != i)) {
            return false;
        }
    }
    return true;
}
bool is_square(dd::CMat& A) {
    luint rows = A.size();

    for (luint i = 0; i < rows; i++) {
        if (A[i].size() != rows) {
            return false;
        }
    }
    return true;
}
bool is_square(vector<CCSparseVector>& A) {
    luint rows = A.size();

    for (luint i = 0; i < rows; i++) {
        if (A[i].dimension() != rows) {
            return false;
        }
    }
    return true;
}
dd::CMat sparse_to_dense(vector<CCSparseVector>& A) {
    dd::CMat result = dd::CMat(A.size());
    for (luint i = 0; i < A.size(); i++) {
        result[i] = dd::CVec(A[i].dimension());
        std::unordered_set<luint>::iterator it = A[i].nonzero_iterator();
        while (it != A[i].nonzero_iterator_end()) {
            result[i][*it] = A[i][*it];
            it++;
        }
    }
    return result;
}
dd::CMat ComplexValue_to_complex(vector<vector<dd::ComplexValue>>& A) {
    dd::CMat result = dd::CMat(A.size());
    for (luint i = 0; i < A.size(); i++) {
        result[i] = dd::CVec(A[i].size());
        for (luint j = 0; j < A[i].size(); j++) {
            if(!A[i][j].approximatelyZero()) {
                result[i][j] = CC(A[i][j].r, A[i][j].i);
            }
        }
    }
    return result;
}
dd::CMat matmul(dd::CMat& A, dd::CMat& B) {
    luint rowsA = A.size(), rowsB = B.size();

    if (rowsA == 0 || rowsB == 0) { throw logic_error("Matrix without rows?"); }
    luint colsA = A[0].size(), colsB = B[0].size();

    if (colsA != rowsB) { throw logic_error("Incorrect size of matrices."); }

    // WE ASSUME ALL ELEMENTS IN A AND B HAVE THE SAME SIZE
    dd::CMat result = dd::CMat(rowsA);
    for (luint i = 0; i < rowsA; i++) {
        result[i] = dd::CVec(colsB);
        for (luint j = 0; j < colsB; j++) {
            result[i][j] = CC(0);
            for (luint k = 0; k < colsA; k ++) {
                result[i][j] += A[i][k]*B[k][j];
            }
        }
    }

    return result;
}
dd::CMat matmul(vector<CCSparseVector>& A, dd::CMat& B) {
    luint rowsA = A.size(), rowsB = B.size();

    if (rowsA == 0 || rowsB == 0) { throw logic_error("Matrix without rows?"); }
    luint colsA = A[0].dimension(), colsB = B[0].size();

    if (colsA != rowsB) { throw logic_error("Incorrect size of matrices."); }

    // WE ASSUME ALL ELEMENTS IN A AND B HAVE THE SAME SIZE
    dd::CMat result = dd::CMat(rowsA);
    for (luint i = 0; i < rowsA; i++) {
        result[i] = dd::CVec(colsB);
        for (luint j = 0; j < colsB; j++) {
            result[i][j] = CC(0);

            unordered_set<luint>::iterator it = A[i].nonzero_iterator();
            while(it != A[i].nonzero_iterator_end()) {
                result[i][j] += A[i][*it]*B[*it][j];
                it++;
            }
        }
    }

    return result;
}
dd::CMat matmul(dd::CMat& A, vector<CCSparseVector>& B) {
    luint rowsA = A.size(), rowsB = B.size();

    if (rowsA == 0 || rowsB == 0) { throw logic_error("Matrix without rows?"); }
    luint colsA = A[0].size(), colsB = B[0].dimension();

    if (colsA != rowsB) { throw logic_error("Incorrect size of matrices."); }

    // WE ASSUME ALL ELEMENTS IN A AND B HAVE THE SAME SIZE
    dd::CMat result = dd::CMat(rowsA);
    for (luint i = 0; i < rowsA; i++) {
        result[i] = dd::CVec(colsB);
        for (luint j = 0; j < colsB; j++) {
            result[i][j] = CC(0);
            for (luint k = 0; k < colsA; k ++) {
                result[i][j] += A[i][k]*B[k][j];
            }
        }
    }

    return result;
}
dd::CMat matmul(vector<CCSparseVector>& A, vector<CCSparseVector>& B) {
    luint rowsA = A.size(), rowsB = B.size();

    if (rowsA == 0 || rowsB == 0) { throw logic_error("Matrix without rows?"); }
    luint colsA = A[0].dimension(), colsB = B[0].dimension();

    if (colsA != rowsB) { throw logic_error("Incorrect size of matrices."); }

    // WE ASSUME ALL ELEMENTS IN A AND B HAVE THE SAME SIZE
    dd::CMat result = dd::CMat(rowsA);
    for (luint i = 0; i < rowsA; i++) {
        result[i] = dd::CVec(colsB);
        for (luint j = 0; j < colsB; j++) {
            result[i][j] = CC(0);

            unordered_set<luint>::iterator it = A[i].nonzero_iterator();
            while(it != A[i].nonzero_iterator_end()) {
                result[i][j] += A[i][*it]*B[*it][j];
                it++;
            }
        }
    }

    return result;
}

dd::CMat identity_matrix(luint size) {
    dd::CMat result = dd::CMat(size);
    for (luint i = 0; i < size; i++) {
        result[i] = dd::CVec(size);
        result[i][i] = CC(1);
    }
    return result;
}

dd::CMat matrix_power(dd::CMat& M, luint t) {
    if (! is_square(M)) {
        throw logic_error("The matrix must be square");
    }
    dd::CMat R, B, I = identity_matrix(M.size());
    luint i;
    if(t & 1) {        //Handle odd values of t (this saves a multiplication later)
        R = M;
        t = t & ~1UL;    //Clear the least significant bit of t
    } else {
        R = I;
    }
    i=1;
    B=M;                //B will always be M^i, where i is a power of 2
    while (t!=0) {
        i = i*2;         //Advance i to the next power of 2
        B = matmul(B,B);         //B was M^(i/2) and is now M^i

        if(t & i) {       //i is of the form 2^j. Is the j-th bit of t set?
            R = matmul(R,B);      //Multiply the result with B=A^i
            t = t & ~i;   //Clear the j-th bit of t
        }
    }

    return R;
}
dd::CMat matrix_power(vector<CCSparseVector>& M, luint t) {
    dd::CMat denseM = sparse_to_dense(M);
    return matrix_power(denseM, t);
}

dd::vEdge* conjugate_edge(dd::vEdge*, std::unique_ptr<dd::Package<>>&);
dd::vNode* conjugate_node(dd::vNode*, std::unique_ptr<dd::Package<>>&);

dd::vEdge* conjugate_edge(dd::vEdge* v, std::unique_ptr<dd::Package<>>& package) {
    dd::vEdge* result = new dd::vEdge{conjugate_node(v->p, package), package->cn.conj(v->w)};
    return result;
}
dd::vNode* conjugate_node(dd::vNode* p, std::unique_ptr<dd::Package<>>& package) {
    if (dd::vNode::isTerminal(p)) {
        return p;
    }

    // The node is not terminal
    dd::vNode* result = new dd::vNode();
    dd::vEdge* zero_out = nullptr;
    dd::vEdge* one_out = nullptr;
    
    // Creating the edges of the new node
    if (p->e[0].w.approximatelyZero()) { // We create a zero edge
        zero_out = new dd::vEdge{dd::vNode::getTerminal(), package->cn.lookup(0)};
    }
    if (p->e[1].w.approximatelyZero()) { // We create a zero edge
        one_out = new dd::vEdge{dd::vNode::getTerminal(), package->cn.lookup(0)};
    }

    if (zero_out == nullptr) { // zero edge was not zero
        zero_out = conjugate_edge(&p->e[0], package);
        if (&p->e[0] == &p->e[1]) {
            one_out = zero_out;
        }
    } 
    if (one_out == nullptr) {
        one_out = conjugate_edge(&p->e[1], package);
    }

    dd::vNode* next = nullptr;
    // Conjugating also the reference if exits
    if (p->next != nullptr) {
        next = conjugate_node(p->next, package);
    }

    // We assign all the elements
    result->e[0] = *zero_out;
    result->e[1] = *one_out;
    result->next = next;
    result->v = p->v;

    return result;    
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
    std::unique_ptr<dd::Package<>> package = std::make_unique<dd::Package<>>(this->qbits);

    for (std::pair<dd::vEdge,CC> this_part : this->components) {
        for (std::pair<dd::vEdge,CC> other_part : vector.components) {
            dd::ComplexValue value = package->innerProduct(this_part.first, other_part.first);
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
DDVector* DDVector::conjugate() {
    DDVector* copy = new DDVector(*this);
    copy->conjugate_in();

    return copy;
}
void DDVector::conjugate_in() {
    std::unordered_map<dd::vEdge, CC> new_components;
    std::unique_ptr<dd::Package<>> package = std::make_unique<dd::Package<>>(this->nQbits());
    for (std::pair<dd::vEdge, CC> ppair : this->components) {
        dd::vEdge* conj_edge = conjugate_edge(&ppair.first, package);
        new_components[*conj_edge] = CC(ppair.second.real(), -ppair.second.imag());
    }
    this->components.clear();
    for (std::pair<dd::vEdge, CC> ppair : new_components) {
        this->components[ppair.first] = ppair.second;
    }
}
DDVector* DDVector::normalize() {
    DDVector* copy = new DDVector(*this);
    copy->normalize_in();

    return copy;
}
void DDVector::normalize_in() {
    CC inv_norm = CC(1/this->norm());
    this->operator*=(inv_norm);
}

DDVector DDVector::apply_circuit(const dd::mEdge& circuit) {
    std::unique_ptr<dd::Package<>> package = std::make_unique<dd::Package<>>(this->qbits);

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
    // std::unique_ptr<dd::Package<>> new_package = std::make_unique<dd::Package<>>(this->nQbits);
    dd::vEdge* result = new dd::vEdge(dd::simulate<>(&M, *v, this->package));
    return result;
}
dd::vEdge* FullDDSubspace::scale(dd::vEdge* v, dd::ComplexValue c) {
    dd::Complex new_value = this->package->cn.lookup(c * dd::ComplexValue(v->w));
    dd::vEdge* output = new dd::vEdge{v->p, new_value}; // Creating the edge with the new weight
    return output;
}
dd::vEdge* FullDDSubspace::add(dd::vEdge* u, dd::vEdge* v) {
    dd::vEdge* result = new dd::vEdge(this->package->add(*u, *v));
    return result;
}
dd::ComplexValue FullDDSubspace::inner_product(dd::vEdge* u, dd::vEdge* v) {
    return this->package->innerProduct(*u, *v);
}
dd::vEdge* FullDDSubspace::conjugate(dd::vEdge* v) {
    dd::vEdge* conj = conjugate_edge(v, this->package);
    return conj;
}