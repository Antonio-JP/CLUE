#include "experiments/CUTExperiment.hpp"

/*** CODE FOR CLASS UNDIRECTED_GRAPH ***/

// Private section
unordered_map<luint,vector<luint>> UndirectedGraph::compute_possible_values() {
    if (this->__modified) { // Not yet computed
        // cerr << "[CUT Graph]\tComputing possible cut values of a graph..." << endl;
        luint nVertices = this->n_vertices();
        for (luint i = 0; i < static_cast<luint>(pow(2, nVertices)); i++) {
            luint new_value = this->cut_value(boost::dynamic_bitset<>(nVertices, i));
            if (!this->possible_values.contains(new_value)) {
                this->possible_values[new_value] = vector<luint>();
            }
            this->possible_values[new_value].push_back(i);
        }
        // cerr << "[CUT Graph]\tComputing possible cut values of a graph -> Done" << endl;
        this->__modified = false;
    }
    return this->possible_values;
}

// Builders for Undirected graphs
UndirectedGraph::UndirectedGraph(luint n_vertices, luint eIterations, ExperimentType eType, dd::Package<>* ePackage) 
    : Experiment("CUT", "H", eIterations, eType, ePackage) {
    // We add the first vertices indicated by the argument n_vertices
    for (luint i = 0; i < n_vertices; i++) {
        this->add_vertex(i);
    }
}

// Methods to modify the graph structure
bool UndirectedGraph::add_vertex(luint label) {
    if (this->vertices.contains(label)) { return false; } // vertex was already there
    // We add the vertex to the list of vertices and the map of edges
    this->vertices.insert(label);
    this->edges[label] = unordered_set<luint>();

    // If this vertex is the default vertex, we update the next valid vertex
    if (label == this->next_to_add) {
        luint candidate = label + 1;
        while (this->vertices.contains(candidate)) { candidate++; }
        this->next_to_add = candidate;
    }

    this->__modified = true;
    this->__nVertices++;
    return true;
}

bool UndirectedGraph::add_edge(luint src, luint trg) {
    if ((!this->vertices.contains(src)) || (!this->vertices.contains(trg))) {
        throw domain_error("[Graph] vertex for edge not present in graph");
    }

    if (this->edges[src].contains(trg)) { return false; } // edge was already present
    // We add the edge to both the sets of `src` and `trg`
    this->edges[src].insert(trg);
    this->edges[trg].insert(src);

    this->__modified = true;
    this->__nEdges++;
    return true;
}

// Static methods
/*static*/ UndirectedGraph* UndirectedGraph::random(luint n_vertices, double density, ExperimentType eType, dd::Package<>* ePackage) {
    luint iterations = static_cast<luint>(ceil(pow(2., static_cast<double>(n_vertices)/2.)));
    UndirectedGraph *G = new UndirectedGraph(n_vertices, iterations, eType, ePackage);

    vector<pair<luint,luint>> possible_edges;
    std::unordered_set<luint>::iterator v_end = G->vertices.end();
    for (std::unordered_set<luint>::iterator first = G->vertices.begin(); first != v_end; first++) {
        for (std::unordered_set<luint>::iterator second = std::next(first); second != v_end; second++) {
            possible_edges.push_back(pair<luint,luint>(*first, *second));
        }
    }
 
    for (pair<luint,luint> edge : possible_edges) {
        if ((static_cast<double>(rand()) / (RAND_MAX)) < density) {
            G->add_edge(edge.first, edge.second);
        }
    }

    return G;
}
/*static*/ UndirectedGraph* UndirectedGraph::random(luint n_vertices, luint n_edges, ExperimentType eType, dd::Package<>* ePackage) {
    luint iterations = static_cast<luint>(ceil(pow(2., static_cast<double>(n_vertices)/2.)));
    UndirectedGraph *G = new UndirectedGraph(n_vertices, iterations, eType, ePackage);

    vector<pair<luint,luint>> possible_edges;
    std::unordered_set<luint>::iterator v_end = G->vertices.end();
    for (std::unordered_set<luint>::iterator first = G->vertices.begin(); first != v_end; first++) {
        for (std::unordered_set<luint>::iterator second = std::next(first); second != v_end; second++) {
            possible_edges.push_back(pair<luint,luint>(*first, *second));
        }
    }

    for (luint i = 0; i < n_edges; i++) {
        luint to_add = static_cast<luint>(rand())%possible_edges.size();
        pair<luint,luint> edge = possible_edges[to_add];
        possible_edges.erase(possible_edges.begin()+static_cast<long int>(to_add));
        G->add_edge(edge.first,edge.second);
    }

    return G;
}

// Graph methods
vector<CCSparseVector> UndirectedGraph::adjacency_matrix() {
    vector<CCSparseVector> result = vector<CCSparseVector>(this->n_vertices(), this->n_vertices());
    luint i = 0;
    for (luint src : this->vertices) {
        luint j = 0;
        for (luint trg: this->vertices) {
            if (this->edges[src].contains(trg)) {
                result[i].set_value(j, CC(1));
            }
            j++;
        }
        i++;
    }

    return result;
}

// Methods to extract information
luint UndirectedGraph::cut_value(boost::dynamic_bitset<> cut) {
    std::unordered_set<luint> zero_set, one_set;
    for (luint i = 0; i < this->n_vertices(); i++) {
        if (cut[i]) { one_set.insert(i); }
        else { zero_set.insert(i); }
    }

    luint count = 0U;
    for (luint src : one_set) {
        for (luint trg: this->edges[src]) {
            if (zero_set.contains(trg)) { count++; }
        }
    }
    return count;
}

string UndirectedGraph::to_string() {
    stringstream output;
    if (this->n_vertices() == 0) { output << "Emtpy Graph"; }
    else {
        output << "\"Vertices (";
        std::unordered_set<luint>::iterator v = this->vertices.begin();
        output << *v; v++;
        while (v != this->vertices.end()) {
            output << ", " << *v;
            v++;
        }
        output << ") -- (Edges: [";
        for (luint src : this->vertices) {
            for (luint trg: this->edges[src]) {
                if (src < trg) {
                    output << "(" << src << ", " << trg << "),";
                }
            }
        }
        output << "])\"";
    }
    
    return output.str();
}

// Implementation of virtual methods from Experiment
luint UndirectedGraph::correct_size() {
    this->compute_possible_values();
    return this->possible_values.size();
}
array<dd::CMat, 2U> UndirectedGraph::direct() {
    luint d = this->correct_size(); // This computes possible_values
    luint full_size = static_cast<luint>(pow(2, this->n_vertices()));
    dd::CMat L = dd::CMat(d), U = dd::CMat(d); 
    luint i = 0;
    for (std::pair<luint,vector<luint>> pair : this->possible_values) {
        L[i] = dd::CVec(full_size);
        CC value = CC(1./sqrt(pair.second.size()));
        for (luint j : pair.second) { L[i][j] = value; }
        U[i] = dd::CVec(d);
        U[i][i] = exp(CC(0, -static_cast<double>(pair.first)));
        i++;
    }
    return {L,U}; // TODO
}
vector<CCSparseVector> UndirectedGraph::matrix() {
    luint full_size = static_cast<luint>(pow(2, this->n_vertices()));
    vector<CCSparseVector> result = vector<CCSparseVector>(full_size, full_size);
    for (luint i = 0; i < full_size; i++) {
        luint value = this->cut_value(boost::dynamic_bitset(this->n_vertices(), i));
        result[i].set_value(i, CC(static_cast<double>(value)));
    }
    return result;
}
dd::CMat UndirectedGraph::matrix_B(dd::CMat& Uhat) {
    if (is_diagonal(Uhat)) {
        // Not yet implemented: we return an identity
        dd::CMat result = dd::CMat(Uhat.size());
        for (luint i = 0; i < Uhat.size(); i++) {
            result[i] = dd::CVec(Uhat.size());
            result[i][i] = CC(1.);
        }
        return result;
        // cerr << "Found a diagonal matrix: \n" << matrix_to_string(Uhat) << endl;
        // throw logic_error("Begin matrix not implemented for diagonal matrices");
    }
    luint num_edges = this->n_edges();
    dd::CMat result = dd::CMat(Uhat.size());
    for (luint i = 0; i < Uhat.size(); i++) {
        result[i] = dd::CVec(Uhat[i].size());
        result[i][i] = CC(1./static_cast<double>(num_edges));
    }
    result[0][0] = CC(static_cast<double>(num_edges));

    return result;
}
qc::QuantumComputation* UndirectedGraph::quantum(double par_val) {
    auto* circuit = new qc::QuantumComputation(this->n_vertices());

    for (luint src : this->vertices) {
        auto s_bit = static_cast<qc::Qubit>(src);
        for (luint trg : this->edges[src]) {
            auto t_bit = static_cast<qc::Qubit>(trg);
            if (s_bit < t_bit) { // we avoid repeating as (u,v) = (v,u)
                circuit->cp(-par_val, {s_bit, qc::Control::Type::Neg}, t_bit);
                circuit->cp(-par_val, {t_bit, qc::Control::Type::Neg}, s_bit);
            }
        }
    }

    return circuit;
}
qc::QuantumComputation* UndirectedGraph::quantum_B(double par_val) {
    // Create a circuit with the appropriate number of qbits
    auto* circuit = new qc::QuantumComputation(this->n_vertices());

    // We create the new gates
    for (luint i = 0; i < this->n_vertices(); i++) {
        auto i_bit = static_cast<qc::Qubit>(i);
        circuit->h(i_bit);
        circuit->p(-par_val, i_bit);
        circuit->x(i_bit);
        circuit->p(par_val, i_bit);
        circuit->x(i_bit);
        circuit->h(i_bit);
    }
    
    return circuit;
}
UndirectedGraph* UndirectedGraph::change_exec_type(ExperimentType new_type) {
    UndirectedGraph* result = new UndirectedGraph(this->n_vertices(), this->iterations, new_type, this->package);
    for (luint src : this->vertices) {
        for (luint trg : this->edges[src]) {
            result->add_edge(src, trg);
        }
    }

    return result;
}
