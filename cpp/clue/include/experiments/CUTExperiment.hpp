#ifndef CLUE_EX_CUT
#define CLUE_EX_CUT

#include "Experiment.hpp"
#include "boost/dynamic_bitset.hpp"

using namespace std;

/**
 * Class for an Undirected and Unweighted Graph 
 * 
 * The vertices of the graph will be identified with an unsigned long int. 
 * Custom indices can be set on new vertices.
 * Adding edges require their vertices to be already in the graph.
*/
class UndirectedGraph : public Experiment {
    private:
        unordered_map<luint, unordered_set<luint>> edges;
        unordered_set<luint> vertices;
        luint next_to_add = 0U;
        unordered_map<luint,vector<luint>> possible_values;
        bool __modified = true;
        luint __nEdges = 0U;
        luint __nVertices = 0U;

        unordered_map<luint,vector<luint>> compute_possible_values();
    public:
        UndirectedGraph(luint, luint, ExperimentType); // Builds graphs with given number of vertices
        UndirectedGraph(luint eIterations, ExperimentType eType) : UndirectedGraph(0, eIterations, eType) { }
        ~UndirectedGraph() = default;

        bool add_vertex(luint);
        bool add_vertex() { return this->add_vertex(this->next_to_add); }
        bool add_edge(luint, luint);

        static UndirectedGraph* random(luint, double, ExperimentType = ExperimentType::DIRECT); // Random builder with probability on edges
        static UndirectedGraph* random(luint, luint, ExperimentType = ExperimentType::DIRECT); // Random builder with fixed number of edges

        luint n_vertices() { return this->__nVertices; }
        luint n_edges() { return this->__nEdges; }
        vector<CCSparseVector> adjacency_matrix();

        /* Method to count the cut value from a cut */
        luint cut_value(boost::dynamic_bitset<>);
        /* Method to transform the formula into a string */
        string to_string();

        /* Virtual methods from Experiment */
        luint size() { return this->vertices.size(); }
        luint correct_size();
        luint bound_size() { return this->n_edges(); }
        array<dd::CMat, 2U> direct();
        vector<CCSparseVector> matrix();
        dd::CMat matrix_B(dd::CMat&);
        qc::QuantumComputation* quantum(double);
        qc::QuantumComputation* quantum_B(double);
        UndirectedGraph* change_exec_type(ExperimentType);
};

#endif