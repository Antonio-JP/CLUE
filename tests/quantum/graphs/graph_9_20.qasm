OPENQASM 2.0;
include "qelib1.inc";
gate gate_E(param0) q0,q1 { x q0; cp(0.001) q0,q1; x q0; x q1; cp(0.001) q1,q0; x q1; }
qreg q[9];
gate_E(0.001) q[3],q[7];
gate_E(0.001) q[5],q[7];
gate_E(0.001) q[8],q[3];
gate_E(0.001) q[1],q[6];
gate_E(0.001) q[2],q[5];
gate_E(0.001) q[4],q[5];
gate_E(0.001) q[3],q[6];
gate_E(0.001) q[8],q[5];
gate_E(0.001) q[0],q[1];
gate_E(0.001) q[2],q[4];
gate_E(0.001) q[1],q[2];
gate_E(0.001) q[0],q[4];
gate_E(0.001) q[2],q[7];
gate_E(0.001) q[1],q[5];
gate_E(0.001) q[3],q[5];
gate_E(0.001) q[8],q[4];
gate_E(0.001) q[8],q[1];
gate_E(0.001) q[0],q[3];
gate_E(0.001) q[2],q[3];
gate_E(0.001) q[1],q[7];


// Description of the graph:
//	 * Vertices: [0, 1, 2, 3, 4, 5, 6, 7, 8]
//	 * Edges: [(3, 7), (5, 7), (8, 3), (1, 6), (2, 5), (4, 5), (3, 6), (8, 5), (0, 1), (2, 4), (1, 2), (0, 4), (2, 7), (1, 5), (3, 5), (8, 4), (8, 1), (0, 3), (2, 3), (1, 7)]
//	 * Number of edges: 20


// For creating the entangled state:
//h q[0];
//h q[1];
//h q[2];
//h q[3];
//h q[4];
//h q[5];
//h q[6];
//h q[7];
//h q[8];
//barrier q;