OPENQASM 2.0;
include "qelib1.inc";
gate mcphase(param0) q0,q1,q2 { cp(-0.0005) q1,q2; cx q1,q0; cp(0.0005) q0,q2; cx q1,q0; cp(-0.0005) q0,q2; }
gate gate_FTT(param0) q0,q1,q2 { cx q0,q1; ccx q0,q1,q2; mcphase(-0.001) q0,q1,q2; ccx q0,q1,q2; cx q0,q1; }
gate gate_FFT(param0) q0,q1,q2 { ccx q0,q1,q2; mcphase(-0.001) q0,q1,q2; ccx q0,q1,q2; }
qreg q[5];
gate_FTT(0.001) q[1],q[3],q[0];
gate_FFT(0.001) q[4],q[3],q[2];


// True formula: [(-x_1 v x_3 v x_0) ^ (-x_4 v x_2 v -x_3)]
// Variables appearing: {0, 1, 2, 3, 4}
// Appear all: True

// For creating the entangled state:
//h q[0];
//h q[1];
//h q[2];
//h q[3];
//h q[4];
//barrier q;