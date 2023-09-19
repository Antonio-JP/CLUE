OPENQASM 2.0;
include "qelib1.inc";
gate mcphase(param0) q0,q1,q2 { cp(-0.0005) q1,q2; cx q1,q0; cp(0.0005) q0,q2; cx q1,q0; cp(-0.0005) q0,q2; }
gate gate_FTT(param0) q0,q1,q2 { cx q0,q1; ccx q0,q1,q2; mcphase(-0.001) q0,q1,q2; ccx q0,q1,q2; cx q0,q1; }
gate gate_FFT(param0) q0,q1,q2 { ccx q0,q1,q2; mcphase(-0.001) q0,q1,q2; ccx q0,q1,q2; }
gate gate_TTT(param0) q0,q1,q2 { x q0; cx q0,q1; ccx q0,q1,q2; mcphase(-0.001) q0,q1,q2; ccx q0,q1,q2; cx q0,q1; x q0; }
gate gate_FT(param0) q0,q1 { cx q0,q1; cp(-0.001) q0,q1; cx q0,q1; }
gate gate_FF(param0) q0,q1 { cp(-0.001) q0,q1; }
qreg q[9];
gate_FTT(0.001) q[8],q[3],q[1];
gate_FFT(0.001) q[7],q[6],q[3];
gate_TTT(0.001) q[3],q[7],q[5];
gate_FT(0.001) q[8],q[7];
gate_FTT(0.001) q[4],q[3],q[2];
gate_TTT(0.001) q[6],q[1],q[2];
gate_TTT(0.001) q[6],q[4],q[7];
gate_FFT(0.001) q[3],q[6],q[1];
gate_FFT(0.001) q[1],q[7],q[6];
gate_FFT(0.001) q[5],q[8],q[2];
gate_FT(0.001) q[1],q[6];
gate_FF(0.001) q[1],q[8];


// True formula: [(x_3 v x_1 v -x_8) ^ (x_3 v -x_7 v -x_6) ^ (x_3 v x_7 v x_5) ^ (x_7 v -x_8) ^ (x_3 v -x_4 v x_2) ^ (x_6 v x_1 v x_2) ^ (x_6 v x_4 v x_7) ^ (x_1 v -x_3 v -x_6) ^ (-x_1 v x_6 v -x_7) ^ (-x_5 v x_2 v -x_8) ^ (-x_1 v x_6) ^ (-x_1 v -x_8)]
// Variables appearing: {1, 2, 3, 4, 5, 6, 7, 8}
// Appear all: False

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