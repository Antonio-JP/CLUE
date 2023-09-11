OPENQASM 2.0;
include "qelib1.inc";
gate mcphase(param0) q0,q1,q2 { cp(-pi/2) q1,q2; cx q1,q0; cp(pi/2) q0,q2; cx q1,q0; cp(-pi/2) q0,q2; }
gate gate_TTT q0,q1,q2 { x q0; cx q0,q1; ccx q0,q1,q2; mcphase(-pi) q0,q1,q2; ccx q0,q1,q2; cx q0,q1; x q0; }
gate gate_FFT q0,q1,q2 { ccx q0,q1,q2; mcphase(-pi) q0,q1,q2; ccx q0,q1,q2; }
gate gate_FTT q0,q1,q2 { cx q0,q1; ccx q0,q1,q2; mcphase(-pi) q0,q1,q2; ccx q0,q1,q2; cx q0,q1; }
gate gate_FT q0,q1 { cx q0,q1; cp(-pi) q0,q1; cx q0,q1; }
qreg q[10];
gate_TTT q[3],q[5],q[0];
gate_TTT q[0],q[4],q[8];
gate_FFT q[6],q[3],q[0];
gate_FTT q[9],q[3],q[5];
gate_FT q[1],q[0];
gate_TTT q[1],q[5],q[8];
gate_TTT q[7],q[2],q[8];
gate_FTT q[8],q[0],q[1];
gate_FTT q[0],q[1],q[2];
gate_FFT q[8],q[0],q[3];
gate_FTT q[4],q[3],q[6];
gate_FFT q[8],q[6],q[2];
gate_FTT q[9],q[0],q[2];


// True formula: [(x_3 v x_5 v x_0) ^ (x_0 v x_4 v x_8) ^ (x_0 v -x_6 v -x_3) ^ (x_3 v -x_9 v x_5) ^ (x_0 v -x_1) ^ (x_1 v x_5 v x_8) ^ (x_7 v x_2 v x_8) ^ (x_0 v x_1 v -x_8) ^ (x_1 v x_2 v -x_0) ^ (x_3 v -x_8 v -x_0) ^ (x_3 v x_6 v -x_4) ^ (x_2 v -x_8 v -x_6) ^ (x_0 v -x_9 v x_2)]
// Variables appearing: {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
// Appear all: True
