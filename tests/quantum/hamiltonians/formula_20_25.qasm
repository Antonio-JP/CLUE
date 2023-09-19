OPENQASM 2.0;
include "qelib1.inc";
gate mcphase(param0) q0,q1,q2 { cp(-0.0005) q1,q2; cx q1,q0; cp(0.0005) q0,q2; cx q1,q0; cp(-0.0005) q0,q2; }
gate gate_FTT(param0) q0,q1,q2 { cx q0,q1; ccx q0,q1,q2; mcphase(-0.001) q0,q1,q2; ccx q0,q1,q2; cx q0,q1; }
gate gate_FT(param0) q0,q1 { cx q0,q1; cp(-0.001) q0,q1; cx q0,q1; }
gate gate_FFT(param0) q0,q1,q2 { ccx q0,q1,q2; mcphase(-0.001) q0,q1,q2; ccx q0,q1,q2; }
gate gate_TTT(param0) q0,q1,q2 { x q0; cx q0,q1; ccx q0,q1,q2; mcphase(-0.001) q0,q1,q2; ccx q0,q1,q2; cx q0,q1; x q0; }
gate gate_FFF(param0) q0,q1,q2 { mcphase(-0.001) q0,q1,q2; }
qreg q[20];
gate_FTT(0.001) q[10],q[14],q[9];
gate_FT(0.001) q[6],q[14];
gate_FTT(0.001) q[7],q[1],q[12];
gate_FFT(0.001) q[0],q[6],q[3];
gate_TTT(0.001) q[12],q[9],q[15];
gate_FTT(0.001) q[8],q[9],q[4];
gate_FTT(0.001) q[7],q[9],q[13];
gate_FFT(0.001) q[10],q[6],q[4];
gate_FFT(0.001) q[12],q[3],q[8];
gate_FFT(0.001) q[8],q[14],q[4];
gate_FFT(0.001) q[18],q[5],q[12];
gate_FFT(0.001) q[15],q[11],q[5];
gate_FFT(0.001) q[15],q[14],q[13];
gate_FTT(0.001) q[3],q[14],q[8];
gate_FTT(0.001) q[19],q[11],q[6];
gate_FFF(0.001) q[4],q[18],q[16];
gate_FFT(0.001) q[19],q[11],q[17];
gate_FFT(0.001) q[9],q[0],q[8];
gate_FTT(0.001) q[0],q[11],q[19];
gate_FFF(0.001) q[9],q[13],q[14];
gate_FTT(0.001) q[15],q[14],q[8];
gate_FTT(0.001) q[15],q[1],q[7];
gate_FTT(0.001) q[14],q[12],q[7];
gate_FTT(0.001) q[0],q[17],q[4];
gate_FFT(0.001) q[10],q[5],q[0];


// True formula: [(x_14 v x_9 v -x_10) ^ (x_14 v -x_6) ^ (x_1 v x_12 v -x_7) ^ (x_3 v -x_0 v -x_6) ^ (x_12 v x_9 v x_15) ^ (x_9 v x_4 v -x_8) ^ (-x_7 v x_9 v x_13) ^ (-x_10 v x_4 v -x_6) ^ (-x_12 v -x_3 v x_8) ^ (x_4 v -x_8 v -x_14) ^ (x_12 v -x_18 v -x_5) ^ (-x_15 v x_5 v -x_11) ^ (-x_15 v x_13 v -x_14) ^ (x_14 v x_8 v -x_3) ^ (x_11 v x_6 v -x_19) ^ (-x_4 v -x_18 v -x_16) ^ (x_17 v -x_19 v -x_11) ^ (-x_9 v x_8 v -x_0) ^ (x_11 v x_19 v -x_0) ^ (-x_9 v -x_13 v -x_14) ^ (x_14 v -x_15 v x_8) ^ (x_1 v x_7 v -x_15) ^ (x_12 v x_7 v -x_14) ^ (x_17 v x_4 v -x_0) ^ (x_0 v -x_10 v -x_5)]
// Variables appearing: {0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}
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
//h q[9];
//h q[10];
//h q[11];
//h q[12];
//h q[13];
//h q[14];
//h q[15];
//h q[16];
//h q[17];
//h q[18];
//h q[19];
//barrier q;