OPENQASM 2.0;
include "qelib1.inc";
gate mcphase(param0) q0,q1,q2 { cp(-pi/2) q1,q2; cx q1,q0; cp(pi/2) q0,q2; cx q1,q0; cp(-pi/2) q0,q2; }
gate gate_TTF_3_141592653589793_ q0,q1,q2 { ccx q0,q1,q2; mcphase(-pi) q0,q1,q2; ccx q0,q1,q2; }
gate gate_FT_3_141592653589793_ q0,q1 { cx q0,q1; cp(-pi) q0,q1; cx q0,q1; }
gate gate_FTT_3_141592653589793_ q0,q1,q2 { cx q0,q1; ccx q0,q1,q2; mcphase(-pi) q0,q1,q2; ccx q0,q1,q2; cx q0,q1; }
gate gate_TFF_3_141592653589793_ q0,q1,q2 { x q0; cx q0,q1; ccx q0,q1,q2; mcphase(-pi) q0,q1,q2; ccx q0,q1,q2; cx q0,q1; x q0; }
gate gate_FFF_3_141592653589793_ q0,q1,q2 { mcphase(-pi) q0,q1,q2; }
gate gate_TT_3_141592653589793_ q0,q1 { x q0; cx q0,q1; cp(-pi) q0,q1; cx q0,q1; x q0; }
qreg q[10];
gate_TTF_3_141592653589793_ q[4],q[6],q[5];
gate_FT_3_141592653589793_ q[1],q[5];
gate_FTT_3_141592653589793_ q[3],q[6],q[4];
gate_TFF_3_141592653589793_ q[3],q[1],q[4];
gate_FT_3_141592653589793_ q[2],q[9];
gate_FFF_3_141592653589793_ q[7],q[6],q[3];
gate_TTF_3_141592653589793_ q[1],q[3],q[2];
gate_FTT_3_141592653589793_ q[5],q[4],q[8];
gate_TT_3_141592653589793_ q[0],q[5];
gate_FTT_3_141592653589793_ q[3],q[0],q[9];
gate_TTF_3_141592653589793_ q[9],q[6],q[5];
