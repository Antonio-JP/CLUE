// Benchmark was created by MQT Bench on 2022-12-15
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 0.2.2
// Qiskit version: {'qiskit-terra': '0.22.3', 'qiskit-aer': '0.11.1', 'qiskit-ignis': '0.7.0', 'qiskit-ibmq-provider': '0.19.2', 'qiskit': '0.39.3', 'qiskit-nature': '0.5.1', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.4.0', 'qiskit-machine-learning': '0.5.0'}

OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg meas[3];
creg meas0[3];
u2(0.41951949,-pi) q[0];
u2(0.41620669,-pi) q[1];
rzz(-0.420917333908502) q[0],q[1];
u2(0.41905329,-pi) q[2];
rzz(-0.421016123405307) q[0],q[2];
u3(2.2348228,1.2558831,-pi/2) q[0];
rzz(-0.420940441831552) q[1],q[2];
u3(2.2348228,1.2087537,-pi/2) q[1];
rzz(-5.98815838177421) q[0],q[1];
u3(2.2348228,1.2492507,-pi/2) q[2];
rzz(-5.98956380537088) q[0],q[2];
u3(0.56042125,-1.2358007,pi/2) q[0];
rzz(-5.98848712542991) q[1],q[2];
u3(0.56042125,-1.2880621,pi/2) q[1];
rzz(-6.64023274758061) q[0],q[1];
u3(0.56042125,-1.2431553,pi/2) q[2];
rzz(-6.6417912133385) q[0],q[2];
rx(-4.06512402388918) q[0];
rzz(-6.64059728943955) q[1],q[2];
rx(-4.06512402388918) q[1];
rx(-4.06512402388918) q[2];
barrier q[0],q[1],q[2];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
barrier q[0],q[1],q[2];
measure q[0] -> meas0[0];
measure q[1] -> meas0[1];
measure q[2] -> meas0[2];
