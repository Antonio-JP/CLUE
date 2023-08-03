// Benchmark was created by MQT Bench on 2022-12-15
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 0.2.2
// Qiskit version: {'qiskit-terra': '0.22.3', 'qiskit-aer': '0.11.1', 'qiskit-ignis': '0.7.0', 'qiskit-ibmq-provider': '0.19.2', 'qiskit': '0.39.3', 'qiskit-nature': '0.5.1', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.4.0', 'qiskit-machine-learning': '0.5.0'}

OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg meas[3];
ry(0.836963112359282) q[0];
ry(0.571083728275385) q[1];
cx q[0],q[1];
ry(0.0636784390368105) q[2];
cx q[0],q[2];
ry(0.281039345928896) q[0];
cx q[1],q[2];
ry(0.0809992491384852) q[1];
cx q[0],q[1];
ry(0.546810273029464) q[2];
cx q[0],q[2];
ry(0.0555633185371374) q[0];
cx q[1],q[2];
ry(0.80715686876096) q[1];
cx q[0],q[1];
ry(0.288211461126035) q[2];
cx q[0],q[2];
ry(0.351866819664277) q[0];
cx q[1],q[2];
ry(0.131713846323771) q[1];
ry(0.152013197833987) q[2];
barrier q[0],q[1],q[2];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
