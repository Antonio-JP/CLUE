// Benchmark was created by MQT Bench on 2022-12-15
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 0.2.2
// Qiskit version: {'qiskit-terra': '0.22.3', 'qiskit-aer': '0.11.1', 'qiskit-ignis': '0.7.0', 'qiskit-ibmq-provider': '0.19.2', 'qiskit': '0.39.3', 'qiskit-nature': '0.5.1', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.4.0', 'qiskit-machine-learning': '0.5.0'}

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
ry(0.881714101527054) q[0];
ry(0.541321326560634) q[1];
cx q[0],q[1];
ry(0.424215438773392) q[2];
cx q[0],q[2];
cx q[1],q[2];
ry(0.254957685363639) q[3];
cx q[0],q[3];
ry(0.0415131418570927) q[0];
cx q[1],q[3];
ry(0.261329019810102) q[1];
cx q[0],q[1];
cx q[2],q[3];
ry(0.896908502089368) q[2];
cx q[0],q[2];
cx q[1],q[2];
ry(0.019165648191707) q[3];
cx q[0],q[3];
ry(0.440696800613162) q[0];
cx q[1],q[3];
ry(0.561042789494927) q[1];
cx q[0],q[1];
cx q[2],q[3];
ry(0.220530909281286) q[2];
cx q[0],q[2];
cx q[1],q[2];
ry(0.980837531996973) q[3];
cx q[0],q[3];
ry(0.895874387849242) q[0];
cx q[1],q[3];
ry(0.24504613179086) q[1];
cx q[2],q[3];
ry(0.878703422689043) q[2];
ry(0.79625283711645) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
