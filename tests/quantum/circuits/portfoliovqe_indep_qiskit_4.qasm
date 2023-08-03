// Benchmark was created by MQT Bench on 2022-12-15
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 0.2.2
// Qiskit version: {'qiskit-terra': '0.22.3', 'qiskit-aer': '0.11.1', 'qiskit-ignis': '0.7.0', 'qiskit-ibmq-provider': '0.19.2', 'qiskit': '0.39.3', 'qiskit-nature': '0.5.1', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.4.0', 'qiskit-machine-learning': '0.5.0'}

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
ry(4.38902787056145) q[0];
ry(-0.773458764355351) q[1];
cz q[0],q[1];
ry(0.176876102258767) q[2];
cz q[0],q[2];
cz q[1],q[2];
ry(-0.481235900285375) q[3];
cz q[0],q[3];
ry(4.11674195428327) q[0];
cz q[1],q[3];
ry(3.55192487880099) q[1];
cz q[0],q[1];
cz q[2],q[3];
ry(0.891578107475493) q[2];
cz q[0],q[2];
cz q[1],q[2];
ry(2.18965537790507) q[3];
cz q[0],q[3];
ry(5.15512553806975) q[0];
cz q[1],q[3];
ry(-2.03561082582181) q[1];
cz q[0],q[1];
cz q[2],q[3];
ry(3.92551945027863) q[2];
cz q[0],q[2];
cz q[1],q[2];
ry(3.56703428476674) q[3];
cz q[0],q[3];
ry(0.877898827626657) q[0];
cz q[1],q[3];
ry(-3.67671218635784) q[1];
cz q[2],q[3];
ry(-3.31120487021356) q[2];
ry(3.07010267540377) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
