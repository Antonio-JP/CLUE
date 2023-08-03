// Benchmark was created by MQT Bench on 2022-12-15
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 0.2.2
// Qiskit version: {'qiskit-terra': '0.22.3', 'qiskit-aer': '0.11.1', 'qiskit-ignis': '0.7.0', 'qiskit-ibmq-provider': '0.19.2', 'qiskit': '0.39.3', 'qiskit-nature': '0.5.1', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.4.0', 'qiskit-machine-learning': '0.5.0'}

OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg meas[5];
ry(1.7376841) q[0];
ry(1.459802) q[1];
cx q[1],q[0];
ry(1.3166149) q[0];
cx q[1],q[0];
ry(3*pi/8) q[2];
cry(0) q[0],q[2];
cry(0) q[1],q[2];
x q[1];
x q[3];
x q[4];
ccx q[1],q[4],q[3];
x q[1];
cx q[3],q[2];
u(0.29425236,0,0) q[2];
cx q[3],q[2];
u3(0.29425236,-pi,-pi) q[2];
cx q[3],q[2];
u(-0.11449191,0,0) q[2];
cx q[3],q[2];
u(0.11449191,0,0) q[2];
ccx q[3],q[0],q[2];
cx q[3],q[2];
u(0.11449191,0,0) q[2];
cx q[3],q[2];
u(-0.11449191,0,0) q[2];
ccx q[3],q[0],q[2];
cx q[3],q[2];
u(-0.22898382,0,0) q[2];
cx q[3],q[2];
u(0.22898382,0,0) q[2];
ccx q[3],q[1],q[2];
cx q[3],q[2];
u(0.22898382,0,0) q[2];
cx q[3],q[2];
u(-0.22898382,0,0) q[2];
ccx q[3],q[1],q[2];
x q[1];
ccx q[1],q[4],q[3];
x q[1];
x q[3];
x q[4];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
measure q[4] -> meas[4];
