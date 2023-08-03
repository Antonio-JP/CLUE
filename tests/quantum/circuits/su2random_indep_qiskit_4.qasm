// Benchmark was created by MQT Bench on 2022-12-15
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 0.2.2
// Qiskit version: {'qiskit-terra': '0.22.3', 'qiskit-aer': '0.11.1', 'qiskit-ignis': '0.7.0', 'qiskit-ibmq-provider': '0.19.2', 'qiskit': '0.39.3', 'qiskit-nature': '0.5.1', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.4.0', 'qiskit-machine-learning': '0.5.0'}

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
u3(0.49359606,0.65578004,0) q[0];
u3(0.15851076,0.12050369,0) q[1];
cx q[0],q[1];
u3(0.73686279,0.88603985,0) q[2];
cx q[0],q[2];
cx q[1],q[2];
u3(0.5616796,0.0063691558,0) q[3];
cx q[0],q[3];
u3(0.30048709,0.49862556,0) q[0];
cx q[1],q[3];
u3(0.66071535,0.3046931,0) q[1];
cx q[0],q[1];
cx q[2],q[3];
u3(0.9547313,0.24638304,0) q[2];
cx q[0],q[2];
cx q[1],q[2];
u3(0.88076822,0.14141282,0) q[3];
cx q[0],q[3];
u3(0.90058873,0.14505921,0) q[0];
cx q[1],q[3];
u3(0.9751974,0.68638092,0) q[1];
cx q[0],q[1];
cx q[2],q[3];
u3(0.96590763,0.64059651,0) q[2];
cx q[0],q[2];
cx q[1],q[2];
u3(0.017968323,0.62482142,0) q[3];
cx q[0],q[3];
u3(0.24963773,0.99315677,0) q[0];
cx q[1],q[3];
u3(0.057983376,0.54311719,0) q[1];
cx q[2],q[3];
u3(0.58966328,0.42377176,0) q[2];
u3(0.64987671,0.78179897,0) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
