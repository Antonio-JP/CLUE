// Benchmark was created by MQT Bench on 2022-12-15
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 0.2.2
// Qiskit version: {'qiskit-terra': '0.22.3', 'qiskit-aer': '0.11.1', 'qiskit-ignis': '0.7.0', 'qiskit-ibmq-provider': '0.19.2', 'qiskit': '0.39.3', 'qiskit-nature': '0.5.1', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.4.0', 'qiskit-machine-learning': '0.5.0'}

OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg meas[3];
u3(0.61463479,0.37505493,0) q[0];
u3(0.54365968,0.54745756,0) q[1];
cx q[0],q[1];
u3(0.63359363,0.65545404,0) q[2];
cx q[0],q[2];
u3(0.40305041,0.47316479,0) q[0];
cx q[1],q[2];
u3(0.40676977,0.47117934,0) q[1];
cx q[0],q[1];
u3(0.15296443,0.38306456,0) q[2];
cx q[0],q[2];
u3(0.7496615,0.82291163,0) q[0];
cx q[1],q[2];
u3(0.39807897,0.8667263,0) q[1];
cx q[0],q[1];
u3(0.91021005,0.69816765,0) q[2];
cx q[0],q[2];
u3(0.19259252,0.5496988,0) q[0];
cx q[1],q[2];
u3(0.21696966,0.94059241,0) q[1];
u3(0.35115908,0.82452038,0) q[2];
barrier q[0],q[1],q[2];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
