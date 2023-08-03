// Benchmark was created by MQT Bench on 2022-12-15
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 0.2.2
// Qiskit version: {'qiskit-terra': '0.22.3', 'qiskit-aer': '0.11.1', 'qiskit-ignis': '0.7.0', 'qiskit-ibmq-provider': '0.19.2', 'qiskit': '0.39.3', 'qiskit-nature': '0.5.1', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.4.0', 'qiskit-machine-learning': '0.5.0'}

OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
creg meas[2];
ry(2.28643189871332) q[0];
ry(-1.986944517459) q[1];
cx q[0],q[1];
ry(1.09267111867709) q[0];
ry(-1.28315816917591) q[1];
cx q[0],q[1];
ry(0.644936496959995) q[0];
ry(-2.32308793116108) q[1];
cx q[0],q[1];
ry(1.46213217253153) q[0];
ry(2.4135628694577) q[1];
barrier q[0],q[1];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
