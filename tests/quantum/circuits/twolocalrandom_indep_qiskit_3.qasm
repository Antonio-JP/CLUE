// Benchmark was created by MQT Bench on 2022-12-15
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 0.2.2
// Qiskit version: {'qiskit-terra': '0.22.3', 'qiskit-aer': '0.11.1', 'qiskit-ignis': '0.7.0', 'qiskit-ibmq-provider': '0.19.2', 'qiskit': '0.39.3', 'qiskit-nature': '0.5.1', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.4.0', 'qiskit-machine-learning': '0.5.0'}

OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg meas[3];
ry(0.945141397462588) q[0];
ry(0.758693271995721) q[1];
cx q[0],q[1];
ry(0.113930364498725) q[2];
cx q[0],q[2];
ry(0.659958167729383) q[0];
cx q[1],q[2];
ry(0.495388568874548) q[1];
cx q[0],q[1];
ry(0.48058770505592) q[2];
cx q[0],q[2];
ry(0.117918828788378) q[0];
cx q[1],q[2];
ry(0.639481369171448) q[1];
cx q[0],q[1];
ry(0.794817036328988) q[2];
cx q[0],q[2];
ry(0.551034098377576) q[0];
cx q[1],q[2];
ry(0.706950794203591) q[1];
ry(0.0594946190689315) q[2];
barrier q[0],q[1],q[2];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
