// Benchmark was created by MQT Bench on 2022-12-15
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 0.2.2
// Qiskit version: {'qiskit-terra': '0.22.3', 'qiskit-aer': '0.11.1', 'qiskit-ignis': '0.7.0', 'qiskit-ibmq-provider': '0.19.2', 'qiskit': '0.39.3', 'qiskit-nature': '0.5.1', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.4.0', 'qiskit-machine-learning': '0.5.0'}

OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg meas[3];
ry(3.84165532650502) q[0];
ry(-3.44795970350275) q[1];
cz q[0],q[1];
ry(-0.599575587811508) q[2];
cz q[0],q[2];
ry(3.23220644941661) q[0];
cz q[1],q[2];
ry(-1.15012568688631) q[1];
cz q[0],q[1];
ry(-1.08550862252372) q[2];
cz q[0],q[2];
ry(3.51901603754477) q[0];
cz q[1],q[2];
ry(6.73026434236355) q[1];
cz q[0],q[1];
ry(-5.14875148093778) q[2];
cz q[0],q[2];
ry(-2.12660116372912) q[0];
cz q[1],q[2];
ry(-0.394449825925003) q[1];
ry(3.21037473955161) q[2];
barrier q[0],q[1],q[2];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
