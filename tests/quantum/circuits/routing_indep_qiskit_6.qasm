// Benchmark was created by MQT Bench on 2022-12-15
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 0.2.2
// Qiskit version: {'qiskit-terra': '0.22.3', 'qiskit-aer': '0.11.1', 'qiskit-ignis': '0.7.0', 'qiskit-ibmq-provider': '0.19.2', 'qiskit': '0.39.3', 'qiskit-nature': '0.5.1', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.4.0', 'qiskit-machine-learning': '0.5.0'}

OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg meas[6];
ry(0.0119544964187103) q[0];
ry(-2.56513408222355) q[1];
ry(0.485645025987857) q[2];
ry(-2.82453001030404) q[3];
ry(1.7953016103785) q[4];
ry(0.141182531355062) q[5];
cx q[4],q[5];
cx q[3],q[4];
cx q[2],q[3];
cx q[1],q[2];
cx q[0],q[1];
ry(1.39525820318052) q[0];
ry(-0.997689575281629) q[1];
ry(0.366263938801539) q[2];
ry(-1.10514090564729) q[3];
ry(1.76726071192933) q[4];
ry(1.39316663875762) q[5];
cx q[4],q[5];
cx q[3],q[4];
cx q[2],q[3];
cx q[1],q[2];
cx q[0],q[1];
ry(-1.48758310874485) q[0];
ry(1.53252523867519) q[1];
ry(1.97802271535056) q[2];
ry(2.42087122214994) q[3];
ry(-0.40296440665647) q[4];
ry(1.69396275174505) q[5];
cx q[4],q[5];
cx q[3],q[4];
cx q[2],q[3];
cx q[1],q[2];
cx q[0],q[1];
ry(-2.0028886728069) q[0];
ry(0.166220166411288) q[1];
ry(-2.63092405326779) q[2];
ry(2.22367481594017) q[3];
ry(-1.59845441468004) q[4];
ry(-1.46356841205796) q[5];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
measure q[4] -> meas[4];
measure q[5] -> meas[5];
