// Benchmark was created by MQT Bench on 2022-12-15
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 0.2.2
// Qiskit version: {'qiskit-terra': '0.22.3', 'qiskit-aer': '0.11.1', 'qiskit-ignis': '0.7.0', 'qiskit-ibmq-provider': '0.19.2', 'qiskit': '0.39.3', 'qiskit-nature': '0.5.1', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.4.0', 'qiskit-machine-learning': '0.5.0'}

OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg meas[6];
ry(0.151280675134311) q[0];
ry(0.506347858781068) q[1];
cx q[0],q[1];
ry(0.256177473304505) q[2];
cx q[0],q[2];
cx q[1],q[2];
ry(0.3074331924179) q[3];
cx q[0],q[3];
cx q[1],q[3];
cx q[2],q[3];
ry(0.373490442207367) q[4];
cx q[0],q[4];
cx q[1],q[4];
cx q[2],q[4];
cx q[3],q[4];
ry(0.818512009169715) q[5];
cx q[0],q[5];
ry(0.222884802123806) q[0];
cx q[1],q[5];
ry(0.691868680193051) q[1];
cx q[0],q[1];
cx q[2],q[5];
ry(0.205534200837773) q[2];
cx q[0],q[2];
cx q[1],q[2];
cx q[3],q[5];
ry(0.508595787840571) q[3];
cx q[0],q[3];
cx q[1],q[3];
cx q[2],q[3];
cx q[4],q[5];
ry(0.542983787546901) q[4];
cx q[0],q[4];
cx q[1],q[4];
cx q[2],q[4];
cx q[3],q[4];
ry(0.313780002714196) q[5];
cx q[0],q[5];
ry(0.499896342536338) q[0];
cx q[1],q[5];
ry(0.129423162574339) q[1];
cx q[0],q[1];
cx q[2],q[5];
ry(0.733681007626881) q[2];
cx q[0],q[2];
cx q[1],q[2];
cx q[3],q[5];
ry(0.263108580489538) q[3];
cx q[0],q[3];
cx q[1],q[3];
cx q[2],q[3];
cx q[4],q[5];
ry(0.602220297367877) q[4];
cx q[0],q[4];
cx q[1],q[4];
cx q[2],q[4];
cx q[3],q[4];
ry(0.192286573569491) q[5];
cx q[0],q[5];
ry(0.609687671121186) q[0];
cx q[1],q[5];
ry(0.909837763931131) q[1];
cx q[2],q[5];
ry(0.940667478186107) q[2];
cx q[3],q[5];
ry(0.821929951524639) q[3];
cx q[4],q[5];
ry(0.0644004488878683) q[4];
ry(0.544385772261176) q[5];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
// measure q[0] -> meas[0];
// measure q[1] -> meas[1];
// measure q[2] -> meas[2];
// measure q[3] -> meas[3];
// measure q[4] -> meas[4];
// measure q[5] -> meas[5];
