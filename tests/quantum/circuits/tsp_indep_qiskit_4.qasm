// Benchmark was created by MQT Bench on 2022-12-15
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 0.2.2
// Qiskit version: {'qiskit-terra': '0.22.3', 'qiskit-aer': '0.11.1', 'qiskit-ignis': '0.7.0', 'qiskit-ibmq-provider': '0.19.2', 'qiskit': '0.39.3', 'qiskit-nature': '0.5.1', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.4.0', 'qiskit-machine-learning': '0.5.0'}

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
ry(6.39895244257109) q[0];
ry(-2.65545332677397) q[1];
cz q[0],q[1];
ry(2.24399645113615) q[0];
ry(4.30074251662019) q[2];
cz q[1],q[2];
ry(-4.13127412867375) q[1];
cz q[0],q[1];
ry(-2.72510760634829) q[0];
ry(-5.00137370937142) q[3];
cz q[2],q[3];
ry(3.11789536237741) q[2];
cz q[1],q[2];
ry(7.01896016901562) q[1];
cz q[0],q[1];
ry(0.92556202682689) q[0];
ry(4.04358952341942) q[3];
cz q[2],q[3];
ry(3.4220449620582) q[2];
cz q[1],q[2];
ry(2.52311099453506) q[1];
cz q[0],q[1];
ry(-3.6733347632896) q[0];
ry(-3.0969112176004) q[3];
cz q[2],q[3];
ry(3.39875737412341) q[2];
cz q[1],q[2];
ry(5.16959104687041) q[1];
cz q[0],q[1];
ry(-3.11211705088501) q[0];
ry(6.9671370300791) q[3];
cz q[2],q[3];
ry(-3.76619289597362) q[2];
cz q[1],q[2];
ry(4.60923563707103) q[1];
ry(5.61939718326072) q[3];
cz q[2],q[3];
ry(-4.1423518633464) q[2];
ry(5.67284961632846) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
