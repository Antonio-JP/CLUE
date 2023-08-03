// Benchmark was created by MQT Bench on 2022-12-15
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 0.2.2
// Qiskit version: {'qiskit-terra': '0.22.3', 'qiskit-aer': '0.11.1', 'qiskit-ignis': '0.7.0', 'qiskit-ibmq-provider': '0.19.2', 'qiskit': '0.39.3', 'qiskit-nature': '0.5.1', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.4.0', 'qiskit-machine-learning': '0.5.0'}

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
creg meas1[4];
u2(0.024029451,-pi) q[0];
u2(-0.014097935,-pi) q[1];
rzz(-5.63854173437551) q[0],q[1];
u2(-0.021434588,-pi) q[2];
rzz(-5.63761712760369) q[0],q[2];
rzz(-5.63855070693831) q[1],q[2];
u2(-0.003560167,-pi) q[3];
rzz(-5.63856045196477) q[0],q[3];
u3(0.76707631,-1.5729124,pi/2) q[0];
rzz(-5.63853639122568) q[1],q[3];
u3(0.76707631,-1.5695548,pi/2) q[1];
rzz(0.49654376455305) q[0],q[1];
rzz(-5.63872626127154) q[2],q[3];
u3(0.76707631,-1.5689087,pi/2) q[2];
rzz(0.496462341421177) q[0],q[2];
rzz(0.496544554698857) q[1],q[2];
u3(0.76707631,-1.5704828,pi/2) q[3];
rzz(0.496545412869702) q[0],q[3];
u3(1.8240681,1.5540439,-pi/2) q[0];
rzz(0.496543294022219) q[1],q[3];
u3(1.8240681,1.5806249,-pi/2) q[1];
rzz(3.93097781444086) q[0],q[1];
rzz(0.49656001444246) q[2],q[3];
u3(1.8240681,1.5857397,-pi/2) q[2];
rzz(3.9303332136064) q[0],q[2];
rzz(3.93098406977193) q[1],q[2];
u3(1.8240681,1.5732783,-pi/2) q[3];
rzz(3.93099086363542) q[0],q[3];
rx(-4.90826769756198) q[0];
rzz(3.9309740893991) q[1],q[3];
rx(-4.90826769756198) q[1];
rzz(3.93110645960634) q[2],q[3];
rx(-4.90826769756198) q[2];
rx(-4.90826769756198) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas1[0];
measure q[1] -> meas1[1];
measure q[2] -> meas1[2];
measure q[3] -> meas1[3];
