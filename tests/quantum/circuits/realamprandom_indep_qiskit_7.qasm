// Benchmark was created by MQT Bench on 2022-12-15
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 0.2.2
// Qiskit version: {'qiskit-terra': '0.22.3', 'qiskit-aer': '0.11.1', 'qiskit-ignis': '0.7.0', 'qiskit-ibmq-provider': '0.19.2', 'qiskit': '0.39.3', 'qiskit-nature': '0.5.1', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.4.0', 'qiskit-machine-learning': '0.5.0'}

OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg meas[7];
ry(0.783682709352997) q[0];
ry(0.659768782870392) q[1];
cx q[0],q[1];
ry(0.655514283359269) q[2];
cx q[0],q[2];
cx q[1],q[2];
ry(0.751855011871936) q[3];
cx q[0],q[3];
cx q[1],q[3];
cx q[2],q[3];
ry(0.11276168580456) q[4];
cx q[0],q[4];
cx q[1],q[4];
cx q[2],q[4];
cx q[3],q[4];
ry(0.0330204136625072) q[5];
cx q[0],q[5];
cx q[1],q[5];
cx q[2],q[5];
cx q[3],q[5];
cx q[4],q[5];
ry(0.697857471840325) q[6];
cx q[0],q[6];
ry(0.942143441444552) q[0];
cx q[1],q[6];
ry(0.436907412462656) q[1];
cx q[0],q[1];
cx q[2],q[6];
ry(0.989183385317161) q[2];
cx q[0],q[2];
cx q[1],q[2];
cx q[3],q[6];
ry(0.63396522145249) q[3];
cx q[0],q[3];
cx q[1],q[3];
cx q[2],q[3];
cx q[4],q[6];
ry(0.13439356153129) q[4];
cx q[0],q[4];
cx q[1],q[4];
cx q[2],q[4];
cx q[3],q[4];
cx q[5],q[6];
ry(0.818720622344007) q[5];
cx q[0],q[5];
cx q[1],q[5];
cx q[2],q[5];
cx q[3],q[5];
cx q[4],q[5];
ry(0.998452330083134) q[6];
cx q[0],q[6];
ry(0.383130249629609) q[0];
cx q[1],q[6];
ry(0.594238834266465) q[1];
cx q[0],q[1];
cx q[2],q[6];
ry(0.407173230160596) q[2];
cx q[0],q[2];
cx q[1],q[2];
cx q[3],q[6];
ry(0.991509901292922) q[3];
cx q[0],q[3];
cx q[1],q[3];
cx q[2],q[3];
cx q[4],q[6];
ry(0.490493392798058) q[4];
cx q[0],q[4];
cx q[1],q[4];
cx q[2],q[4];
cx q[3],q[4];
cx q[5],q[6];
ry(0.512578486101571) q[5];
cx q[0],q[5];
cx q[1],q[5];
cx q[2],q[5];
cx q[3],q[5];
cx q[4],q[5];
ry(0.121390517445704) q[6];
cx q[0],q[6];
ry(0.301107294648713) q[0];
cx q[1],q[6];
ry(0.916199465345135) q[1];
cx q[2],q[6];
ry(0.802760350690373) q[2];
cx q[3],q[6];
ry(0.997234060884546) q[3];
cx q[4],q[6];
ry(0.776840379558342) q[4];
cx q[5],q[6];
ry(0.0319601670081443) q[5];
ry(0.684727357314778) q[6];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6];
// measure q[0] -> meas[0];
// measure q[1] -> meas[1];
// measure q[2] -> meas[2];
// measure q[3] -> meas[3];
// measure q[4] -> meas[4];
// measure q[5] -> meas[5];
// measure q[6] -> meas[6];
