# Quantum circuits

## Meeting of 2023/07/31

1. How good is the improvement of performance?
2. Decided to send to TACAS: change from engineering to more Theoretical Computer Science perspective
3. CDC: maybe send to TAC journal.
4. Prepare new table: just one with all the examples that reduce
  - Check the examples are correct in terms of the input
  - Distinguish between state preparations, algorithms and gates.

## Meeting of 2023/06/09

1. Prepare table to see data in "case_studies.ipynb" (see image)
2. Rerun the benchmark examples.

## Meeting of 2023/05/31

1. Run all the benchmarks from the website up to 7 q-bits
2. Produce case studies with several sizes to see scalability
3. Increase the content on the draft.

## Things TO DO

1. Write some theory for lumping on Quantum circuit
2. Develop the case for Grover's algorithm
   * See that the good/bad observables provide the good lumping.
   * Recover the bad/good observables from the other.
   * See why any observable behave so nicely.
   * Can we recover the good/bad observables from any observable as input to CLUE?
3. Develop the case for Shor's algorithm
   * Generate the unitary matrix without measuring
   * See that we get a simple lumping from the good observable
   * What does happen in other observables?
4. Develop the case for Quantum Fourier Transform
   * Is there a reason the $0$ and $2^{n-1}$ observables behave so nice?
5. Develop the case for Quantum Phase ...
6. Get and run the examples from meeting on 2023/05/12

## Meeting on 2023/05/12

* Run for Deutsch-Jozsa + HHL + Shor + Phase Estimation
* Try to interpret the output of some observables
* Think about the probabilities on the observables

## Other links

* [Main Benchmark website](https://www.cda.cit.tum.de/mqtbench/)
* [Qiskit Tutorial](https://github.com/Qiskit/qiskit-tutorials/blob/master/tutorials/simulators/1_aer_provider.ipynb)
* [Qiskit documentation](https://qiskit.org/documentation/)
* [How get unitary from Quantum Circuit](https://quantumcomputinguk.org/tutorials/how-to-obtain-the-unitary-matrix-of-a-circuit-in-qiskit-with-code)