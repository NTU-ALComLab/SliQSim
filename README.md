# SliQSim - A BDD-based Quantum Circuit Simulator

## Introduction
`SliQSim` is a BDD-based quantum circuit simulator implemented in C/C++ on top of [CUDD](http://web.mit.edu/sage/export/tmp/y/usr/share/doc/polybori/cudd/cuddIntro.html) package. In `SliQSim`, a bit-slicing technique based on BDDs is used to represent quantum state vectors. For more details of the simulator, please refer to the [paper](https://arxiv.org/abs/2007.09304).

## Build
To build the simulator, one needs to first configure `CUDD`:
```commandline
cd cudd
./configure --enable-dddmp --enable-obj --enable-shared --enable-static 
cd ..
```
Next, build the binary file `SliQSim`:
```commandline
make
```

## Execution
The circuit format being simulated is `OpenQASM` used by IBM's [Qiskit](https://github.com/Qiskit/qiskit), and the gate set supported in this simulator now contains Pauli-X (x), Pauli-Y (y), Pauli-Z (z), Hadamard (h), Phase and its inverse (s and sdg), π/8 and its inverse (t and tdg), Rotation-X with phase π/2 (rx(pi/2)), Rotation-Y with phase π/2 (ry(pi/2)), Controlled-NOT (cx), Controlled-Z (cz), Toffoli (ccx and mcx), SWAP (swap), and Fredkin (cswap). One can find some example benchmarks in `examples` folder. 

For simulation types, we provide "sampling", "all_amplitude", and "query" simulation options. The help message states the details:

```commandline
$ ./SliQSim --help
Options:
--help                produce help message
--sim_qasm arg        simulate qasm file string
--seed [=arg(=1)]     seed for random number generator
--print_info          print simulation statistics such as runtime, memory, etc.
--type arg (=0)       the simulation type being executed.
                      0: sampling mode (default option), where the sampled outcomes will be provided. 
                      1: all_amplitude mode, where the final state vector will be shown. 
                      2: query mode, where only the values of properties defined in obs_file will be provided.
--shots arg (=1)      the number of outcomes being sampled in "sampling mode" .
--obs_file arg        self-defined measurement operation file string (if any).
--r arg (=32)         integer bit size.
--reorder arg (=1)    allow variable reordering or not.
                      0: disable reordering.
                      1: enable reordering (default option).
--alloc arg (=1)      allocate new BDDs when overflow is detected.
                      0: do not allocate new BDDs. This may lead to numerical errors.
                      1: allocate new BDDs (default option).

```
To use the sampling mode (default), it is required to have measurement operations included in the qasm file. Conversely, in all_amplitude mode, measurement operations are generally omitted, but if they are present in the qasm file, the final state vector will collapse based on the measurement result. It is important to note that all_amplitude mode is not recommended for simulations involving a large number of qubits, as it could result in a significantly long runtime.

For example, simulating `examples/bell_state_measure.qasm`, which is a 2-qubit bell state circuit with measurement gates at the end, with the sampling mode simulation option can be executed by
```commandline
./SliQSim --sim_qasm examples/bell_state_measure.qasm --type 0 --shots 1024
```

Then the sampled results will be shown:
```commandline
{ "counts": { "11": 542, "00": 482 } }
```

If option `--print_info` is used, simulation statistics such as runtime and memory usage will also be provided: 
```commandline
  Runtime: 0.014433 seconds
  Peak memory usage: 12611584 bytes
  #Applied gates: 2
  Max #nodes: 13
  Precision of integers: 32
  Accuracy loss: 2.22045e-16
```

To demonstrate the all_amplitude mode simulation, we use `examples/bell_state.qasm`, which is the same circuit as in the sampling mode simulation example except that the measurement gates are removed:
```commandline
./SliQSim --sim_qasm examples/bell_state.qasm --type 1
```

This will show the resulting state vector:
```commandline
{"statevector": ["0.707107", "0", "0", "0.707107"] }
```

`SliQSim` also supports various types of queries for user-defined properties. The users can use the expressions listed below to specify their desired properties.

* `bf {formula}`: It returns the probability that the measurement result satisfies the specified formula.
* `hweq/hwneq/hwgt/hwlt {qubits} {num}`: It returns the probability that the Hamming weight of specified qubits is equal to, non-equal to, greater than, lower than the specified number num.
* `inteq/intneq/intgt/intlt {qubits} {num}`: It returns the probability that the binary integer represented by the specified qubits is equal to, non-equal to, greater than, lower than the specified number num.
* `expt {qubits} {Pauli(-value)_string}`: It returns the expectation value of the observable defined by the specified Pauli string with optionally specified measured-values (0 or 1) for each individual Pauli operator over the specified qubits. Specifically, the returned value is calculated by the product of the probability of obtaining the specified measured values and the expectation value of the Pauli operators without measured values with respect to the post-measurement state.
* `weightedsum {weight} {expression} {...} endweighteedsum`: It returns the weighted sum of values of the specified expressions. Each weight and its corresponding expression are stated in an independent line, and the lines are clipped by keywords weightedsum and endweightedsum. The expression can be `bf`, `hweq`, `hwneq`, `hwgt`, `hwlt`, `inteq`, `intneq`, `intgt`, `intlt`, or `expt`.
* `assign {var_name} {expression}`: It does not return values, but rather stores the specified expression in the specified variable name var_name, which can be further utilized in `bf` expression. The expression can be `bf`, `hweq`, `hwneq`, `hwgt`, `hwlt`, `inteq`, `intneq`, `intgt`, or `intlt`.
* `between/outof/leq/geq {threshold}`: It can be specified before the expressions: `bf`, `hweq`, `hwneq`, `hwgt`, `hwlt`, `inteq`, `intneq`, `intgt`, `intlt`, `expt`, `weightedsum`. The predicate returns
true or false according to whether or not the probability returned by its subsequent expression function is between, out of, less than or equal to, or greater than or equal to the specified range.
* `amp {compt_basis}`: It returns the probability amplitude (as a complex number) of the specified computational basis compt_basis.
* `dist {qubits}`: It returns the exact spectrum of the probability distribution upon measuring the specified qubits.

To demonstrate the query mode simulation, we use `examples/grover_10.qasm` and `examples/demo.obs`:
```commandline
./SliQSim --sim_qasm examples/grover_10.qasm --obs_file examples/demo.obs --type 2
```

This will show the query results of each user-defined properties:
```commandline
"bf tmp_var^q[2] | tmp_var":
        0.541386
"hweq q[0] q[2] 1":
        0.5
"hwgt q[0] q[2] q[1] 2":
        0.0206929
"inteq q[0] q[1] q[2] 4":
        0.0206929
"intneq q[0] q[1] q[2] 5":
        0.979307
"expt q[0] q[1] q[2] zzz":
        0
"expt q[0] q[1] q[2] yxy":
        -5.55112e-17
"weightedsum
0.750000 expt q[0] q[1] q[2] zzz
0.150000 expt q[0] q[1] q[2] izi
-0.500000 expt q[0] q[1] q[2] zzi
endweightedsum":
        -0.125169
"0.5 <= hwgt q[0] q[2] q[1] 2 <= 0.8":
        false
"amp 1001001010":
        -0.012715
"dist q[0] q[1] q[2]":
        000: 0, 001: 0.000275976, 010: 0.437921, 011: 0.0206929, 100: 0, 101: 0.000275976, 110: 0.437921, 111: 0.0206929
```

One may also execute our simulator as a backend option of Qiskit through [SliQSim Qiskit Interface](https://github.com/NTU-ALComLab/SliQSim-Qiskit-Interface), which supports "sampling" and "all_amplitude" simulation options now.


## Citation
Please cite the following paper if you use our simulator for your research:

<summary>
  <a href="https://ieeexplore.ieee.org/document/9586191">Y.-H. Tsai, J.-H. R. Jiang, and C.-S. Jhang, “Bit-slicing the Hilbert space:  Scaling up accurate quantum circuit simulation,” in <em>Design Automation Conference (DAC)</em>, 2021, pp. 439–444.</a>
</summary>

```bibtex
@INPROCEEDINGS{9586191,
  author={Tsai, Yuan-Hung and Jiang, Jie-Hong R. and Jhang, Chiao-Shan},
  booktitle={Design Automation Conference (DAC)}, 
  title={Bit-Slicing the Hilbert Space: Scaling Up Accurate Quantum Circuit Simulation}, 
  year={2021},
  pages={439-444},
  doi={10.1109/DAC18074.2021.9586191}
}
```

## Contact
If you have any questions or suggestions, feel free to [create an issue](https://github.com/NTU-ALComLab/SliQSim/issues), or contact us through matthewyhtsai@gmail.com.