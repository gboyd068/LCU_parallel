# LCU Parallel

This repository contains code accompanying the paper [Low-Overhead Parallelisation of LCU via Commuting Operators](https://arxiv.org/abs/2312.00696), describing implementation of parallelized versions of the Linear Combination of Unitaries (LCU) primitive. LCU is a powerful tool in quantum computing for simulating quantum systems and solving certain optimization problems.

The repo contains python code in `src/preprocessing`, designed to take a Hamiltonian expressed as a list of Pauli terms and express it as a series of parallelizable groups, along with Clifford transformations allowing terms within the groups to be performed in parallel. The preprocessing code also contains methods for taking a Clifford transformation expressed in tableau form and producing a canonical form which can be then performed in constant depth.

`src/qs` provides Q# code written for the purposes of resource estimation that implements components of the circuit for applying parallelized LCU on a quantum computer

The repo also contains jupyter notebooks and data relevant to the plots produced in the paper.

# Acknowledgements

The repo makes use of [OpenFermion](https://github.com/quantumlib/OpenFermion), [Stim](https://github.com/quantumlib/Stim) and [graph-tool](https://git.skewed.de/count0/graph-tool)

