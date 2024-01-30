from functools import reduce
import bisect
import openfermion as of
import stim
import numpy as np


def qubitop_to_stim_pauli_strings(qubit_operator, n_qubits):
    """Convert an OpenFermion QubitOperator to a list of stim PauliStrings."""
    pauli_strings = []
    for term in qubit_operator.terms:
        string = ["_"] * n_qubits
        for qubit, pauli in term:
            string[qubit] = pauli
        string = reduce(lambda a,b: a+b, string)
        pauli_strings.append(stim.PauliString(string))
    return pauli_strings

def clifford_idx_from_pauli_index(pidx, group_sizes):
    """Convert a pauli index to a clifford index"""
    cidx = bisect.bisect(np.cumsum(group_sizes), pidx)
    idx_within_group = int(pidx - np.sum(group_sizes[:cidx]))
    return cidx, idx_within_group


if __name__ == "__main__":
    # term = of.ops.QubitOperator('X0 Y1') + of.ops.QubitOperator('X0 Y2')
    # print(qubitop_to_stim_pauli_strings(term, 3))
    group_sizes = [10, 5,3,2]
    print(clifford_idx_from_pauli_index(18, group_sizes))