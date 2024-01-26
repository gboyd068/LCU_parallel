from functools import reduce
import openfermion as of
import stim


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

if __name__ == "__main__":
    term = of.ops.QubitOperator('X0 Y1') + of.ops.QubitOperator('X0 Y2')
    print(qubitop_to_stim_pauli_strings(term, 3))