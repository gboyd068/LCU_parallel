import numpy as np
import stim
import galois
BinaryMatrix = galois.GF(2)


def convert_tableau_to_symplectics(tableau):
    """cursed conversions"""
    return list(map(BinaryMatrix, map(lambda x: x.astype(int), tableau.to_numpy())))


def hadamard_all_sym(n_qubits):
    """Return the symplectic matrix corresponding to a layer of Hadamard gates on all qubits."""
    tab = stim.Tableau(n_qubits)
    H = stim.Tableau.from_named_gate("H")
    for i in range(n_qubits):
        tab.append(H, [i])
    return convert_tableau_to_symplectics(tab)


def phase_all_sym(n_qubits):
    """Return the symplectic matrix corresponding to a layer of phase gates on all qubits."""
    tab = stim.Tableau(n_qubits)
    S = stim.Tableau.from_named_gate("S")
    for i in range(n_qubits):
        tab.append(S, [i])
    return convert_tableau_to_symplectics(tab)


def cnot_sym(M):
    """Given an n x n matrix M corresponding to the action of a circuit of CNOT gates, return the symplectic matrix corresponding to the circuit,
    which is [[M, 0], [0, M.transpose()**-1]]"""
    n = M.shape[0]
    M = BinaryMatrix(M)
    if not np.linalg.det(M):
        raise ValueError("Input matrix cannot be converted into a CNOT circuit as it is not invertible.")
    return BinaryMatrix(np.block([[M, BinaryMatrix.Zeros((n,n))], [BinaryMatrix.Zeros((n, n)), np.linalg.inv(M.transpose())]]))
    

def partial_phase_sym(sym_mat):
    """Return the symplectic matrix corresponding to a layer of phase gates on the qubits."""
    pass

def partial_hadamard_sym(sym_mat):
    """Return the symplectic matrix corresponding to a layer of Hadamard gates on the qubits."""
    pass


matrix = BinaryMatrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0 ,0], [0, 0, 0, 1]])
print(cnot_sym(matrix))
print(hadamard_all_sym(4))
print(phase_all_sym(4))