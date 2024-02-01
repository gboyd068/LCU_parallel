import numpy as np
import stim
import galois
BinaryMatrix = galois.GF(2)
from sympy import Matrix

def get_blocks(sym_mat):
    """Return the blocks of a symplectic matrix."""
    n = sym_mat.shape[0]//2
    return BinaryMatrix(sym_mat[:n, :n]), BinaryMatrix(sym_mat[:n, n:]), BinaryMatrix(sym_mat[n:, :n]), BinaryMatrix(sym_mat[n:, n:])

def convert_tableau_to_symplectic(tableau):
    """cursed conversions"""
    binmats = list(map(BinaryMatrix, map(lambda x: x.astype(int), tableau.to_numpy())))[:-2]
    return BinaryMatrix(np.block([[binmats[0], binmats[1]], [binmats[2], binmats[3]]]))

def hadamard_all_sym(n_qubits):
    """Return the symplectic matrix corresponding to a layer of Hadamard gates on all qubits."""
    tab = stim.Tableau(n_qubits)
    H = stim.Tableau.from_named_gate("H")
    for i in range(n_qubits):
        tab.append(H, [i])
    return convert_tableau_to_symplectic(tab)


def phase_all_sym(n_qubits):
    """Return the symplectic matrix corresponding to a layer of phase gates on all qubits."""
    tab = stim.Tableau(n_qubits)
    myS = stim.Tableau(1)
    S = stim.Tableau.from_named_gate("S")
    H = stim.Tableau.from_named_gate("H")
    # need to do this gate sequence as there seems to be a difference in convention between S in stim and S in the paper
    myS.append(H, [0])
    myS.append(S, [0])
    myS.append(H, [0])
    for i in range(n_qubits):
        tab.append(myS, [i])
    return convert_tableau_to_symplectic(tab)


def cnot_sym(M):
    """Given an n x n binary matrix M corresponding to the action of a circuit of CNOT gates, return the symplectic matrix corresponding to the circuit,
    which is [[M, 0], [0, M.transpose()**-1]]"""
    n = M.shape[0]
    M = BinaryMatrix(M)
    if not np.linalg.det(M):
        raise ValueError("Input matrix cannot be converted into a CNOT circuit as it is not invertible.")
    return BinaryMatrix(np.block([[M, BinaryMatrix.Zeros((n,n))], [BinaryMatrix.Zeros((n, n)), np.linalg.inv(M.transpose())]]))
    

def sym_plus_diag_to_invertible(A):
    """Given a symmetric binary matrix A, calculate a diagonal matrix D and invertible matrix M such that A+D = M  M^{T}
    according to lemma 7 of 'Improved simulation of stabilizer circuits' https://arxiv.org/abs/quant-ph/0406196"""
    n = A.shape[0]
    M = BinaryMatrix.Identity(n)
    for i in range(n):
        for j in range(n):
            if i <= j:
                continue
            M[i,j] = A[i,j] - np.sum(M[i,:j] * M[j,:j])
    D = M @ M.transpose() - A
    return D, M

def partial_hadamard_sym(sym_mat):
    """Return the symplectic matrix corresponding to a layer of Hadamard gates on the qubits required to make the top right block of sym_mat invertible"""
    n = sym_mat.shape[0]//2
    B = get_blocks(sym_mat)[1]
    rank_B = np.linalg.matrix_rank(B)

    Htab = stim.Tableau(n)
    H = stim.Tableau.from_named_gate("H")

    # find rank_B non-zero linearly indpendent rows of B
    li_idxs = linearly_independent_row_idxs(B)
    print(li_idxs)
    print(B)
    for i in [x for x in range(n) if x not in li_idxs[:rank_B]]:
        Htab.append(H, [i])
    return convert_tableau_to_symplectic(Htab)


def linearly_independent_row_idxs(matrix):
    """Return a list indices for the linearly independent rows of matrix"""
    _, li_idxs = Matrix(matrix).T.rref()
    return li_idxs


if __name__ == "__main__":
    n_qubits = 4 
    tab = stim.Tableau(n_qubits)
    S = stim.Tableau.from_named_gate("S")
    for i in range(3):
        tab.append(S, [i])
    mat = convert_tableau_to_symplectic(tab)


