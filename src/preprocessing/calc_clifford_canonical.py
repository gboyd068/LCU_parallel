import numpy as np
import stim
import galois
BinaryMatrix = galois.GF(2)
from sympy import Matrix as SympyMatrix

def get_blocks(sym_mat):
    """Return the blocks of a symplectic matrix."""
    n = sym_mat.shape[0]//2
    return BinaryMatrix(sym_mat[:n, :n]), BinaryMatrix(sym_mat[:n, n:]), BinaryMatrix(sym_mat[n:, :n]), BinaryMatrix(sym_mat[n:, n:])

def convert_tableau_to_symplectic(tableau):
    binmats = list(map(BinaryMatrix, map(lambda x: x.astype(int), tableau.to_numpy())))[:-2]
    return BinaryMatrix(np.block([[binmats[0], binmats[1]], [binmats[2], binmats[3]]]))

def convert_symplectic_to_tableau(sym_mat):
    n = sym_mat.shape[0]//2
    return stim.Tableau.from_numpy(np.block([[sym_mat[:n, :n], sym_mat[:n, n:]], [sym_mat[n:, :n], sym_mat[n:, n:]]]))

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

def partial_phase_sym(diagonal):
    """takes a diagonal matrix from a symplectic matrix and returns a symplectic matrix that is the phase layer corresponding to that diagonal"""
    n = diagonal.shape[0]
    return BinaryMatrix(np.block([[BinaryMatrix.Identity(n), BinaryMatrix.Zeros((n,n))], [diagonal, BinaryMatrix.Identity(n)]]))


def partial_hadamard_sym(sym_mat):
    """Return the symplectic matrix corresponding to a layer of Hadamard gates on the qubits required to make the top right block of sym_mat invertible"""
    n = sym_mat.shape[0]//2
    B = get_blocks(sym_mat)[1]
    rank_B = np.linalg.matrix_rank(B)

    Htab = stim.Tableau(n)
    H = stim.Tableau.from_named_gate("H")

    # find rank_B non-zero linearly indpendent rows of B
    li_idxs = linearly_independent_row_idxs(B)[:rank_B]
    for i in [x for x in range(n) if x not in li_idxs]:
        Htab.append(H, [i])
    return convert_tableau_to_symplectic(Htab)


def linearly_independent_row_idxs(matrix):
    """Return a list indices for the linearly independent rows of matrix"""
    _, li_idxs = SympyMatrix(matrix).T.rref()
    return li_idxs



def clifford_canonical_form(tableau):
    """Calculate the canonical form of a Clifford circuit given by a stim.Tableau object, using the method
    given in Proctor and Young https://arxiv.org/abs/2310.10882"""
    # WARNING this my not properly track tableau signs!!!
    n_qubits = len(tableau)
    # STEP 1, make the top right block of the symplectic matrix invertible using a hadamard layer on lhs
    mat = convert_tableau_to_symplectic(tableau)
    Hmat = partial_hadamard_sym(mat)
    mat = Hmat @ mat

    # STEP 2, apply a cnot layer on rhs to make top right block identity
    Q = get_blocks(mat)[1]
    CQT = cnot_sym(Q.transpose())
    mat = mat @ CQT

    # STEP 3, apply a phase layer to make the bottom right a product of two invertible matrices
    D2 = get_blocks(mat)[3]
    D, M = sym_plus_diag_to_invertible(D2)
    P1 = partial_phase_sym(D)
    mat = P1 @ mat

    # STEP 4, apply a cnot layer
    Ninv = M.transpose()
    CNinv = cnot_sym(Ninv)
    mat = CNinv @ mat

    # STEP 5, apply another cnot layer to get rid of factors of Ninv
    CNTinv = cnot_sym(Ninv.transpose())
    mat = mat @ CNTinv

    # STEP 6, apply a phase layer to make the bottom right block 0
    Pall = phase_all_sym(n_qubits)
    mat = Pall @ mat

    # STEP 7, apply an all hadamard layer
    Hall = hadamard_all_sym(n_qubits)
    mat = Hall @ mat

    # STEP 8, apply a phase layer to make bottom left a product of two invertible matrices
    A6 = get_blocks(mat)[2]
    D, MTinv = sym_plus_diag_to_invertible(A6)
    P2 = partial_phase_sym(D)
    mat = P2 @ mat

    # STEP 9, apply a cnot layer
    M = np.linalg.inv(MTinv.transpose())
    CM = cnot_sym(M)
    mat = mat @ CM

    # STEP 10, apply an all phase layer
    mat = mat @ Pall

    # STEP 11, cnot layer to reduce to identity!
    CMinv = cnot_sym(MTinv.transpose())
    mat = CMinv @ mat

    # get values to return

    # combined cnots layer matrix (eq 29)
    L = np.linalg.inv(Q.transpose()) @ np.linalg.inv(Ninv.transpose()) @ MTinv.transpose()
    M = np.linalg.inv(MTinv.transpose())
    N = np.linalg.inv(Ninv)

    CL = cnot_sym(L)
    CM = cnot_sym(M)
    CN = cnot_sym(N)
    return CL, Pall, CM, P2, Hall, Pall, CN, P1, Hmat, mat 








if __name__ == "__main__":
    n_qubits = 4 
    tab = stim.Tableau(n_qubits)
    S = stim.Tableau.from_named_gate("S")
    for i in range(3):
        tab.append(S, [i])
    mat = convert_tableau_to_symplectic(tab)


