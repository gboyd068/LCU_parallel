import numpy as np
import stim
import galois

BinaryMatrix = galois.GF(2)
from src.preprocessing.calc_clifford_canonical import (
    sym_plus_diag_to_invertible,
    partial_hadamard_sym,
    convert_tableau_to_symplectic,
    get_blocks,
    linearly_independent_row_idxs,
    clifford_canonical_form,
)


def test_sym_plus_diag_to_invertible():
    """test that the output of sym_plus_diag_to_invertible is correct for a
    series of random inputs"""
    n = 5
    for _ in range(20):
        A = BinaryMatrix.Random((n, n))
        A = A + A.transpose()
        D, M = sym_plus_diag_to_invertible(A)
        assert (M @ M.transpose() - D == A).all()
        assert (M @ np.linalg.inv(M) == BinaryMatrix.Identity(n)).all()


# test_cases =   [[[1,0,0],[1,0,0],[0,1,0]],
#                 ,[[0, 0, 0],[1, 0, 0],[0, 0, 1]]
#                 ,[[1, 0, 0], [1, 0, 0],[0, 0, 1]]
#                 ,[[1, 1, 0], [1, 0, 1], [0, 1, 1]]]


def test_partial_hadamard_sym():
    """test that the output of partial_hadamard_sym is correct for a series of
    inputs
    REDO this for some particular inputs"""
    n = 2
    for tableau in stim.Tableau.iter_all(n, unsigned=True):
        mat = convert_tableau_to_symplectic(tableau)
        blocks = get_blocks(mat)
        print("rank B:", np.linalg.matrix_rank(blocks[1]))
        print("rank D:", np.linalg.matrix_rank(blocks[3]))
        print(mat)
        Hmat = partial_hadamard_sym(mat)
        # assert top right block of product is invertible
        assert np.linalg.det(get_blocks(Hmat @ mat)[1])


def test_clifford_canonical_form():
    for tableau in stim.Tableau.iter_all(3, unsigned=True):
        print(clifford_canonical_form(tableau, test=True))
