import stim
import galois
import numpy as np
from src.preprocessing.calc_clifford_canonical import convert_tableau_to_symplectic

BinaryMatrix = galois.GF(2)


def matrix_to_stabilizers(mat):
    """converts the binary matrix giving stabilzers to a list of stim.PauliString"""
    return [
        stim.PauliString.from_numpy(
            xs=np.array(row[: len(row) // 2], dtype=bool),
            zs=np.array(row[len(row) // 2 :], dtype=bool),
        )
        for row in mat
    ]


def split_stabs_to_x_and_z(stabilizers):
    """assuming CSS stabilizers, splits the stabilizers into x and z type
    stabilizers"""
    xs = []
    zs = []
    for stab in stabilizers:
        x, z = stab.to_numpy()
        xtype = x.any()
        ztype = z.any()
        if xtype and ztype:
            raise ValueError("Stabilizer is not CSS")
        if xtype:
            xs.append(stab)
        if ztype:
            zs.append(stab)
    return xs, zs


def complete_stabilizers(stabilizers):
    return split_stabs_to_x_and_z(matrix_to_stabilizers(convert_tableau_to_symplectic(stim.Tableau.from_stabilizers(stabilizers, allow_underconstrained=True))))


def x_resource_state_stabilizers(n):
    """returns the stabilizers for the x resource state, although we dont really need these since it 
    is so easy to produce the state"""
    stabs = []
    # get the x stabilizers
    for i in range(n):
        ps = stim.PauliString(2 * n)
        ps[i] = 1
        ps[i + n] = 1
        stabs.append(ps)

    return complete_stabilizers(stabs)


def z_resource_state_stabilizers(n, U):
    """returns the stabilizers for the z resource state, where U is a BinaryMatrix
    giving the cnot layer operation"""

    matrix = BinaryMatrix.Identity(n) @ np.linalg.inv(U.transpose()) + BinaryMatrix.Identity(n)

    stabs = []
    for i, row in enumerate(matrix):
        ps = stim.PauliString(2 * n)
        for j, val in enumerate(row):
            if val:
                ps[j] = 1
        stabs.append(ps)
    return complete_stabilizers(stabs)