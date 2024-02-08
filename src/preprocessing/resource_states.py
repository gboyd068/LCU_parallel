import stim
import galois
import numpy as np

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


def cstage_state1(n_qubits):
    """Produces the stabilizers of the first resource
    state for the cnot stage, n_qubits is the number of qubits
    in the system register"""
    mat = BinaryMatrix.Zeros((2 * 4 * n_qubits, 4 * n_qubits))
    mat[:n_qubits, :n_qubits] = BinaryMatrix.Identity(n_qubits)
    mat[n_qubits : 2 * n_qubits, :n_qubits] = BinaryMatrix.Identity(n_qubits)
    mat[2 * n_qubits : 4 * n_qubits, 2 * n_qubits : 4 * n_qubits] = (
        BinaryMatrix.Identity(2 * n_qubits)
    )

    mat[4 * n_qubits : 5 * n_qubits, n_qubits : 2 * n_qubits] = BinaryMatrix.Identity(
        n_qubits
    )
    mat[5 * n_qubits : 6 * n_qubits, n_qubits : 2 * n_qubits] = BinaryMatrix.Identity(
        n_qubits
    )
    return matrix_to_stabilizers(mat.transpose())


def cstage_state2(n_qubits, M):
    """Produces the stabilizers of the 2nd resource
    state for the cnot stage, n_qubits is the number of qubits
    in the system register, M is the binary matrix that corresponds to
    the cnot layer"""
    mat = BinaryMatrix.Zeros((2 * 4 * n_qubits, 4 * n_qubits))
    mat[2 * n_qubits : 3 * n_qubits, 2 * n_qubits : 3 * n_qubits] = np.linalg.inv(
        M
    ) + BinaryMatrix.Identity(n_qubits)
    mat[6 * n_qubits : 7 * n_qubits, :n_qubits] = BinaryMatrix.Identity(n_qubits)
    mat[7 * n_qubits :, :n_qubits] = np.linalg.inv(
        M.transpose()
    ) + BinaryMatrix.Identity(n_qubits)
    mat[4 * n_qubits : 6 * n_qubits, 2 * n_qubits : 4 * n_qubits] = (
        BinaryMatrix.Identity(2 * n_qubits)
    )
    return matrix_to_stabilizers(mat.transpose())


def cstage_state3(n_qubits):
    """Produces the stabilizers of the 3rd resource
    state for the cnot stage, n_qubits is the number of qubits
    in the system register"""
    mat = BinaryMatrix.Zeros((2 * 4 * n_qubits, 4 * n_qubits))
    mat[:n_qubits, :n_qubits] = BinaryMatrix.Identity(n_qubits)
    mat[2 * n_qubits : 4 * n_qubits, 2 * n_qubits :] = BinaryMatrix.Identity(
        2 * n_qubits
    )

    mat[5 * n_qubits : 6 * n_qubits, n_qubits : 2 * n_qubits] = BinaryMatrix.Identity(
        n_qubits
    )
    return matrix_to_stabilizers(mat.transpose())


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


xs, zs = split_stabs_to_x_and_z(cstage_state1(10))
print(xs)
print(zs)
