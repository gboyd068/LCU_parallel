import numpy as np
import stim


def clifford_transform_group_to_zs(operator_group, m_qubits):
    """Takes a commuting groups of oeprators and calculate the clifford transformation a stim.Tableau object
    that take those Paulis and transform them to Z operators on m_qubits"""

    groups, idxs = greedy_search(operator_group, check_pauli_linear_independence, m_qubits)

    tableaus = []
    for group in groups:
        # I think the fact that redundant stabilisers are ignored means we cant use this
        # might have to figure out how to use from_conjugated_generators or just use the matirx form...
        tableaus.append(stim.Tableau.from_stabilizers(stabilizers=group,
                allow_underconstrained=True, allow_redundant=False))
        
    group_sizes = [len(group) for group in groups]
    return tableaus, group_sizes


def clifford_transform_multiple_groups_to_zs(operator_groups, m_qubits):
    """Takes a list of commuting groups of oeprators and calculate the clifford transformation a stim.Tableau object
    that take those Paulis and transform them to Z operators on m_qubits"""
    tableaus = []
    group_sizes = []
    for operator_group in operator_groups:
        tabs, gsizes = clifford_transform_group_to_zs(operator_group, m_qubits)
        tableaus +=tabs
        group_sizes += gsizes
    return tableaus, group_sizes


def gaussian_elimination_binary(matrix):
    """tested once, should maybe test some more just in case"""
    m = len(matrix[0])  # Number of columns
    n = len(matrix)     # Number of rows

    for j in range(m):
        row_i = -1
        # Search for 1 in column j
        for i in range(n):
            if matrix[i][j] == 1:
                row_i = i
                break
        
        if row_i != -1:
            # Mark row i (optional, depending on what you want to do with marked rows)
            
            # Add column j to column k where A_ik is 1 and k != j
            for k in range(m):
                if k != j and matrix[row_i][k] == 1:
                    for i in range(n):
                        matrix[i][k] ^= matrix[i][j]  # XOR operation
    return matrix



def generate_symplectic(group):
    """Takes a list of Pauli operators and returns a symplectic matrix for the transformation that takes Z operators to those paulis"""
    AB = np.stack(list(map(lambda x: x.to_numpy(), group)))

    A = AB[:, :AB.shape[1]//2]
    B = AB[:, AB.shape[1]//2:]
    print(A)
    print(B)
    
    # find null vectors of A via gaussian elimination
    At = A.transpose().reshape(A.shape[2], A.shape[0])
    Atpad = np.array(np.concatenate((At, np.zeros((At.shape[0], At.shape[0]-At.shape[1]))), axis=1), dtype=bool)
    print(gaussian_elimination_binary(Atpad))

# gaussian_elimination_binary(np.array([[1,1,0,0],[1,1,0,1,],[0,1,1,1,],[0,0,1,0],[0,0,0,1]]))

# generate_symplectic([stim.PauliString("XYX"), stim.PauliString("Z_Z")])
    

def greedy_search(operators, condition, m_qubits):
    groups = []
    remaining_operators = operators.copy()

    while remaining_operators:
        current_group = [remaining_operators.pop(0)]
        i = 0

        while i < len(remaining_operators):
            if condition(current_group + [remaining_operators[i]]):
                current_group.append(remaining_operators.pop(i))
            else:
                i += 1
            if len(current_group) == m_qubits:
                break

        groups.append(current_group)

    return groups


def check_pauli_linear_independence(pauli_list):
    """Takes a list of stim.PauliString objects and checks if they are linearly independent by attempting to call stim.Tableau.from_stabilizers on them and returning
    false if it raises an exception
    Doing it this way because doing this with stim is probably faster than any python implementation"""
    try:
        stim.Tableau.from_stabilizers(pauli_list, allow_underconstrained=True, allow_redundant=False)
        return True
    except ValueError as e:
        # should probably make sure that the error is not the anticommuting value error
        return False
    