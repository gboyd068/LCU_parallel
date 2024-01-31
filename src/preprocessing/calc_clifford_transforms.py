import stim


def clifford_transform_group_to_zs(operator_group, group_idxs, m_qubits):
    """Takes a commuting groups of oeprators and calculate the clifford transformation a stim.Tableau object
    that take those Paulis and transform them to Z operators on m_qubits"""

    groups, pauli_idxs = greedy_search(operator_group, group_idxs, check_pauli_linear_independence, m_qubits)

    tableaus = []
    for group in groups:
        # I think the fact that redundant stabilisers are ignored means we cant use this
        # might have to figure out how to use from_conjugated_generators or just use the matirx form...
        tableaus.append(stim.Tableau.from_stabilizers(stabilizers=group,
                allow_underconstrained=True, allow_redundant=False))
        
    return tableaus, pauli_idxs


def clifford_transform_multiple_groups_to_zs(operator_groups, list_group_idxs, m_qubits):
    """Takes a list of commuting groups of oeprators and calculate the clifford transformation a stim.Tableau object
    that take those Paulis and transform them to Z operators on m_qubits"""
    tableaus = []
    group_sizes = []
    for i, operator_group in enumerate(operator_groups):
        group_idxs = list_group_idxs[i]
        tabs, pauli_idxs = clifford_transform_group_to_zs(operator_group, group_idxs, m_qubits)
        tableaus += tabs
        group_sizes += pauli_idxs
    return tableaus, group_sizes
    

def greedy_search(operators, group_idxs, condition, m_qubits):
    groups = []
    idxs = []
    remaining_operators = operators.copy()
    remaining_idxs = list(group_idxs.copy())

    while remaining_operators:
        current_group = [remaining_operators.pop(0)]
        current_idxs = [remaining_idxs.pop(0)]
        i = 0

        while i < len(remaining_operators):
            if condition(current_group + [remaining_operators[i]]):
                current_group.append(remaining_operators.pop(i))
                current_idxs.append(remaining_idxs.pop(i))
            else:
                i += 1
            if len(current_group) == m_qubits:
                break

        groups.append(current_group)
        idxs.append(current_idxs)

    return groups, idxs


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
        

