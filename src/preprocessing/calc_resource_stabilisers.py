import stim

def clifford_transform_group_to_zs(operator_group, m_qubits):
    """Takes a commuting groups of oeprators and calculate the clifford transformation a stim.Tableau object
    that take those Paulis and transform them to Z operators on m_qubits"""
    # split the operator group into paritions of size m_qubits
    sized_groups = [operator_group[i:i+m_qubits] for i in range(0, len(operator_group), m_qubits)]
    tableaus = []
    for group in sized_groups:
        print(group)
        # I think the fact that redundant stabilisers are ignored means we cant use this
        # might have to figure out how to use from_conjugated_generators or just use the matirx form...
        tableaus.append(stim.Tableau.from_stabilizers(stabilizers=group,
                allow_underconstrained=True, allow_redundant=True))
        
    group_sizes = [len(group) for group in sized_groups]
    return tableaus, group_sizes

def clifford_transform_multiple_groups_to_zs(operator_groups, m_qubits):
    """Takes a list of commuting groups of oeprators and calculate the clifford transformation a stim.Tableau object
    that take those Paulis and transform them to Z operators on m_qubits"""
    tableaus = []
    group_sizes = []
    for operator_group in operator_groups:
        tabs, gsizes = clifford_transform_group_to_zs(operator_group, m_qubits)
        tableaus += tabs
        group_sizes += gsizes
    return tableaus, group_sizes

