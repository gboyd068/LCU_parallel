import openfermion as of
from openfermionpyscf import run_pyscf
import stim
from src.preprocessing.commuting_groups import commuting_groups
from src.preprocessing.calc_clifford_transforms import clifford_transform_group_to_zs, clifford_transform_multiple_groups_to_zs
from src.utils import clifford_idx_from_pauli_index, qubitop_to_stim_pauli_strings
import pytest


# test_data = [([stim.PauliString("XX"), stim.PauliString("YY"), stim.PauliString("ZZ")], 2)]

# @pytest.mark.parametrize("operator_group, m_qubits", test_data)
# def test_clifford_transform_group_to_zs(operator_group, m_qubits):
#     tableaus, pauli_idxs = clifford_transform_group_to_zs(operator_group, m_qubits)

#     # check correct number of operators in total
#     assert sum([len(idxs) for idxs in pauli_idxs]) == len(operator_group)


def test_clifford_transform_multiple_groups_to_zs():
    geometry = [('H', (0., 1., 0.)), ('H', (0.,0. , 1.)), ('O', (0.,0. , 0.))]
    molecule = of.MolecularData(geometry, 'sto-3g', 1)

    # Set up the PySCF molecule object
    mol = run_pyscf(molecule)

    # Generate the molecular Hamiltonian
    hamiltonian = of.transforms.get_fermion_operator(
        mol.get_molecular_hamiltonian())
    n_qubits = mol.n_qubits

    hamiltonian = of.transforms.jordan_wigner(hamiltonian)
    operators = qubitop_to_stim_pauli_strings(hamiltonian, n_qubits)
    operator_groups, group_idxs = commuting_groups(hamiltonian, n_qubits)
    m_qubits = n_qubits
    tableaus, pauli_idxs = clifford_transform_multiple_groups_to_zs(operator_groups, group_idxs, m_qubits)

    # check that the tableaus transform the operators to Zs
    identity_tableau = stim.Tableau(m_qubits)
    for tableau_idx, ts in enumerate(pauli_idxs):
        for place_within_tableau, hamiltonian_idx in enumerate(ts):
            if tableaus[tableau_idx] == identity_tableau:
                continue
            assert tableaus[tableau_idx].z_output(place_within_tableau) == operators[hamiltonian_idx]
    

    # check there are the correct number of operators
    assert sum([len(idxs) for idxs in pauli_idxs]) == len(operators)
        