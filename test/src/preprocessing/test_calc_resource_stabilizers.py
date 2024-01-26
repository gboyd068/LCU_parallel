import openfermion as of
from openfermionpyscf import run_pyscf
import stim
from src.preprocessing.commuting_groups import commuting_groups
from src.preprocessing.calc_resource_stabilisers import clifford_transform_group_to_zs, clifford_transform_multiple_groups_to_zs
import pytest


test_data = [([stim.PauliString("XX"), stim.PauliString("YY"), stim.PauliString("ZZ")], 2)]

@pytest.mark.parametrize("operator_group, m_qubits", test_data)
def test_clifford_transform_group_to_zs(operator_group, m_qubits):
    tableaus, group_sizes = clifford_transform_group_to_zs(operator_group, m_qubits)

    # check that the tableaus transform the operators to Zs
    for i, g in enumerate(operator_group):
        assert tableaus[i//m_qubits].z_output(i%m_qubits) == operator_group[i]
    
    # check correct number of operators in total
    assert sum(group_sizes) == len(operator_group)



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

    operator_groups, group_idxs = commuting_groups(hamiltonian, n_qubits)
    m_qubits = n_qubits
    tableaus, group_sizes = clifford_transform_multiple_groups_to_zs(operator_groups, m_qubits)


    # check that the tableaus transform the operators to Zs
    for i, group in enumerate(operator_groups):
        for j, g in enumerate(group):
            assert tableaus[i//m_qubits].z_output(i%m_qubits) == group[j]
    
    # check correct number of operators in total
    # assert sum(group_sizes) == len(operator_group)