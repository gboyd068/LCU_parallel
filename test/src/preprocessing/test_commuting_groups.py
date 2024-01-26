import openfermion as of
from openfermionpyscf import run_pyscf
import pytest
from src.preprocessing.commuting_groups import commuting_groups

@pytest.mark.parametrize("geometry", [[('H', (0., 1., 0.)), ('H', (0.,0. , 1.)), ('O', (0.,0. , 0.))],
                            [('Li', (0., 0., 0.)), ('H', (0.,0. , 1.))]])
def test_commuting_groups(geometry):
    """test that commuting_groups returns a list of commuting groups of operators"""
    # get the MolecularData object
    molecule = of.MolecularData(geometry, 'sto-3g', 1)

    # Set up the PySCF molecule object
    mol = run_pyscf(molecule)

    # Generate the molecular Hamiltonian
    hamiltonian = of.transforms.get_fermion_operator(
        mol.get_molecular_hamiltonian())
    n_qubits = mol.n_qubits

    hamiltonian = of.transforms.jordan_wigner(hamiltonian)

    operator_groups, group_idxs = commuting_groups(hamiltonian, n_qubits)

    zero_op = of.QubitOperator.zero()
    for group in operator_groups:
        for i, term1 in enumerate(group):
            for j, term2 in enumerate(group):
                if i==j:
                    continue
                assert of.utils.commutator(term1,term2) == zero_op
