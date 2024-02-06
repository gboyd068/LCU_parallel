import openfermion as of
from openfermionpyscf import run_pyscf
import stim
import pytest
from src.preprocessing.commuting_groups import commuting_groups


@pytest.mark.parametrize(
    "geometry",
    [
        [("H", (0.0, 1.0, 0.0)), ("H", (0.0, 0.0, 1.0)), ("O", (0.0, 0.0, 0.0))],
        [("Li", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, 1.0))],
    ],
)
def test_commuting_groups(geometry):
    """test that commuting_groups returns a list of commuting groups of operators"""
    # get the MolecularData object
    molecule = of.MolecularData(geometry, "sto-3g", 1)

    # Set up the PySCF molecule object
    mol = run_pyscf(molecule)

    # Generate the molecular Hamiltonian
    hamiltonian = of.transforms.get_fermion_operator(mol.get_molecular_hamiltonian())
    n_qubits = mol.n_qubits

    hamiltonian = of.transforms.jordan_wigner(hamiltonian)

    operator_groups, group_idxs = commuting_groups(hamiltonian, n_qubits)

    for group in operator_groups:
        for i, term1 in enumerate(group):
            for j, term2 in enumerate(group):
                if i == j:
                    continue
                assert term1.commutes(term2)
