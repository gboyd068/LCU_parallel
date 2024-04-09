import numpy as np
import matplotlib.pyplot as plt
import openfermion as of
from openfermionpyscf import run_pyscf
import stim
from src.preprocessing.calc_clifford_transforms import clifford_transform_multiple_groups_to_zs
from src.preprocessing.commuting_groups import commuting_groups
from src.utils import qubitop_to_stim_pauli_strings
from openfermion.chem import geometry_from_pubchem
molecule_names = ["molecular hydrogen", "lithium hydride", "water", "ammonia", "methane", "oxygen", "molecular nitrogen", "ethane", "disodium", 
                  # "XVOFZWCCFLVFRR-UHFFFAOYSA-N" # CrO oxochromium
                  ]
geometries = [geometry_from_pubchem(name) for name in molecule_names]
moleculeData = []
for i, geometry in enumerate(geometries):
    basis = "sto-3g"
    molecule = of.MolecularData(geometry, basis, 1)
    mol = run_pyscf(molecule)
    moleculeData.append((mol.n_qubits, mol, molecule_names[i], basis))
for i, geometry in enumerate(geometries[:-3]):
    basis = "6-31g"
    molecule = of.MolecularData(geometry, basis, 1)
    mol = run_pyscf(molecule)
    moleculeData.append((mol.n_qubits, mol, molecule_names[i], basis))

hamiltonianData = []
for md in moleculeData:
    print(md[2], md[3])
    ham = of.transforms.jordan_wigner(of.transforms.get_fermion_operator(md[1].get_molecular_hamiltonian()))
    print(len(ham.terms))
    hamiltonianData.append((md[0], len(ham.terms) ,  ham, md[2], md[3]))
hamiltonianData = sorted(hamiltonianData, key=lambda x: x[1])
print(hamiltonianData)

#     molecules.append(mol)
# ns = [mol.n_qubits for mol in molecules]
# hamiltonians = [of.transforms.jordan_wigner(of.transforms.get_fermion_operator(mol.get_molecular_hamiltonian())) for mol in molecules]
# fillings = []
# for i, hamiltonian in enumerate(hamiltonians):
#     n = ns[i]
#     print(n, len(hamiltonian.terms))
#     operator_groups, group_idxs = commuting_groups(hamiltonian, n)
#     print("calculated operator groups")
#     tableaus, pauli_idxs = clifford_transform_multiple_groups_to_zs(operator_groups, group_idxs, n)
#     fillings.append(np.mean([len(idx) for idx in pauli_idxs]) / n)

# plt.scatter(ns, fillings)