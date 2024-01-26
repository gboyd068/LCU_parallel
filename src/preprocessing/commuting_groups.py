import openfermion as of
from openfermionpyscf import run_pyscf
from graph_tool import Graph
from graph_tool.topology import sequential_vertex_coloring
from graph_tool.spectral import adjacency
import numpy as np
from scipy.sparse import csr_matrix
from src.utils import qubitop_to_stim_pauli_strings




def commutativity_graph(hamiltonian, n_qubits):
    """Return a graph where edges are between terms that commute."""
    pauli_strings = qubitop_to_stim_pauli_strings(hamiltonian, n_qubits)
    graph = Graph(len(pauli_strings), directed=False)
    edges = []
    for i, term1 in enumerate(pauli_strings):
        for j, term2 in enumerate(pauli_strings[:i]):
            if i==j or j > i:
                continue
            if term1.commutes(term2):
                edges.append((i, j))
    graph.add_edge_list(edges)
    return graph


def calc_min_cliques(graph):
    """Return a list of minimal cliques in a graph."""
    # Get the complement graph
    # REPLACE THIS WITH FASTER ADJACENCY MATRIX METHOD AND REMOVE NETWORKX 
    adj = np.asarray(adjacency(graph).todense() == 0)
    np.fill_diagonal(adj, False)
    comp_sparse = np.dstack(csr_matrix(adj).nonzero())[0]
    complement_graph = Graph(g=comp_sparse, directed=False)
    coloring = sequential_vertex_coloring(complement_graph)
    return coloring.a
    

def commuting_groups(hamiltonian, n_qubits=None):
    """Return a list of commuting groups of paulis in the form of stim PauliStrings, as well as the indexes of terms in the Hamiltonian that belong to each group."""
    if n_qubits is None:
        n_qubits = hamiltonian.n_qubits
    # Get the minimal cliques
    clique_array = calc_min_cliques(commutativity_graph(hamiltonian, n_qubits))
    # count number of appearances of each integer in array
    clique_idxs, clique_counts = np.unique(clique_array, return_counts=True)
    # get the commuting groups
    pauli_strings = qubitop_to_stim_pauli_strings(hamiltonian, n_qubits)
    group_idxs = [np.argwhere(clique_array == idx).flatten() for idx in clique_idxs]
    operator_groups = []
    for i, group in enumerate(group_idxs):
        operator_groups.append([pauli_strings[idx] for idx in group])
    return operator_groups, group_idxs