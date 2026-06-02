import json
import numpy as np
import networkx as nx
import ot
from Bio.PDB import PDBParser

from scipy.optimize import   linear_sum_assignment
import numpy as np


# =========================================================
# 1. Build residue graph from PDB
# =========================================================
def pdb_to_residue_graph(pdb_path, threshold=8.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("enzyme", pdb_path)

    model = structure[0]
    chain = list(model.get_chains())[0]

    residues = [res for res in chain if res.has_id("CA")]

    G = nx.Graph()

    # add nodes
    for i, res in enumerate(residues):
        G.add_node(
            i,
            resname=res.get_resname(),
            resid=res.get_id()[1]
        )

    # add edges (CA-CA distance threshold)
    for i, res_i in enumerate(residues):
        ca_i = res_i["CA"].get_coord()

        for j in range(i + 1, len(residues)):
            ca_j = residues[j]["CA"].get_coord()

            if np.linalg.norm(ca_i - ca_j) <= threshold:
                G.add_edge(i, j)

    return G


# =========================================================
# 2. Amino acid encoding
# =========================================================
AMINO_ACIDS = {
    "ALA": 0, "ARG": 1, "ASN": 2, "ASP": 3,
    "CYS": 4, "GLN": 5, "GLU": 6, "GLY": 7,
    "HIS": 8, "ILE": 9, "LEU": 10, "LYS": 11,
    "MET": 12, "PHE": 13, "PRO": 14, "SER": 15,
    "THR": 16, "TRP": 17, "TYR": 18, "VAL": 19
}


# =========================================================
# 3. Graph distance matrix (for FGW)
# =========================================================
def graph_distance_matrix(G):
    D = nx.floyd_warshall_numpy(G)
    D = np.asarray(D, dtype=float)

    finite_max = np.max(D[np.isfinite(D)])
    D[np.isinf(D)] = 2 * finite_max

    return D
    
    
def approx_ged(G1, G2):

    nodes1 = list(G1.nodes(data=True))
    nodes2 = list(G2.nodes(data=True))

    n, m = len(nodes1), len(nodes2)

    # cost matrix (node substitutions)
    C = np.zeros((n, m))

    for i, (_, d1) in enumerate(nodes1):
        for j, (_, d2) in enumerate(nodes2):

            if d1["resname"] == d2["resname"]:
                C[i, j] = 0
            else:
                C[i, j] = 1

    # Hungarian matching
    row_ind, col_ind = linear_sum_assignment(C)

    node_cost = C[row_ind, col_ind].sum()

    # --- edge mismatch penalty ---
    # map nodes G1 -> G2
    mapping = {nodes1[i][0]: nodes2[j][0] for i, j in zip(row_ind, col_ind)}

    edge_cost = 0

    for (u, v) in G1.edges():
        if u in mapping and v in mapping:
            u2, v2 = mapping[u], mapping[v]

            if not G2.has_edge(u2, v2):
                edge_cost += 1
        else:
            edge_cost += 1

    return node_cost + edge_cost


# =========================================================
# 4. Node cost matrix (for FGW)
# =========================================================
def node_cost_matrix(G1, G2):
    nodes1 = list(G1.nodes(data=True))
    nodes2 = list(G2.nodes(data=True))

    n, m = len(nodes1), len(nodes2)
    C = np.zeros((n, m), dtype=float)

    for i, (_, d1) in enumerate(nodes1):
        a1 = AMINO_ACIDS.get(d1["resname"], -1)

        for j, (_, d2) in enumerate(nodes2):
            a2 = AMINO_ACIDS.get(d2["resname"], -2)

            C[i, j] = 0.0 if a1 == a2 else 1.0

    return C


# =========================================================
# 5. FGW computation
# =========================================================
def compute_fgw(G1, G2, alpha=0.5):

    D1 = graph_distance_matrix(G1)
    D2 = graph_distance_matrix(G2)

    C = node_cost_matrix(G1, G2)

    p = np.ones(len(G1)) / len(G1)
    q = np.ones(len(G2)) / len(G2)

    fgw_dist = ot.gromov.fused_gromov_wasserstein2(
        C,
        D1,
        D2,
        p,
        q,
        loss_fun="square_loss",
        alpha=alpha
    )

    return fgw_dist


# =========================================================
# 6. GED computation (NetworkX A* approximation)
# =========================================================
def compute_ged(G1, G2):

    def node_subst_cost(n1_attr, n2_attr):
        return 0 if n1_attr["resname"] == n2_attr["resname"] else 1

    def node_del_cost(n_attr):
        return 1

    def node_ins_cost(n_attr):
        return 1

    def edge_subst_cost(e1, e2):
        return 0

    def edge_del_cost(e):
        return 1

    def edge_ins_cost(e):
        return 1

    ged = nx.graph_edit_distance(
        G1,
        G2,
        node_subst_cost=node_subst_cost,
        node_del_cost=node_del_cost,
        node_ins_cost=node_ins_cost,
        edge_subst_cost=edge_subst_cost,
        edge_del_cost=edge_del_cost,
        edge_ins_cost=edge_ins_cost
    )

    return ged

# =========================================================
# 7. Main execution
# =========================================================
if __name__ == "__main__":

    graph1_path = "../New-enzymes/B1-501/enzyme_mimicB1-501EN0a/B1-501-enzyme_mimicB1-501EN0a-ranked_0.pdb"
    graph2_path = "../New-enzymes/B1-501/enzyme_mimicB1-501EN0a/B1-501-enzyme_mimicB1-501EN0a-ranked_2.pdb"

    # Build graphs once
    G1 = pdb_to_residue_graph(graph1_path)
    G2 = pdb_to_residue_graph(graph2_path)

    print("Nodes G1:", len(G1), "Nodes G2:", len(G2))

    # -------------------------
    # FGW similarity
    # -------------------------
    fgw_dist = compute_fgw(G1, G2, alpha=0.5)
    print("\nFused Gromov-Wasserstein distance:", fgw_dist)

    # -------------------------
    # GED similarity
    # -------------------------
    ged_dist = approx_ged(G1, G2)
    print("Graph Edit Distance (approx):", ged_dist)
