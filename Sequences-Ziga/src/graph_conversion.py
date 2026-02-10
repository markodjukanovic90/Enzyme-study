import os
import json
import numpy as np
import networkx as nx
import pandas as pd
from Bio.PDB import PDBParser


## generate a residue graph 
def pdb_to_residue_graph(pdb_path, threshold=8.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("enzyme", pdb_path)

    model = structure[0]
    chain = list(model.get_chains())[0]
    
    #CA -- central residue
    residues = [res for res in chain if res.has_id("CA")]

    G = nx.Graph()
    for i, res in enumerate(residues):
        G.add_node(
            i,
            resname=res.get_resname(),
            resid=res.get_id()[1]
        )

    for i, res_i in enumerate(residues):
        ca_i = res_i["CA"].get_coord()
        for j in range(i + 1, len(residues)):
            ca_j = residues[j]["CA"].get_coord()
            # threshold (8A)
            if np.linalg.norm(ca_i - ca_j) <= threshold:
                G.add_edge(i, j)

    return G



def save_graph_json(G, out_path):
    data = nx.node_link_data(G)
    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)


def extract_graph_features(G, title):
    n = G.number_of_nodes()
    m = G.number_of_edges()

    features = {
        "graph_id": title,
        "num_nodes": n,
        "num_edges": m,
        "density": nx.density(G),
        "avg_degree": (2 * m / n) if n > 0 else 0,
        "min_degree": min(dict(G.degree()).values()) if n > 0 else 0,
        "max_degree": max(dict(G.degree()).values()) if n > 0 else 0,
        "avg_clustering": nx.average_clustering(G),
        "num_components": nx.number_connected_components(G),
    }

    # Degree centrality (cheap)
    deg_cent = nx.degree_centrality(G)
    features["avg_degree_centrality"] = np.mean(list(deg_cent.values())) if n > 0 else 0
    features["max_degree_centrality"] = max(deg_cent.values()) if n > 0 else 0

    # Path-based features only if connected
    if nx.is_connected(G) and n > 1:
        features["diameter"] = nx.diameter(G)
        features["avg_shortest_path"] = nx.average_shortest_path_length(G)
    else:
        features["diameter"] = None
        features["avg_shortest_path"] = None
    
    return features


#TODO: RF-B1-503_top_50_rmsd_pdbs, RF-B1-503_top_50_seq_score_pdbs
for number in range(10):
    ensyme_type=f"B1-501-379-0{number}" # zero move: 
    #ensyme_type="B1-502"
    pdb_dir = f"../New-enzymes/{ensyme_type}/"
    graph_dir = f"graphs_json/{ensyme_type}/"

    os.makedirs(graph_dir, exist_ok=True)

    feature_rows = []

    for root, _, files in os.walk(pdb_dir):
        for fname in files:
            if not fname.endswith(".pdb"):
                continue

            pdb_path = os.path.join(root, fname)

            # create a unique ID including subdir
            rel_path = os.path.relpath(pdb_path, pdb_dir)
            graph_id = os.path.splitext(rel_path)[0].replace(os.sep, "_")
            if graph_id.split("_")[-1] != "0":
                continue

            print(f"Processing {graph_id}")

            G = pdb_to_residue_graph(pdb_path)

            # save graph JSON (mirror structure flattened)
            json_path = os.path.join(graph_dir, f"{graph_id}.json")
            save_graph_json(G, json_path)

            # extract features
            feats = extract_graph_features(G, graph_id)
            feature_rows.append(feats)
        
    df = pd.DataFrame(feature_rows)
    df.to_csv(f"features/enzyme_graph_features {ensyme_type}.csv", index=False)

