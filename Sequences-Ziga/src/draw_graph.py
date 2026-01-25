import json
import argparse
import networkx as nx
import matplotlib.pyplot as plt

INPUT=""


def load_graph(json_file):
    with open(json_file, "r") as f:
        data = json.load(f)

    G = nx.Graph()

    # add nodes
    for node in data["nodes"]:
        G.add_node(node["id"], label=node.get("label", node["id"]))

    # add edges
    for edge in data["links"]:
        G.add_edge(
            edge["source"],
            edge["target"],
            weight=edge.get("weight", 1.0)
        )

    return G


def draw_graph(G, input_file):
    pos = nx.spring_layout(G, seed=42)

    labels = nx.get_node_attributes(G, "label")
    edge_labels = nx.get_edge_attributes(G, "weight")

    pos = nx.spring_layout(G, seed=42)

    nx.draw(
        G,
        pos,
        with_labels=False,     # usually better for large enzyme graphs
        node_size=50,          # <<< much smaller nodes
        linewidths=0.2
    )
    
    plt.savefig( input_file.split("/")[-1].split(".")[0] + ".png")
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Draw graph from JSON")
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to input graph JSON file"
    )

    args = parser.parse_args()
    
    if INPUT == "":
       INPUT = args.input
       
    G = load_graph(INPUT)
    draw_graph(G, args.input)


if __name__ == "__main__":
    main()

