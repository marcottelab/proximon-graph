import networkx as nx
import argparse


def filter_graph_by_weight(
    input_gml: str,
    output_gml: str,
    weight_key: str = "abs_weight",
    min_weight: float = 3.0,
):
    # Load graph
    G = nx.read_gml(input_gml)

    print(f"[loaded] nodes={G.number_of_nodes()} edges={G.number_of_edges()}")

    # Create new graph
    H = nx.DiGraph()

    # Add nodes
    H.add_nodes_from(G.nodes(data=True))

    kept_edges = 0

    for u, v, data in G.edges(data=True):
        weight = data.get(weight_key, 0)

        if weight >= min_weight:
            H.add_edge(u, v, **data)
            kept_edges += 1

    # Remove isolated nodes (important)
    H.remove_nodes_from(list(nx.isolates(H)))

    print(f"[filtered] kept_edges={kept_edges}")
    print(f"[result] nodes={H.number_of_nodes()} edges={H.number_of_edges()}")

    # Save
    nx.write_gml(H, output_gml)
    print(f"[wrote] {output_gml}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_gml", required=True)
    parser.add_argument("--output_gml", required=True)
    parser.add_argument("--min_weight", type=float, default=3.0)
    parser.add_argument(
        "--weight_key",
        default="abs_weight",
        choices=["abs_weight", "norm_weight"],
        help="Which edge weight to filter on",
    )

    args = parser.parse_args()

    filter_graph_by_weight(
        input_gml=args.input_gml,
        output_gml=args.output_gml,
        weight_key=args.weight_key,
        min_weight=args.min_weight,
    )


if __name__ == "__main__":
    main()
    
    
# EXAMPLE USAGE
# python filter_graph_by_weight.py \
#   --input_gml full_graph.gml \
#   --output_gml high_weight.gml \
#   --min_weight 3