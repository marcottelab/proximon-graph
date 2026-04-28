#!/usr/bin/env python3
"""Extract marker-containing components from a GML graph and export Cytoscape JSON."""

import argparse
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import networkx as nx
import pandas as pd


NODE_INT_FIELDS = {"cluster_size"}
NODE_FLOAT_FIELDS = {"size"}
EDGE_INT_FIELDS = {"abs_weight", "num_genomes", "weight_raw", "genome_support"}
EDGE_FLOAT_FIELDS = {"norm_weight", "weight_norm"}

GENE_NAME_KEYS = [
    "all_gene_names",
    "gene_names",
    "member_gene_names",
    "members_gene_names",
    "member_genes",
    "genes",
    "gene_name",
    "name",
]


def read_table(path, sep="\t"):
    """Load a delimited table, allowing either tab or whitespace separation."""
    if sep == "whitespace":
        return pd.read_csv(path, sep=r"\s+", engine="python")

    table = pd.read_csv(path, sep=sep)
    if len(table.columns) == 1 and sep == "\t":
        return pd.read_csv(path, sep=r"\s+", engine="python")
    return table


def normalize_value(value):
    """Convert a table or graph value to a clean string."""
    if value is None:
        return ""
    if pd.isna(value):
        return ""
    return str(value).strip()


def split_maybe_listlike(value, split_pattern=r"\s*[;,|]\s*"):
    """Split list-like graph attributes stored as strings."""
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        out = []
        for item in value:
            out.extend(split_maybe_listlike(item, split_pattern=split_pattern))
        return out

    text = normalize_value(value)
    if not text:
        return []

    if re.search(split_pattern, text):
        return [part.strip() for part in re.split(split_pattern, text) if part.strip()]
    return [text]


def split_members(value):
    """Split the members attribute into raw protein ids."""
    # Member lists are semicolon-delimited, but individual member ids may
    # legitimately contain "|" as a namespace separator.
    return split_maybe_listlike(value, split_pattern=r"\s*[;,]\s*")


def unique_join(values):
    """Return a deterministic semicolon-separated string of unique values."""
    seen = set()
    ordered = []
    for value in values:
        cleaned = normalize_value(value)
        if not cleaned:
            continue
        key = cleaned.lower()
        if key in seen:
            continue
        seen.add(key)
        ordered.append(cleaned)
    return ";".join(ordered)


def collect_gene_names(data):
    """Gather candidate gene names from common node attributes."""
    names = []
    for key in GENE_NAME_KEYS:
        names.extend(split_maybe_listlike(data.get(key)))
    return [name for name in names if normalize_value(name)]


def choose_preferred_label(node_id, data):
    """Pick a readable fallback label for non-marker context nodes."""
    gene_names = collect_gene_names(data)
    for gene_name in sorted(gene_names, key=lambda item: (len(item), item.lower())):
        if gene_name.lower().startswith(("lac", "trp")):
            return gene_name

    for key in ("gene_name", "name", "product", "top_product", "annotation"):
        value = normalize_value(data.get(key))
        if value:
            return value
    return str(node_id)


def _coerce_int_maybe(value):
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, str):
        text = value.strip()
        if re.fullmatch(r"[+-]?\d+", text):
            return int(text)
        if re.fullmatch(r"[+-]?\d+\.0+", text):
            return int(float(text))
    return value


def _coerce_float_maybe(value):
    if isinstance(value, bool):
        return float(value)
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        text = value.strip()
        try:
            return float(text)
        except ValueError:
            return value
    return value


def coerce_graph_numeric_attributes(graph):
    """Coerce known numeric attributes after GML import."""
    for _, data in graph.nodes(data=True):
        for key in list(data.keys()):
            if key in NODE_INT_FIELDS:
                data[key] = _coerce_int_maybe(data[key])
            elif key in NODE_FLOAT_FIELDS:
                data[key] = _coerce_float_maybe(data[key])

    for _, _, data in graph.edges(data=True):
        for key in list(data.keys()):
            if key in EDGE_INT_FIELDS:
                data[key] = _coerce_int_maybe(data[key])
            elif key in EDGE_FLOAT_FIELDS:
                data[key] = _coerce_float_maybe(data[key])

    return graph


def load_marker_rows(path, id_column, sep="\t"):
    """Read the marker table and index rows by raw_id-like identifier."""
    table = read_table(path, sep=sep)
    required = {id_column, "label", "role"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(
            f"marker table missing columns: {sorted(missing)}. Found: {list(table.columns)}"
        )

    rows_by_id = defaultdict(list)
    columns = list(table.columns)
    for _, row in table.iterrows():
        marker_id = normalize_value(row[id_column])
        if not marker_id:
            continue

        cleaned_row = {}
        for column in columns:
            value = normalize_value(row[column])
            if value:
                cleaned_row[column] = value
        rows_by_id[marker_id].append(cleaned_row)

    return rows_by_id, columns


def build_graph_index(graph, members_attr="members"):
    """Map raw protein ids to graph node ids using node ids and member lists."""
    protein_to_nodes = defaultdict(set)
    for node_id, data in graph.nodes(data=True):
        protein_to_nodes[normalize_value(node_id)].add(node_id)
        for member in split_members(data.get(members_attr, "")):
            protein_to_nodes[member].add(node_id)
    return protein_to_nodes


def summarize_column_hits(rows, column):
    """Collapse column values from multiple matching marker rows."""
    return unique_join(row.get(column, "") for row in rows)


def tag_marker_nodes(graph, rows_by_id, query_columns, members_attr="members", prefix="marker"):
    """Annotate matched nodes and return matched seed nodes plus unmatched raw ids."""
    protein_to_nodes = build_graph_index(graph, members_attr=members_attr)
    matched_ids_by_node = defaultdict(set)
    matched_rows_by_node = defaultdict(list)
    unmatched_ids = []

    for protein_id, rows in rows_by_id.items():
        matched_nodes = protein_to_nodes.get(protein_id, set())
        if not matched_nodes:
            unmatched_ids.append(protein_id)
            continue

        for node_id in matched_nodes:
            matched_ids_by_node[node_id].add(protein_id)
            matched_rows_by_node[node_id].extend(rows)

    for node_id, data in graph.nodes(data=True):
        matched_ids = sorted(matched_ids_by_node.get(node_id, set()))
        matched_rows = matched_rows_by_node.get(node_id, [])

        data[f"{prefix}_hit"] = 1 if matched_ids else 0
        data[f"{prefix}_seed"] = 1 if matched_ids else 0
        data[f"{prefix}_raw_ids"] = ";".join(matched_ids)
        data[f"{prefix}_count"] = len(matched_ids)

        for column in query_columns:
            column_attr = f"{prefix}_{column}"
            data[column_attr] = summarize_column_hits(matched_rows, column) if matched_rows else ""

        labels = data.get(f"{prefix}_label", "")
        roles = data.get(f"{prefix}_role", "")
        data[f"{prefix}_display_label"] = labels if labels else choose_preferred_label(node_id, data)
        data[f"{prefix}_display_role"] = roles
        data[f"{prefix}_class"] = "marker_hit" if matched_ids else "context"

    seed_nodes = {
        node_id
        for node_id, data in graph.nodes(data=True)
        if int(data.get(f"{prefix}_hit", 0)) == 1
    }
    return seed_nodes, sorted(unmatched_ids)


def iter_components(graph):
    """Yield connected or weakly connected components, largest first."""
    if graph.is_directed():
        components = list(nx.weakly_connected_components(graph))
    else:
        components = list(nx.connected_components(graph))
    components.sort(key=len, reverse=True)
    return components


def extract_components_for_seeds(graph, seed_nodes):
    """Return the subgraph containing every component that has at least one seed node."""
    if not seed_nodes:
        return graph.__class__()

    selected_nodes = set()
    component_lookup = {}
    component_sizes = {}
    for component_index, component_nodes in enumerate(iter_components(graph), start=1):
        component_nodes = set(component_nodes)
        if component_nodes.intersection(seed_nodes):
            selected_nodes.update(component_nodes)
            component_sizes[component_index] = len(component_nodes)
            for node_id in component_nodes:
                component_lookup[node_id] = component_index

    subgraph = graph.subgraph(selected_nodes).copy()
    for node_id, data in subgraph.nodes(data=True):
        if node_id in component_lookup:
            component_id = component_lookup[node_id]
            data["marker_component_id"] = component_id
            data["marker_component_size"] = component_sizes.get(component_id, 0)
        else:
            data["marker_component_id"] = 0
            data["marker_component_size"] = 0
    return subgraph


def extract_n_hop_subgraph(graph, seed_nodes, hops=1, directed="either"):
    """Return the graph induced by nodes within an n-hop neighborhood of the seed nodes."""
    if not seed_nodes:
        return graph.__class__()

    if hops <= 0:
        return graph.subgraph(seed_nodes).copy()

    if directed not in {"out", "in", "either"}:
        raise ValueError("directed must be one of: out, in, either")

    if directed == "either":
        neighborhood_graph = graph.to_undirected(as_view=True) if graph.is_directed() else graph
        neighborhood = set(seed_nodes)
        frontier = set(seed_nodes)
        for _ in range(hops):
            next_frontier = set()
            for node_id in frontier:
                next_frontier.update(neighborhood_graph.neighbors(node_id))
            next_frontier -= neighborhood
            neighborhood |= next_frontier
            frontier = next_frontier
        return graph.subgraph(neighborhood).copy()

    neighborhood = set(seed_nodes)
    frontier = set(seed_nodes)
    for _ in range(hops):
        next_frontier = set()
        for node_id in frontier:
            if directed == "out":
                next_frontier.update(target for _, target in graph.out_edges(node_id))
            else:
                next_frontier.update(source for source, _ in graph.in_edges(node_id))
        next_frontier -= neighborhood
        neighborhood |= next_frontier
        frontier = next_frontier

    return graph.subgraph(neighborhood).copy()


def annotate_selected_seed_distances(graph, seed_nodes, marker_prefix="marker", directed="either"):
    """Annotate nodes in a selected subgraph with hop distance to the nearest seed node."""
    distance_attr = f"{marker_prefix}_seed_distance"
    is_seed_attr = f"{marker_prefix}_seed"

    if graph.number_of_nodes() == 0:
        return graph

    selected_seeds = set(seed_nodes).intersection(graph.nodes())
    if not selected_seeds:
        for _, data in graph.nodes(data=True):
            data[distance_attr] = -1
            data[is_seed_attr] = 0
        return graph

    if directed == "either" or not graph.is_directed():
        traversal_graph = graph.to_undirected(as_view=True) if graph.is_directed() else graph
    elif directed == "out":
        traversal_graph = graph
    elif directed == "in":
        traversal_graph = graph.reverse(copy=False)
    else:
        raise ValueError("directed must be one of: out, in, either")

    distances = nx.multi_source_shortest_path_length(traversal_graph, selected_seeds)
    for node_id, data in graph.nodes(data=True):
        data[distance_attr] = int(distances.get(node_id, -1))
        data[is_seed_attr] = 1 if node_id in selected_seeds else 0

    return graph


def select_subgraph(graph, seed_nodes, selection_mode="component", hops=1, directed="either"):
    """Select either full seed-containing components or an n-hop neighborhood."""
    if selection_mode == "component":
        return extract_components_for_seeds(graph, seed_nodes)
    if selection_mode == "hops":
        return extract_n_hop_subgraph(graph, seed_nodes, hops=hops, directed=directed)
    raise ValueError("selection_mode must be one of: component, hops")


def bbox(points):
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    return min(xs), min(ys), max(xs), max(ys)


def bbox_from_positions(positions):
    if not positions:
        return 0.0, 0.0, 0.0, 0.0
    return bbox(list(positions.values()))


def normalize_to_origin(positions):
    if not positions:
        return {}
    xs = [xy[0] for xy in positions.values()]
    ys = [xy[1] for xy in positions.values()]
    min_x, min_y = min(xs), min(ys)
    return {node_id: (xy[0] - min_x, xy[1] - min_y) for node_id, xy in positions.items()}


def translate_positions(positions, dx, dy):
    return {node_id: (xy[0] + dx, xy[1] + dy) for node_id, xy in positions.items()}


def scale_positions(positions, scale):
    return {node_id: (xy[0] * scale, xy[1] * scale) for node_id, xy in positions.items()}


def fit_to_aspect_ratio(positions, target_aspect=1.0):
    """Lightly rescale a packed layout so it is not too tall or too wide."""
    if not positions:
        return {}

    min_x, min_y, max_x, max_y = bbox_from_positions(positions)
    width = max(max_x - min_x, 1e-9)
    height = max(max_y - min_y, 1e-9)
    current_aspect = width / height
    if abs(current_aspect - target_aspect) < 0.05:
        return positions

    if current_aspect > target_aspect:
        factor_y = current_aspect / target_aspect
        return {node_id: (x, y * factor_y) for node_id, (x, y) in positions.items()}

    factor_x = target_aspect / current_aspect
    return {node_id: (x * factor_x, y) for node_id, (x, y) in positions.items()}


def component_edge_density(graph, component_nodes):
    """Compute component density for packing order."""
    component = set(component_nodes)
    node_count = len(component)
    if node_count <= 1:
        return 0.0

    subgraph = graph.subgraph(component)
    if graph.is_directed():
        possible_edges = node_count * (node_count - 1)
    else:
        possible_edges = node_count * (node_count - 1) / 2.0
    return subgraph.number_of_edges() / max(possible_edges, 1.0)


def pack_components(component_entries, graph, padding=220.0, target_aspect=1.2):
    """Pack component layouts left-to-right, top-to-bottom."""
    if not component_entries:
        return {}

    enriched = []
    for entry in component_entries:
        positions = normalize_to_origin(entry["pos"])
        min_x, min_y, max_x, max_y = bbox_from_positions(positions)
        width = max(max_x - min_x, 1.0)
        height = max(max_y - min_y, 1.0)
        nodes = list(entry["nodes"])
        enriched.append(
            {
                "pos": positions,
                "nodes": nodes,
                "node_count": len(nodes),
                "width": width,
                "height": height,
                "area": width * height,
                "density": component_edge_density(graph, nodes),
            }
        )

    enriched.sort(
        key=lambda item: (item["node_count"], item["density"], item["area"]),
        reverse=True,
    )

    total_area = sum((item["width"] + padding) * (item["height"] + padding) for item in enriched)
    max_width = max(item["width"] for item in enriched)
    wrap_width = max(
        max_width + padding,
        math.sqrt(max(total_area, 1.0) * max(target_aspect, 1e-3)) * 1.25,
    )

    merged = {}
    x_cursor = 0.0
    y_cursor = 0.0
    row_height = 0.0

    for item in enriched:
        step_width = item["width"] + padding
        if x_cursor > 0 and (x_cursor + step_width) > wrap_width:
            x_cursor = 0.0
            y_cursor += row_height + padding
            row_height = 0.0

        merged.update(translate_positions(item["pos"], x_cursor, y_cursor))
        x_cursor += step_width
        row_height = max(row_height, item["height"])

    merged = normalize_to_origin(merged)
    return normalize_to_origin(fit_to_aspect_ratio(merged, target_aspect=target_aspect))


def compute_component_positions(
    graph,
    seed=42,
    k_small=0.9,
    iter_small=300,
    scale=450.0,
    pair_scale=150.0,
    padding=220.0,
    target_aspect=1.2,
    large_component_min_nodes=1500,
    large_component_layout="sfdp",
    large_component_scale=1.0,
):
    """Lay out each connected component independently, then pack them compactly."""
    component_entries = []
    for component_nodes in iter_components(graph):
        subgraph = graph.subgraph(component_nodes).copy()
        layout_graph = subgraph.to_undirected() if subgraph.is_directed() else subgraph

        if subgraph.number_of_nodes() == 1:
            node_id = next(iter(subgraph.nodes()))
            pos = {node_id: (0.0, 0.0)}
        elif subgraph.number_of_nodes() == 2:
            node_a, node_b = list(subgraph.nodes())
            pos = scale_positions({node_a: (0.0, 0.0), node_b: (1.0, 0.0)}, pair_scale)
        elif (
            large_component_layout == "sfdp"
            and subgraph.number_of_nodes() >= large_component_min_nodes
        ):
            pos = compute_sfdp_positions(layout_graph)
            pos = normalize_to_origin(pos)
            if large_component_scale != 1.0:
                pos = scale_positions(pos, large_component_scale)
        else:
            pos = nx.spring_layout(layout_graph, seed=seed, k=k_small, iterations=iter_small)
            pos = scale_positions(normalize_to_origin(pos), scale)

        component_entries.append({"pos": pos, "nodes": list(component_nodes)})

    return pack_components(component_entries, graph, padding=padding, target_aspect=target_aspect)


def to_cytoscapejs_json(graph, positions, marker_prefix="marker"):
    """Export a Cytoscape.js-compatible JSON payload with explicit positions."""
    elements = {"nodes": [], "edges": []}

    for node_id, data in graph.nodes(data=True):
        marker_label = normalize_value(data.get(f"{marker_prefix}_label"))
        marker_role = normalize_value(data.get(f"{marker_prefix}_role"))
        preferred_label = choose_preferred_label(node_id, data)
        display_label = marker_label if marker_label else preferred_label

        node = {
            "data": {
                "id": str(node_id),
                "label": display_label,
                "display_label": display_label,
                "context_label": preferred_label,
                "marker_label": marker_label,
                "marker_role": marker_role,
                "shape": "ellipse",
            }
        }

        for key, value in data.items():
            if isinstance(value, (str, int, float, bool)) or value is None:
                node["data"][key] = value

        if "cluster_size" in data:
            node["data"]["cluster_size"] = _coerce_int_maybe(data["cluster_size"])
        if "size" in data:
            node["data"]["size"] = _coerce_float_maybe(data["size"])
        elif "cluster_size" in data:
            node["data"]["size"] = float(_coerce_int_maybe(data["cluster_size"]))

        if node_id in positions:
            x, y = positions[node_id]
            node["position"] = {"x": float(x), "y": float(y)}

        elements["nodes"].append(node)

    for source, target, data in graph.edges(data=True):
        edge = {"data": {"id": f"{source}__{target}", "source": str(source), "target": str(target)}}
        for key, value in data.items():
            if isinstance(value, (str, int, float, bool)) or value is None:
                edge["data"][key] = value
        elements["edges"].append(edge)

    return {"elements": elements}


def compute_sfdp_positions(graph):
    """Lay out a graph with Graphviz sfdp."""
    try:
        from networkx.drawing.nx_agraph import graphviz_layout
    except Exception as exc:
        raise RuntimeError(
            "Could not import pygraphviz / nx_agraph graphviz_layout."
        ) from exc

    positions = graphviz_layout(graph, prog="sfdp")
    return {node_id: (float(x), float(y)) for node_id, (x, y) in positions.items()}


def write_positions_tsv(path, positions):
    """Write node coordinates as a TSV for downstream inspection."""
    with Path(path).open("w", encoding="utf-8") as handle:
        handle.write("id\tx\ty\n")
        for node_id, (x, y) in positions.items():
            handle.write(f"{node_id}\t{x}\t{y}\n")


def write_unmatched_ids(path, unmatched_ids):
    """Write unmatched marker raw ids one per line."""
    with Path(path).open("w", encoding="utf-8") as handle:
        for marker_id in unmatched_ids:
            handle.write(f"{marker_id}\n")


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Extract either marker-containing connected components or an n-hop marker neighborhood "
            "from a GML graph, annotate marker nodes with label/role metadata, and export Cytoscape JSON."
        )
    )
    ap.add_argument("--gml", required=True, help="Input GML graph")
    ap.add_argument("--markers_tsv", required=True, help="Marker TSV with raw_id, label, and role columns")
    ap.add_argument("--out_prefix", required=True, help="Output prefix for .cyjs.json and .positions.tsv")
    ap.add_argument("--id_column", default="raw_id", help="Marker-table column containing raw protein ids")
    ap.add_argument("--sep", default="\t", help="Table separator. Use '\\t' for TSV or 'whitespace'")
    ap.add_argument("--members_attr", default="members", help="Node attribute listing member protein ids")
    ap.add_argument("--marker_prefix", default="marker", help="Prefix for added node attributes")
    ap.add_argument(
        "--selection_mode",
        choices=["component", "hops"],
        default="component",
        help="Extract whole connected components containing marker seeds, or an n-hop seed neighborhood",
    )
    ap.add_argument(
        "--hops",
        type=int,
        default=1,
        help="Neighborhood size when --selection_mode hops is used",
    )
    ap.add_argument(
        "--directed",
        choices=["either", "out", "in"],
        default="either",
        help="Neighborhood directionality for hop expansion; 'either' treats the graph as undirected",
    )
    ap.add_argument("--write_gml", action="store_true", help="Also write the extracted marker-component subgraph as GML")
    ap.add_argument("--write_tagged_full_gml", action="store_true", help="Also write the fully tagged input graph as GML")
    ap.add_argument("--write_unmatched_ids", action="store_true", help="Also write unmatched raw_ids as a text file")
    ap.add_argument("--seed", type=int, default=42, help="Random seed for component spring layouts")
    ap.add_argument("--k_small", type=float, default=0.9, help="spring_layout k for extracted components")
    ap.add_argument("--iter_small", type=int, default=300, help="spring_layout iterations for extracted components")
    ap.add_argument("--scale", type=float, default=450.0, help="Scale factor for component layouts")
    ap.add_argument("--pair_scale", type=float, default=150.0, help="Spacing for 2-node components")
    ap.add_argument("--padding", type=float, default=220.0, help="Padding between packed components")
    ap.add_argument("--target_aspect", type=float, default=1.2, help="Desired width/height ratio for packed layouts")
    ap.add_argument(
        "--large_component_min_nodes",
        type=int,
        default=1500,
        help="Use the large-component layout mode at or above this component size",
    )
    ap.add_argument(
        "--large_component_layout",
        choices=["spring", "sfdp"],
        default="sfdp",
        help="Layout mode for very large extracted components",
    )
    ap.add_argument(
        "--large_component_scale",
        type=float,
        default=1.0,
        help="Optional rescaling factor applied after large-component layout",
    )
    args = ap.parse_args()

    out_prefix = Path(args.out_prefix)
    graph = nx.read_gml(args.gml)
    graph = coerce_graph_numeric_attributes(graph)

    rows_by_id, query_columns = load_marker_rows(args.markers_tsv, args.id_column, sep=args.sep)
    seed_nodes, unmatched_ids = tag_marker_nodes(
        graph,
        rows_by_id,
        query_columns,
        members_attr=args.members_attr,
        prefix=args.marker_prefix,
    )

    if args.write_tagged_full_gml:
        tagged_path = out_prefix.with_name(out_prefix.name + "_tagged_full").with_suffix(".gml")
        nx.write_gml(graph, tagged_path)
        print(f"Wrote tagged full graph: {tagged_path}")

    if args.write_unmatched_ids:
        unmatched_path = out_prefix.with_name(out_prefix.name + "_unmatched_raw_ids").with_suffix(".txt")
        write_unmatched_ids(unmatched_path, unmatched_ids)
        print(f"Wrote unmatched raw ids: {unmatched_path}")

    if not seed_nodes:
        raise SystemExit(
            f"No matching marker nodes were found for {len(rows_by_id)} query ids. "
            f"Unmatched ids: {len(unmatched_ids)}"
        )

    subgraph = select_subgraph(
        graph,
        seed_nodes,
        selection_mode=args.selection_mode,
        hops=args.hops,
        directed=args.directed,
    )
    if subgraph.number_of_nodes() == 0:
        raise SystemExit("The extracted marker subgraph is empty.")

    annotate_selected_seed_distances(
        subgraph,
        seed_nodes,
        marker_prefix=args.marker_prefix,
        directed=args.directed,
    )

    positions = compute_component_positions(
        subgraph,
        seed=args.seed,
        k_small=args.k_small,
        iter_small=args.iter_small,
        scale=args.scale,
        pair_scale=args.pair_scale,
        padding=args.padding,
        target_aspect=args.target_aspect,
        large_component_min_nodes=args.large_component_min_nodes,
        large_component_layout=args.large_component_layout,
        large_component_scale=args.large_component_scale,
    )
    cyjs = to_cytoscapejs_json(subgraph, positions, marker_prefix=args.marker_prefix)

    json_path = out_prefix.with_suffix(".cyjs.json")
    positions_path = out_prefix.with_suffix(".positions.tsv")
    json_path.write_text(json.dumps(cyjs), encoding="utf-8")
    write_positions_tsv(positions_path, positions)
    print(f"Wrote: {json_path}")
    print(f"Wrote: {positions_path}")

    if args.write_gml:
        gml_path = out_prefix.with_suffix(".gml")
        nx.write_gml(subgraph, gml_path)
        print(f"Wrote: {gml_path}")

    hit_attr = f"{args.marker_prefix}_hit"
    component_ids = {data.get("marker_component_id", 0) for _, data in subgraph.nodes(data=True)}
    print(
        f"Matched {len(seed_nodes)} graph nodes from {len(rows_by_id)} marker ids; "
        f"unmatched ids={len(unmatched_ids)}"
    )
    print(
        f"Export subgraph has {subgraph.number_of_nodes()} nodes, {subgraph.number_of_edges()} edges, "
        f"and {sum(1 for _, data in subgraph.nodes(data=True) if int(data.get(hit_attr, 0)) == 1)} marked nodes"
    )
    if args.selection_mode == "component":
        print(f"Selected {len([cid for cid in component_ids if cid])} marker-containing components")
    else:
        print(
            f"Selected a {args.hops}-hop marker neighborhood "
            f"(direction={args.directed})"
        )


if __name__ == "__main__":
    main()
