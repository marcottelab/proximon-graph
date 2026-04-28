#!/usr/bin/env python3
import argparse
import json
import math
import re
from math import ceil, sqrt
from pathlib import Path

import networkx as nx


# ---------------- geometry helpers ----------------
def bbox(pts):
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    return (min(xs), min(ys), max(xs), max(ys))


def bbox_from_positions(pos):
    if not pos:
        return (0.0, 0.0, 0.0, 0.0)
    return bbox(list(pos.values()))


def translate(pos, dx, dy):
    return {n: (xy[0] + dx, xy[1] + dy) for n, xy in pos.items()}


def normalize_to_origin(pos):
    if not pos:
        return {}
    xs = [xy[0] for xy in pos.values()]
    ys = [xy[1] for xy in pos.values()]
    minx, miny = min(xs), min(ys)
    return {n: (xy[0] - minx, xy[1] - miny) for n, xy in pos.items()}


def scale_positions(pos, scale):
    return {n: (xy[0] * scale, xy[1] * scale) for n, xy in pos.items()}


def fit_to_aspect_ratio(pos, target_aspect=1.0):
    """
    Light-touch rescaling so the final drawing is not wildly wide or tall.
    target_aspect = width / height.
    """
    if not pos:
        return {}
    x0, y0, x1, y1 = bbox_from_positions(pos)
    w = max(x1 - x0, 1e-9)
    h = max(y1 - y0, 1e-9)
    cur = w / h
    if abs(cur - target_aspect) < 0.05:
        return pos

    if cur > target_aspect:
        # too wide: stretch y
        factor_y = cur / target_aspect
        return {n: (x, y * factor_y) for n, (x, y) in pos.items()}
    else:
        # too tall: stretch x
        factor_x = target_aspect / cur
        return {n: (x * factor_x, y) for n, (x, y) in pos.items()}




# ---------------- attribute coercion helpers ----------------
NODE_INT_FIELDS = {"cluster_size"}
NODE_FLOAT_FIELDS = {"size"}
EDGE_INT_FIELDS = {"abs_weight", "num_genomes", "weight_raw", "genome_support"}
EDGE_FLOAT_FIELDS = {"norm_weight", "weight_norm"}


def _coerce_int_maybe(v):
    if isinstance(v, bool):
        return int(v)
    if isinstance(v, int):
        return v
    if isinstance(v, float) and v.is_integer():
        return int(v)
    if isinstance(v, str):
        s = v.strip()
        if re.fullmatch(r"[+-]?\d+", s):
            try:
                return int(s)
            except Exception:
                return v
        if re.fullmatch(r"[+-]?\d+\.0+", s):
            try:
                return int(float(s))
            except Exception:
                return v
    return v


def _coerce_float_maybe(v):
    if isinstance(v, bool):
        return float(v)
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, str):
        s = v.strip()
        try:
            return float(s)
        except Exception:
            return v
    return v


def coerce_graph_numeric_attributes(G):
    for n, data in G.nodes(data=True):
        for k in list(data.keys()):
            if k in NODE_INT_FIELDS:
                data[k] = _coerce_int_maybe(data[k])
            elif k in NODE_FLOAT_FIELDS:
                data[k] = _coerce_float_maybe(data[k])

    for u, v, data in G.edges(data=True):
        for k in list(data.keys()):
            if k in EDGE_INT_FIELDS:
                data[k] = _coerce_int_maybe(data[k])
            elif k in EDGE_FLOAT_FIELDS:
                data[k] = _coerce_float_maybe(data[k])
    return G

# ---------------- parsing / labeling helpers ----------------
def _split_maybe_listlike(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        out = []
        for v in value:
            out.extend(_split_maybe_listlike(v))
        return out
    s = str(value).strip()
    if not s:
        return []
    # common separators seen in exported attrs
    parts = re.split(r"\s*[;,|]\s*", s)
    return [p for p in parts if p]


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


def collect_gene_names(data):
    names = []
    for key in GENE_NAME_KEYS:
        if key in data:
            names.extend(_split_maybe_listlike(data.get(key)))
    # de-dup but preserve order
    seen = set()
    out = []
    for n in names:
        n2 = n.strip()
        if not n2:
            continue
        low = n2.lower()
        if low in seen:
            continue
        seen.add(low)
        out.append(n2)
    return out


OPERON_PREFIXES = ("lac", "trp")


def choose_preferred_label(node_id, data):
    gene_names = collect_gene_names(data)
    gene_names_sorted = sorted(gene_names, key=lambda x: (len(x), x.lower()))

    # 1) preferentially use lac* / trp* when present among member gene names
    for prefix in OPERON_PREFIXES:
        pref = [g for g in gene_names_sorted if g.lower().startswith(prefix)]
        if pref:
            return pref[0], gene_names

    # 2) otherwise prefer an explicit gene_name if present
    for key in ("gene_name", "name"):
        val = str(data.get(key, "")).strip()
        if val:
            return val, gene_names

    # 3) otherwise use the most informative product-like label available
    for key in ("product", "top_product", "annotation"):
        val = str(data.get(key, "")).strip()
        if val:
            return val, gene_names

    return str(node_id), gene_names


# ---------------- packing ----------------
def component_edge_density(G, comp_nodes):
    comp = set(comp_nodes)
    n = len(comp)
    if n <= 1:
        return 0.0
    sub = G.subgraph(comp)
    if G.is_directed():
        possible = n * (n - 1)
    else:
        possible = n * (n - 1) / 2.0
    return sub.number_of_edges() / max(possible, 1.0)


def pack_components_by_size_bands(component_entries, G, padding=200.0, target_aspect=1.0):
    """
    Arrange components in a rectangle-filling shelf layout while preserving the
    semantic ordering the user wants:
    - larger connected components first
    - denser components earlier within the same size class
    - fill rows left-to-right under a common wrap width

    This avoids the large staircase gaps caused by strict one-band-per-size
    placement, but still keeps the biggest / most complex structures toward the
    top-left of the figure.
    """
    if not component_entries:
        return {}

    enriched = []
    for entry in component_entries:
        pos = normalize_to_origin(entry["pos"])
        nodes = list(entry["nodes"])
        x0, y0, x1, y1 = bbox(list(pos.values()))
        w = max(x1 - x0, 1.0)
        h = max(y1 - y0, 1.0)
        enriched.append({
            "pos": pos,
            "nodes": nodes,
            "n_nodes": len(nodes),
            "w": w,
            "h": h,
            "density": component_edge_density(G, nodes),
            "area": w * h,
        })

    size_groups = {}
    for item in enriched:
        size_groups.setdefault(item["n_nodes"], []).append(item)

    # preserve grouping by component size for ordering, but flatten into one
    # rectangle-filling sequence so the overall figure becomes compact.
    ordered = []
    for size in sorted(size_groups.keys(), reverse=True):
        items = size_groups[size]
        items.sort(key=lambda d: (d["density"], d["area"], d["w"], d["h"]), reverse=True)
        ordered.extend(items)

    total_area = sum((it["w"] + padding) * (it["h"] + padding) for it in ordered)
    max_item_w = max(it["w"] for it in ordered)
    wrap_w = max(
        max_item_w + padding,
        math.sqrt(max(total_area, 1.0) * max(target_aspect, 1e-3)) * 1.25,
    )

    merged = {}
    x_cursor = 0.0
    y_cursor = 0.0
    row_h = 0.0

    for item in ordered:
        iw = item["w"]
        ih = item["h"]
        step_w = iw + padding

        if x_cursor > 0 and (x_cursor + step_w) > wrap_w:
            x_cursor = 0.0
            y_cursor += row_h + padding
            row_h = 0.0

        placed = translate(item["pos"], x_cursor, y_cursor)
        merged.update(placed)

        x_cursor += step_w
        row_h = max(row_h, ih)

    return merged


# ---------------- cytoscape.js export ----------------
def to_cytoscapejs_json(G, positions):
    elements = {"nodes": [], "edges": []}

    for n, data in G.nodes(data=True):
        node = {"data": {"id": str(n)}}
        for k, v in data.items():
            if isinstance(v, (str, int, float, bool)) or v is None:
                node["data"][k] = v

        preferred_label, gene_names = choose_preferred_label(n, data)
        node["data"]["label"] = preferred_label
        node["data"]["display_label"] = preferred_label
        node["data"]["shape"] = "ellipse"  # Cytoscape.js / Cytoscape both accept ellipse-like circle styling
        node["data"]["all_gene_names"] = ";".join(gene_names)

        if "cluster_size" in data:
            node["data"]["cluster_size"] = _coerce_int_maybe(data["cluster_size"])
        if "size" in data:
            node["data"]["size"] = _coerce_float_maybe(data["size"])
        elif "cluster_size" in data:
            node["data"]["size"] = float(_coerce_int_maybe(data["cluster_size"]))

        if n in positions:
            x, y = positions[n]
            node["position"] = {"x": float(x), "y": float(y)}

        elements["nodes"].append(node)

    for u, v, data in G.edges(data=True):
        e = {"data": {"id": f"{u}__{v}", "source": str(u), "target": str(v)}}
        for k, val in data.items():
            if isinstance(val, (str, int, float, bool)) or val is None:
                e["data"][k] = val
        elements["edges"].append(e)

    return {"elements": elements}


def write_positions_tsv(path: Path, positions: dict):
    with path.open("w", encoding="utf-8") as f:
        f.write("id\tx\ty\n")
        for n, (x, y) in positions.items():
            f.write(f"{n}\t{x}\t{y}\n")


# ---------------- directed-aware component logic ----------------
def iter_components(G):
    if G.is_directed():
        comps = list(nx.weakly_connected_components(G))
    else:
        comps = list(nx.connected_components(G))
    comps.sort(key=len, reverse=True)
    return comps


# ---------------- pruning ----------------
def prune_total_degree_edges(G, max_k, keep_nylon_edges=False, nylon_attr="nylon_hit"):
    if max_k is None or max_k <= 0:
        return G

    H = G.copy()

    def fnum(x):
        try:
            return float(x)
        except Exception:
            return 0.0

    def edge_rank(u, v, d):
        abs_w = fnum(d.get("abs_weight", d.get("weight_raw", 0)))
        ng = fnum(d.get("num_genomes", d.get("genome_support", 0)))
        norm = fnum(d.get("norm_weight", d.get("weight_norm", 0)))
        return (abs_w, ng, norm, str(u), str(v))

    def is_protected_edge(u, v):
        if not keep_nylon_edges:
            return False
        return (
            str(H.nodes[u].get(nylon_attr, 0)) == "1"
            or str(H.nodes[v].get(nylon_attr, 0)) == "1"
        )

    def ek(u, v):
        return (str(u), str(v)) if H.is_directed() else tuple(sorted((str(u), str(v))))

    keep = set()
    protected_incident = {n: 0 for n in H.nodes()}

    for u, v, d in H.edges(data=True):
        if is_protected_edge(u, v):
            keep.add(ek(u, v))
            protected_incident[u] += 1
            protected_incident[v] += 1

    candidates = []
    for u, v, d in H.edges(data=True):
        if ek(u, v) in keep:
            continue
        candidates.append((u, v, d))

    candidates.sort(key=lambda e: edge_rank(e[0], e[1], e[2]), reverse=True)
    incident = dict(protected_incident)

    for u, v, d in candidates:
        if incident[u] >= max_k or incident[v] >= max_k:
            continue
        keep.add(ek(u, v))
        incident[u] += 1
        incident[v] += 1

    H2 = H.__class__()
    H2.add_nodes_from(H.nodes(data=True))
    for u, v, d in H.edges(data=True):
        if ek(u, v) in keep:
            H2.add_edge(u, v, **d)
    return H2


def drop_singletons_and_doublets(G, drop_singletons=False, drop_doublets=False):
    if not (drop_singletons or drop_doublets):
        return G

    H = G.copy()
    comps = iter_components(H)
    to_remove = set()
    for comp in comps:
        n = len(comp)
        if drop_singletons and n == 1:
            to_remove |= set(comp)
        if drop_doublets and n == 2:
            to_remove |= set(comp)

    if to_remove:
        H.remove_nodes_from(to_remove)
    return H


# ---------------- layout ----------------
def compute_component_layout(H, seed, scale, pair_scale, k_small, iter_small, layout_mode="spring"):
    n_nodes = H.number_of_nodes()
    Hpos = H.to_undirected() if H.is_directed() else H

    if n_nodes == 1:
        n = next(iter(H.nodes()))
        return {n: (0.0, 0.0)}
    if n_nodes == 2:
        n1, n2 = list(H.nodes())
        return scale_positions({n1: (0.0, 0.0), n2: (1.0, 0.0)}, pair_scale)

    if layout_mode == "kamada_kawai":
        pos = nx.kamada_kawai_layout(Hpos)
        return scale_positions(normalize_to_origin(pos), scale)

    pos = nx.spring_layout(Hpos, seed=seed, k=k_small, iterations=iter_small)
    return scale_positions(normalize_to_origin(pos), scale)


def compute_spring_positions_per_component(
    G,
    seed,
    k_small,
    iter_small,
    scale,
    padding,
    pair_scale,
    target_aspect=1.0,
):
    component_entries = []
    comps = iter_components(G)

    for comp in comps:
        H = G.subgraph(comp).copy()
        pos = compute_component_layout(
            H,
            seed=seed,
            scale=scale,
            pair_scale=pair_scale,
            k_small=k_small,
            iter_small=iter_small,
        )
        component_entries.append({"pos": pos, "nodes": list(comp)})

    packed = pack_components_by_size_bands(component_entries, G, padding=padding, target_aspect=target_aspect)
    packed = fit_to_aspect_ratio(normalize_to_origin(packed), target_aspect=target_aspect)
    return normalize_to_origin(packed)


def compute_sfdp_positions(H):
    try:
        from networkx.drawing.nx_agraph import graphviz_layout
    except Exception as e:
        raise RuntimeError(
            "Could not import pygraphviz / nx_agraph graphviz_layout. "
            "Install with: conda install -c conda-forge graphviz pygraphviz"
        ) from e

    Hpos = H.to_undirected() if H.is_directed() else H
    pos = graphviz_layout(Hpos, prog="sfdp")
    return {n: (float(x), float(y)) for n, (x, y) in pos.items()}


def build_layout_graph_for_hairball(H, degree_power=1.0):
    """
    Preserve every node and edge in the exported graph, but reduce the layout pull
    of hub-heavy edges so a few chokepoints do not collapse the whole layout.
    """
    base = H.to_undirected() if H.is_directed() else H.copy()
    L = nx.Graph()
    for n, data in base.nodes(data=True):
        L.add_node(n, **data)

    deg = dict(base.degree())

    def fnum(x):
        try:
            return float(x)
        except Exception:
            return 0.0

    for u, v, d in base.edges(data=True):
        raw = fnum(d.get("norm_weight", d.get("weight_norm", d.get("abs_weight", 1.0))))
        hub_penalty = ((deg.get(u, 1) + 1) * (deg.get(v, 1) + 1)) ** (0.5 * degree_power)
        layout_weight = max(raw, 0.05) / max(hub_penalty, 1e-9)
        L.add_edge(u, v, layout_weight=layout_weight, **d)
    return L


def compute_hairball_positions(
    G_hair,
    seed,
    hairball_scale,
    hairball_k,
    hairball_iter,
    target_aspect=1.0,
    mode="hub_relaxed_spring",
):
    if mode == "sfdp":
        hair_positions = compute_sfdp_positions(G_hair)
        hair_positions = normalize_to_origin(hair_positions)
        if hairball_scale != 1.0:
            hair_positions = scale_positions(hair_positions, hairball_scale)
        return fit_to_aspect_ratio(hair_positions, target_aspect=target_aspect)

    layout_graph = build_layout_graph_for_hairball(G_hair)
    pos = nx.spring_layout(
        layout_graph,
        seed=seed,
        k=hairball_k,
        iterations=hairball_iter,
        weight="layout_weight",
    )
    pos = normalize_to_origin(pos)
    pos = scale_positions(pos, hairball_scale)
    return normalize_to_origin(fit_to_aspect_ratio(pos, target_aspect=target_aspect))




def compact_edges_locally(G, pos, passes=25, alpha=0.08, protect_singletons=True):
    """
    Gentle post-layout compaction pass:
    - preserves the broad component arrangement
    - shortens overly long local edges
    - keeps hubs from yanking the whole graph into a point
    """
    if not pos:
        return {}
    H = G.to_undirected() if G.is_directed() else G
    new_pos = {n: (float(x), float(y)) for n, (x, y) in pos.items()}

    for _ in range(max(0, passes)):
        updated = {}
        for n in H.nodes():
            nbrs = [m for m in H.neighbors(n) if m in new_pos]
            if not nbrs:
                if protect_singletons:
                    updated[n] = new_pos[n]
                continue
            x, y = new_pos[n]
            mx = sum(new_pos[m][0] for m in nbrs) / len(nbrs)
            my = sum(new_pos[m][1] for m in nbrs) / len(nbrs)
            updated[n] = (x * (1.0 - alpha) + mx * alpha, y * (1.0 - alpha) + my * alpha)
        new_pos.update(updated)
    return new_pos

# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gml", required=True)
    ap.add_argument("--out_prefix", default="layout")

    ap.add_argument("--drop_singletons", action="store_true")
    ap.add_argument("--drop_doublets", action="store_true")

    ap.add_argument("--max_k", type=int, default=0, help="Limit total incident edges per node to top-K. 0 disables.")
    ap.add_argument("--keep_nylon_edges", action="store_true", help="If nodes have nylon_hit=1, never prune edges touching those nodes.")

    ap.add_argument("--hairball_min_nodes", type=int, default=300)
    ap.add_argument("--hairball_scale", type=float, default=1400.0)
    ap.add_argument("--hairball_layout", choices=["hub_relaxed_spring", "sfdp"], default="hub_relaxed_spring")
    ap.add_argument("--hairball_k", type=float, default=1.8)
    ap.add_argument("--hairball_iter", type=int, default=800)

    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--k_small", type=float, default=0.8)
    ap.add_argument("--iter_small", type=int, default=300)
    ap.add_argument("--scale", type=float, default=450.0)
    ap.add_argument("--padding", type=float, default=180.0)
    ap.add_argument("--pair_scale", type=float, default=150.0)
    ap.add_argument("--target_aspect", type=float, default=1.2, help="Desired width/height ratio for packed layouts.")
    ap.add_argument("--compact_passes", type=int, default=25, help="Neighbor-averaging passes to shorten long local edges after layout.")
    ap.add_argument("--compact_alpha", type=float, default=0.08, help="Strength of each compaction pass (0-1).")

    args = ap.parse_args()
    out_prefix = Path(args.out_prefix)

    G = nx.read_gml(args.gml)
    G = coerce_graph_numeric_attributes(G)

    if args.max_k and args.max_k > 0:
        G = prune_total_degree_edges(G, max_k=args.max_k, keep_nylon_edges=args.keep_nylon_edges)

    G = drop_singletons_and_doublets(
        G,
        drop_singletons=args.drop_singletons,
        drop_doublets=args.drop_doublets,
    )

    if G.number_of_nodes() == 0:
        raise SystemExit("Graph has 0 nodes after filtering; nothing to export.")

    comps = iter_components(G)
    largest = comps[0]
    largest_n = len(largest)

    if largest_n < args.hairball_min_nodes:
        print(f"No hairball detected: largest component has {largest_n} nodes (< {args.hairball_min_nodes}).")
        positions = compute_spring_positions_per_component(
            G,
            seed=args.seed,
            k_small=args.k_small,
            iter_small=args.iter_small,
            scale=args.scale,
            padding=args.padding,
            pair_scale=args.pair_scale,
            target_aspect=args.target_aspect,
        )
        positions = compact_edges_locally(G, positions, passes=args.compact_passes, alpha=args.compact_alpha)
        positions = normalize_to_origin(fit_to_aspect_ratio(positions, target_aspect=args.target_aspect))
        cyjs = to_cytoscapejs_json(G, positions)
        json_path = out_prefix.with_suffix(".cyjs.json")
        tsv_path = out_prefix.with_suffix(".positions.tsv")
        json_path.write_text(json.dumps(cyjs), encoding="utf-8")
        write_positions_tsv(tsv_path, positions)
        print(f"Wrote: {json_path} and {tsv_path}")
        print(f"Export graph has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
        return

    hair_nodes = set(largest)
    G_hair = G.subgraph(hair_nodes).copy()
    G_cruise = G.copy()
    G_cruise.remove_nodes_from(hair_nodes)

    print(f"Hairball detected: largest component has {largest_n} nodes (>= {args.hairball_min_nodes}).")
    print(f"Cruise graph: {G_cruise.number_of_nodes()} nodes, {G_cruise.number_of_edges()} edges")
    print(f"Hairball graph: {G_hair.number_of_nodes()} nodes, {G_hair.number_of_edges()} edges")

    cruise_positions = compute_spring_positions_per_component(
        G_cruise,
        seed=args.seed,
        k_small=args.k_small,
        iter_small=args.iter_small,
        scale=args.scale,
        padding=args.padding,
        pair_scale=args.pair_scale,
        target_aspect=args.target_aspect,
    )
    cruise_positions = compact_edges_locally(G_cruise, cruise_positions, passes=args.compact_passes, alpha=args.compact_alpha)
    cruise_positions = normalize_to_origin(fit_to_aspect_ratio(cruise_positions, target_aspect=args.target_aspect))
    cruise_json = to_cytoscapejs_json(G_cruise, cruise_positions)
    cruise_json_path = out_prefix.with_name(out_prefix.name + "_cruise").with_suffix(".cyjs.json")
    cruise_tsv_path = out_prefix.with_name(out_prefix.name + "_cruise").with_suffix(".positions.tsv")
    cruise_json_path.write_text(json.dumps(cruise_json), encoding="utf-8")
    write_positions_tsv(cruise_tsv_path, cruise_positions)

    hair_positions = compute_hairball_positions(
        G_hair,
        seed=args.seed,
        hairball_scale=args.hairball_scale,
        hairball_k=args.hairball_k,
        hairball_iter=args.hairball_iter,
        target_aspect=args.target_aspect,
        mode=args.hairball_layout,
    )
    hair_positions = compact_edges_locally(G_hair, hair_positions, passes=args.compact_passes, alpha=args.compact_alpha)
    hair_positions = normalize_to_origin(fit_to_aspect_ratio(hair_positions, target_aspect=args.target_aspect))
    hair_json = to_cytoscapejs_json(G_hair, hair_positions)
    hair_json_path = out_prefix.with_name(out_prefix.name + "_hairball").with_suffix(".cyjs.json")
    hair_tsv_path = out_prefix.with_name(out_prefix.name + "_hairball").with_suffix(".positions.tsv")
    hair_json_path.write_text(json.dumps(hair_json), encoding="utf-8")
    write_positions_tsv(hair_tsv_path, hair_positions)

    print(f"Wrote: {cruise_json_path} and {cruise_tsv_path}")
    print(f"Wrote: {hair_json_path} and {hair_tsv_path}")

    # Build an all-in-one ordered layout where the largest component (hairball)
    # sits at the top-left and progressively smaller components trend toward the bottom-right.
    all_component_entries = [{"pos": hair_positions, "nodes": list(G_hair.nodes())}]
    for comp in iter_components(G_cruise):
        subpos = {n: cruise_positions[n] for n in comp if n in cruise_positions}
        all_component_entries.append({"pos": subpos, "nodes": list(comp)})

    G_all = nx.compose(G_cruise, G_hair)
    all_positions = pack_components_by_size_bands(all_component_entries, G_all, padding=max(args.padding, 220.0), target_aspect=args.target_aspect)
    all_positions = normalize_to_origin(fit_to_aspect_ratio(all_positions, target_aspect=args.target_aspect))
    all_json = to_cytoscapejs_json(G_all, all_positions)
    all_json_path = out_prefix.with_name(out_prefix.name + "_allinone").with_suffix(".cyjs.json")
    all_tsv_path = out_prefix.with_name(out_prefix.name + "_allinone").with_suffix(".positions.tsv")
    all_json_path.write_text(json.dumps(all_json), encoding="utf-8")
    write_positions_tsv(all_tsv_path, all_positions)

    print(f"Wrote: {all_json_path} and {all_tsv_path}")


if __name__ == "__main__":
    main()
