import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path
from statistics import mean, median, pstdev
from typing import Dict, Iterable, List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx


COMPONENT_SIZE_BINS = [
    ("1", 1, 1),
    ("2", 2, 2),
    ("3-5", 3, 5),
    ("6-10", 6, 10),
    ("11-20", 11, 20),
    ("21-50", 21, 50),
    ("51-100", 51, 100),
    (">100", 101, math.inf),
]


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def coerce_int(value, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def percentile(values: Sequence[float], q: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(float(v) for v in values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * q
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def basic_distribution_stats(values: Sequence[float]) -> Dict[str, float]:
    if not values:
        return {"mean": 0.0, "median": 0.0, "std": 0.0, "max": 0.0}
    return {
        "mean": float(mean(values)),
        "median": float(median(values)),
        "std": float(pstdev(values)) if len(values) > 1 else 0.0,
        "max": float(max(values)),
    }


def discrete_histogram(values: Sequence[int]) -> Dict[str, int]:
    counts = Counter(int(v) for v in values)
    return {str(key): counts[key] for key in sorted(counts)}


def choose_node_label(node_id: str, attrs: Dict) -> str:
    for key in ("reference_gene_name", "gene_name", "reference_product", "product"):
        value = str(attrs.get(key, "") or "").strip()
        if not value:
            continue
        if value.lower() in {"unknown gene", "unknown product"}:
            continue
        return value
    return str(node_id)


def format_run_label(summary: Dict) -> str:
    params = summary.get("parameters", {})
    if not params:
        return summary.get("run_name", "unknown_run")
    return (
        f"one={params.get('enforce_one_per_genome')} | "
        f"op={params.get('operonic_distance')} | "
        f"s1={params.get('stage1_distance')} | "
        f"s2={params.get('stage2_distance')} | "
        f"s3={params.get('stage3_distance')}"
    )


def component_bin_percentages(component_sizes: Sequence[int], total_nodes: int) -> Dict[str, float]:
    counts = {label: 0 for label, _, _ in COMPONENT_SIZE_BINS}
    if total_nodes <= 0:
        return counts

    for size in component_sizes:
        for label, lower, upper in COMPONENT_SIZE_BINS:
            if lower <= size <= upper:
                counts[label] += size
                break

    return {label: (count / total_nodes) * 100.0 for label, count in counts.items()}


def summarize_graph(graph: nx.Graph, run_name: str, parameters: Optional[Dict] = None, top_hubs: int = 15):
    undirected = graph.to_undirected()

    component_sizes = sorted((len(component) for component in nx.connected_components(undirected)), reverse=True)
    total_nodes = undirected.number_of_nodes()
    total_edges = undirected.number_of_edges()
    cluster_sizes = [coerce_int(graph.nodes[node].get("cluster_size"), 1) for node in graph.nodes()]
    total_input_proteins = int(sum(cluster_sizes))
    compression_ratio = float(total_input_proteins / total_nodes) if total_nodes else 0.0

    degrees = [int(degree) for _, degree in undirected.degree()]
    largest_component_fraction = (
        float(component_sizes[0] / total_nodes) if component_sizes and total_nodes else 0.0
    )
    largest_component_size = int(component_sizes[0]) if component_sizes else 0
    number_of_components = int(len(component_sizes))
    number_of_singleton_nodes = int(sum(size for size in component_sizes if size == 1))
    singleton_fraction = (
        float(number_of_singleton_nodes / total_nodes) if total_nodes else 0.0
    )
    number_of_multi_member_clusters = int(sum(1 for size in cluster_sizes if size > 1))

    hub_rows = []
    for node, degree in sorted(undirected.degree(), key=lambda item: (-item[1], str(item[0]))):
        attrs = graph.nodes[node]
        hub_rows.append(
            {
                "node_id": str(node),
                "label": choose_node_label(str(node), attrs),
                "degree": int(degree),
                "cluster_size": coerce_int(attrs.get("cluster_size"), 1),
                "product": str(attrs.get("product", "") or ""),
                "gene_name": str(attrs.get("gene_name", "") or ""),
            }
        )

    summary = {
        "run_name": run_name,
        "parameters": parameters or {},
        "graph_metrics": {
            "total_input_proteins": total_input_proteins,
            "total_nodes": total_nodes,
            "total_edges": total_edges,
            "compression_ratio": compression_ratio,
            "graph_density": float(nx.density(undirected)) if total_nodes > 1 else 0.0,
            "clustering_coefficient": (
                float(nx.average_clustering(undirected)) if total_nodes > 1 else 0.0
            ),
        },
        "connected_components": {
            "total_components": number_of_components,
            "largest_component_size": largest_component_size,
            "largest_component_fraction": largest_component_fraction,
            "number_of_singleton_nodes": number_of_singleton_nodes,
            "singleton_node_fraction": singleton_fraction,
            "size_histogram": discrete_histogram(component_sizes),
            "node_fraction_percent_by_size_bin": component_bin_percentages(component_sizes, total_nodes),
            "values": component_sizes,
        },
        "proteins_per_node": {
            **basic_distribution_stats(cluster_sizes),
            "histogram": discrete_histogram(cluster_sizes),
            "number_of_multi_member_clusters": number_of_multi_member_clusters,
            "values": cluster_sizes,
        },
        "node_degree": {
            **basic_distribution_stats(degrees),
            "p95": float(percentile(degrees, 0.95)),
            "p99": float(percentile(degrees, 0.99)),
            "histogram": discrete_histogram(degrees),
            "values": degrees,
        },
        "top_hubs": hub_rows[:top_hubs],
    }
    return summary


def save_figure(fig: plt.Figure, out_base: Path):
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=200, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_run_distributions(summary: Dict, plots_dir: Path):
    plots_dir = ensure_dir(plots_dir)

    component_hist = summary["connected_components"]["size_histogram"]
    component_bin_pct = summary["connected_components"]["node_fraction_percent_by_size_bin"]
    cluster_hist = summary["proteins_per_node"]["histogram"]
    degree_hist = summary["node_degree"]["histogram"]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Operon Graph Summary: {summary['run_name']}", fontsize=14)

    axes[0, 0].bar(list(component_hist.keys()), list(component_hist.values()), color="#4C78A8")
    axes[0, 0].set_title("Connected Component Sizes")
    axes[0, 0].set_xlabel("Component Size")
    axes[0, 0].set_ylabel("Number of Components")

    axes[0, 1].bar(list(component_bin_pct.keys()), list(component_bin_pct.values()), color="#F58518")
    axes[0, 1].set_title("Nodes by Component Size Bin")
    axes[0, 1].set_xlabel("Component Size Bin")
    axes[0, 1].set_ylabel("Percent of Nodes")
    axes[0, 1].tick_params(axis="x", rotation=30)

    axes[1, 0].bar(list(cluster_hist.keys()), list(cluster_hist.values()), color="#54A24B")
    axes[1, 0].set_title("Proteins per Node")
    axes[1, 0].set_xlabel("Proteins in Cluster")
    axes[1, 0].set_ylabel("Node Count")

    axes[1, 1].bar(list(degree_hist.keys()), list(degree_hist.values()), color="#E45756")
    axes[1, 1].set_title("Node Degree Distribution")
    axes[1, 1].set_xlabel("Undirected Degree")
    axes[1, 1].set_ylabel("Node Count")

    save_figure(fig, plots_dir / "distributions")

    top_hubs = summary.get("top_hubs", [])[:15]
    if top_hubs:
        fig, ax = plt.subplots(figsize=(14, 6))
        labels = [hub["label"] for hub in top_hubs]
        degrees = [hub["degree"] for hub in top_hubs]
        ax.bar(range(len(top_hubs)), degrees, color="#72B7B2")
        ax.set_title(f"Top Hub Nodes: {summary['run_name']}")
        ax.set_xlabel("Node Annotation")
        ax.set_ylabel("Undirected Degree")
        ax.set_xticks(range(len(top_hubs)))
        ax.set_xticklabels(labels, rotation=50, ha="right")
        save_figure(fig, plots_dir / "top_hubs")


def write_summary_json(summary: Dict, summary_json: Path):
    ensure_dir(summary_json.parent)
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True))


def _histogram_total_count(histogram: Dict[str, int], key_filter=None) -> int:
    total = 0
    for key, count in histogram.items():
        try:
            numeric_key = int(key)
        except (TypeError, ValueError):
            continue
        if key_filter is None or key_filter(numeric_key):
            total += int(count)
    return total


def summary_to_master_row(summary: Dict) -> Dict[str, object]:
    params = summary.get("parameters", {})
    graph_metrics = summary["graph_metrics"]
    components = summary["connected_components"]
    proteins = summary["proteins_per_node"]
    degree = summary["node_degree"]
    top_hub = summary.get("top_hubs", [{}])[0] if summary.get("top_hubs") else {}
    component_values = components.get("values", [])
    component_histogram = components.get("size_histogram", {})
    protein_histogram = proteins.get("histogram", {})

    largest_component_size = components.get("largest_component_size")
    if largest_component_size is None:
        if component_values:
            largest_component_size = max(component_values)
        elif component_histogram:
            largest_component_size = max(int(key) for key in component_histogram.keys())
        else:
            largest_component_size = 0

    number_of_components = components.get("total_components")
    if number_of_components is None:
        if component_values:
            number_of_components = len(component_values)
        else:
            number_of_components = sum(int(count) for count in component_histogram.values())

    number_of_singleton_nodes = components.get("number_of_singleton_nodes")
    if number_of_singleton_nodes is None:
        if component_values:
            number_of_singleton_nodes = sum(1 for value in component_values if int(value) == 1)
        else:
            number_of_singleton_nodes = int(component_histogram.get("1", 0))

    number_of_multi_member_clusters = proteins.get("number_of_multi_member_clusters")
    if number_of_multi_member_clusters is None:
        if proteins.get("values"):
            number_of_multi_member_clusters = sum(1 for value in proteins["values"] if float(value) > 1)
        else:
            number_of_multi_member_clusters = _histogram_total_count(
                protein_histogram,
                key_filter=lambda value: value > 1,
            )

    return {
        "run_name": summary["run_name"],
        "enforce_one_per_genome": params.get("enforce_one_per_genome"),
        "operonic_distance": params.get("operonic_distance"),
        "stage1_distance": params.get("stage1_distance"),
        "stage2_distance": params.get("stage2_distance"),
        "stage3_distance": params.get("stage3_distance"),
        "total_input_proteins": graph_metrics["total_input_proteins"],
        "total_nodes": graph_metrics["total_nodes"],
        "total_edges": graph_metrics["total_edges"],
        "compression_ratio": graph_metrics["compression_ratio"],
        "graph_density": graph_metrics["graph_density"],
        "clustering_coefficient": graph_metrics["clustering_coefficient"],
        "largest_component_size": largest_component_size,
        "number_of_components": number_of_components,
        "number_of_singleton_nodes": number_of_singleton_nodes,
        "largest_component_fraction": components["largest_component_fraction"],
        "singleton_node_fraction": components["singleton_node_fraction"],
        "proteins_per_node_mean": proteins["mean"],
        "proteins_per_node_median": proteins["median"],
        "proteins_per_node_std": proteins["std"],
        "proteins_per_node_max": proteins["max"],
        "number_of_multi_member_clusters": number_of_multi_member_clusters,
        "degree_mean": degree["mean"],
        "degree_median": degree["median"],
        "degree_max": degree["max"],
        "degree_p95": degree["p95"],
        "degree_p99": degree["p99"],
        "top_hub_label": top_hub.get("label", ""),
        "top_hub_degree": top_hub.get("degree", 0),
    }


def write_master_summary_csv(summaries: Sequence[Dict], out_csv: Path):
    rows = [summary_to_master_row(summary) for summary in summaries]
    if not rows:
        return
    ensure_dir(out_csv.parent)
    fieldnames = list(rows[0].keys())
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def value_to_label(value) -> str:
    value = float(value)
    if value.is_integer():
        return str(int(value))
    return str(value)


def metric_value(summary: Dict, metric_key: str):
    mapping = {
        "total_nodes": summary["graph_metrics"]["total_nodes"],
        "total_edges": summary["graph_metrics"]["total_edges"],
        "fraction_in_largest_component": summary["connected_components"]["largest_component_fraction"],
        "singleton_node_fraction": summary["connected_components"]["singleton_node_fraction"],
        "mean_degree": summary["node_degree"]["mean"],
        "max_degree": summary["node_degree"]["max"],
        "proteins_per_node_mean": summary["proteins_per_node"]["mean"],
        "proteins_per_node_median": summary["proteins_per_node"]["median"],
        "proteins_per_node_max": summary["proteins_per_node"]["max"],
        "compression_ratio": summary["graph_metrics"]["compression_ratio"],
        "largest_component_size": summary["connected_components"]["largest_component_size"],
        "number_of_connected_components": summary["connected_components"]["total_components"],
    }
    return float(mapping[metric_key])


def sort_summaries_for_runs(summaries: Sequence[Dict]) -> List[Dict]:
    return sorted(
        summaries,
        key=lambda summary: (
            str(summary.get("parameters", {}).get("enforce_one_per_genome")),
            float(summary.get("parameters", {}).get("operonic_distance", 0)),
            float(summary.get("parameters", {}).get("stage1_distance", 0)),
            float(summary.get("parameters", {}).get("stage2_distance", 0)),
            float(summary.get("parameters", {}).get("stage3_distance", 0)),
        ),
    )


def filter_summaries_by_mode(summaries: Sequence[Dict], enforce_one_per_genome: bool) -> List[Dict]:
    return [
        summary
        for summary in summaries
        if bool(summary.get("parameters", {}).get("enforce_one_per_genome")) == enforce_one_per_genome
    ]


def filter_summaries_by_operonic_distance(summaries: Sequence[Dict], operonic_distance: float) -> List[Dict]:
    return [
        summary
        for summary in summaries
        if float(summary.get("parameters", {}).get("operonic_distance", 0)) == operonic_distance
    ]


def filter_summaries_by_stage1_distance(summaries: Sequence[Dict], stage1_distance: float) -> List[Dict]:
    return [
        summary
        for summary in summaries
        if float(summary.get("parameters", {}).get("stage1_distance", 0)) == stage1_distance
    ]


def build_heatmap_matrix(summaries: Sequence[Dict], metric_key: str):
    stage2_values = sorted(
        {float(summary.get("parameters", {}).get("stage2_distance", 0)) for summary in summaries}
    )
    stage3_values = sorted(
        {float(summary.get("parameters", {}).get("stage3_distance", 0)) for summary in summaries}
    )
    matrix = []
    for stage2 in stage2_values:
        row = []
        for stage3 in stage3_values:
            matched = None
            for summary in summaries:
                params = summary.get("parameters", {})
                if (
                    float(params.get("stage2_distance", 0)) == stage2
                    and float(params.get("stage3_distance", 0)) == stage3
                ):
                    matched = summary
                    break
            row.append(metric_value(matched, metric_key) if matched is not None else float("nan"))
        matrix.append(row)
    return stage2_values, stage3_values, matrix


def annotate_heatmap(ax, matrix):
    for row_index, row in enumerate(matrix):
        for col_index, value in enumerate(row):
            if math.isnan(value):
                label = "NA"
            elif abs(value) >= 100 or float(value).is_integer():
                label = f"{value:.0f}"
            else:
                label = f"{value:.3g}"
            ax.text(col_index, row_index, label, ha="center", va="center", fontsize=8, color="black")


def plot_metric_heatmap(
    summaries: Sequence[Dict],
    metric_key: str,
    metric_title: str,
    output_dir: Path,
    enforce_one_per_genome: bool,
    operonic_distance: float,
    stage1_distance: float,
    delta_from_baseline: bool = False,
):
    if not summaries:
        return

    stage2_values, stage3_values, matrix = build_heatmap_matrix(summaries, metric_key)
    if not matrix:
        return

    plotted_matrix = [list(row) for row in matrix]
    suffix = "delta_" if delta_from_baseline else ""
    if delta_from_baseline:
        baseline = plotted_matrix[0][0]
        plotted_matrix = [
            [
                (value - baseline) if not math.isnan(value) and not math.isnan(baseline) else float("nan")
                for value in row
            ]
            for row in plotted_matrix
        ]

    fig, ax = plt.subplots(figsize=(6, 5))
    cmap = "coolwarm" if delta_from_baseline else "viridis"
    image = ax.imshow(plotted_matrix, aspect="auto", cmap=cmap)
    ax.set_title(
        f"{metric_title} | one_per_genome={str(enforce_one_per_genome).lower()} | "
        f"op={value_to_label(operonic_distance)} | s1={value_to_label(stage1_distance)}"
        + (" | delta vs (1000,1000)" if delta_from_baseline else "")
    )
    ax.set_xlabel("stage3_distance")
    ax.set_ylabel("stage2_distance")
    ax.set_xticks(range(len(stage3_values)))
    ax.set_xticklabels([value_to_label(value) for value in stage3_values])
    ax.set_yticks(range(len(stage2_values)))
    ax.set_yticklabels([value_to_label(value) for value in stage2_values])
    annotate_heatmap(ax, plotted_matrix)
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label(metric_title)
    save_figure(
        fig,
        output_dir
        / (
            f"heatmap_{suffix}{metric_key}_onepergenome_{str(enforce_one_per_genome).lower()}_"
            f"operondist_{value_to_label(operonic_distance)}_"
            f"stage1_{value_to_label(stage1_distance)}"
        ),
    )


def plot_all_heatmaps(summaries: Sequence[Dict], output_dir: Path):
    heatmap_dir = ensure_dir(output_dir / "heatmaps")
    delta_dir = ensure_dir(output_dir / "delta_heatmaps")
    metric_specs = [
        ("total_nodes", "Total Nodes"),
        ("total_edges", "Total Edges"),
        ("fraction_in_largest_component", "Fraction in Largest Component"),
        ("singleton_node_fraction", "Singleton Node Fraction"),
        ("mean_degree", "Mean Degree"),
        ("max_degree", "Max Degree"),
        ("proteins_per_node_mean", "Proteins per Node Mean"),
        ("proteins_per_node_median", "Proteins per Node Median"),
        ("proteins_per_node_max", "Proteins per Node Max"),
        ("compression_ratio", "Compression Ratio"),
        ("largest_component_size", "Largest Component Size"),
        ("number_of_connected_components", "Number of Connected Components"),
    ]

    operonic_distances = sorted(
        {float(summary.get("parameters", {}).get("operonic_distance", 0)) for summary in summaries}
    )
    for operonic_distance in operonic_distances:
        op_summaries = filter_summaries_by_operonic_distance(summaries, operonic_distance)
        stage1_distances = sorted(
            {float(summary.get("parameters", {}).get("stage1_distance", 0)) for summary in op_summaries}
        )
        for stage1_distance in stage1_distances:
            s1_summaries = filter_summaries_by_stage1_distance(op_summaries, stage1_distance)
            for enforce_one_per_genome in (True, False):
                mode_summaries = filter_summaries_by_mode(s1_summaries, enforce_one_per_genome)
                if not mode_summaries:
                    continue
                for metric_key, metric_title in metric_specs:
                    plot_metric_heatmap(
                        mode_summaries,
                        metric_key=metric_key,
                        metric_title=metric_title,
                        output_dir=heatmap_dir,
                        enforce_one_per_genome=enforce_one_per_genome,
                        operonic_distance=operonic_distance,
                        stage1_distance=stage1_distance,
                        delta_from_baseline=False,
                    )
                    plot_metric_heatmap(
                        mode_summaries,
                        metric_key=metric_key,
                        metric_title=metric_title,
                        output_dir=delta_dir,
                        enforce_one_per_genome=enforce_one_per_genome,
                        operonic_distance=operonic_distance,
                        stage1_distance=stage1_distance,
                        delta_from_baseline=True,
                    )


def ecdf_points(values: Sequence[float]):
    if not values:
        return [], []
    ordered = sorted(float(value) for value in values)
    y_values = [(index + 1) / len(ordered) for index in range(len(ordered))]
    return ordered, y_values


def stage3_color_map(summaries: Sequence[Dict]) -> Dict[float, tuple]:
    stage3_values = sorted(
        {float(summary.get("parameters", {}).get("stage3_distance", 0)) for summary in summaries}
    )
    cmap = plt.get_cmap("viridis", max(len(stage3_values), 1))
    return {stage3: cmap(index) for index, stage3 in enumerate(stage3_values)}


def stage2_style_map(summaries: Sequence[Dict]) -> Dict[float, str]:
    stage2_values = sorted(
        {float(summary.get("parameters", {}).get("stage2_distance", 0)) for summary in summaries}
    )
    style_cycle = ["solid", "dashed", "dotted", "dashdot"]
    return {
        stage2: style_cycle[index % len(style_cycle)]
        for index, stage2 in enumerate(stage2_values)
    }


def run_color(summary: Dict, color_map: Dict[float, tuple]):
    stage3 = float(summary.get("parameters", {}).get("stage3_distance", 0))
    return color_map.get(stage3, "#72B7B2")


def run_line_style(summary: Dict, style_map: Dict[float, str]):
    stage2 = float(summary.get("parameters", {}).get("stage2_distance", 0))
    return style_map.get(stage2, "solid")


def plot_distribution_ecdf(
    summaries: Sequence[Dict],
    distribution_key: str,
    distribution_label: str,
    x_label: str,
    output_dir: Path,
):
    output_dir = ensure_dir(output_dir)
    color_map = stage3_color_map(summaries)
    style_map = stage2_style_map(summaries)
    operonic_distances = sorted(
        {float(summary.get("parameters", {}).get("operonic_distance", 0)) for summary in summaries}
    )
    for operonic_distance in operonic_distances:
        op_summaries = filter_summaries_by_operonic_distance(summaries, operonic_distance)
        stage1_distances = sorted(
            {float(summary.get("parameters", {}).get("stage1_distance", 0)) for summary in op_summaries}
        )
        for stage1_distance in stage1_distances:
            s1_summaries = filter_summaries_by_stage1_distance(op_summaries, stage1_distance)
            for enforce_one_per_genome in (True, False):
                mode_summaries = sort_summaries_for_runs(filter_summaries_by_mode(s1_summaries, enforce_one_per_genome))
                if not mode_summaries:
                    continue
                fig, ax = plt.subplots(figsize=(10, 6))
                for summary in mode_summaries:
                    values = summary[distribution_key]["values"] if distribution_key != "connected_components" else summary["connected_components"]["values"]
                    x_values, y_values = ecdf_points(values)
                    ax.plot(
                        x_values,
                        y_values,
                        color=run_color(summary, color_map),
                        linestyle=run_line_style(summary, style_map),
                        linewidth=1.8,
                        label=format_run_label(summary),
                    )
                ax.set_title(
                    f"{distribution_label} ECDF | one_per_genome={str(enforce_one_per_genome).lower()} | "
                    f"op={value_to_label(operonic_distance)} | s1={value_to_label(stage1_distance)}"
                )
                ax.set_xlabel(x_label)
                ax.set_ylabel("ECDF")
                ax.legend(fontsize=8, ncol=1)
                ax.grid(True, axis="both", alpha=0.25)
                filename = (
                    f"ecdf_{distribution_key}_onepergenome_{str(enforce_one_per_genome).lower()}_"
                    f"operondist_{value_to_label(operonic_distance)}_stage1_{value_to_label(stage1_distance)}"
                    .replace("connected_components", "component_sizes")
                )
                save_figure(fig, output_dir / filename)


def plot_distribution_histograms(
    summaries: Sequence[Dict],
    distribution_key: str,
    distribution_label: str,
    x_label: str,
    output_dir: Path,
):
    output_dir = ensure_dir(output_dir)
    color_map = stage3_color_map(summaries)
    style_map = stage2_style_map(summaries)
    operonic_distances = sorted(
        {float(summary.get("parameters", {}).get("operonic_distance", 0)) for summary in summaries}
    )
    for operonic_distance in operonic_distances:
        op_summaries = filter_summaries_by_operonic_distance(summaries, operonic_distance)
        stage1_distances = sorted(
            {float(summary.get("parameters", {}).get("stage1_distance", 0)) for summary in op_summaries}
        )
        for stage1_distance in stage1_distances:
            s1_summaries = filter_summaries_by_stage1_distance(op_summaries, stage1_distance)
            for enforce_one_per_genome in (True, False):
                mode_summaries = sort_summaries_for_runs(filter_summaries_by_mode(s1_summaries, enforce_one_per_genome))
                if not mode_summaries:
                    continue
                stage2_groups = sorted(
                    {float(summary.get("parameters", {}).get("stage2_distance", 0)) for summary in mode_summaries}
                )
                fig, axes = plt.subplots(len(stage2_groups), 1, figsize=(10, max(4 * len(stage2_groups), 4)), sharex=False)
                if len(stage2_groups) == 1:
                    axes = [axes]
                for axis, stage2 in zip(axes, stage2_groups):
                    stage2_summaries = [
                        summary for summary in mode_summaries if float(summary.get("parameters", {}).get("stage2_distance", 0)) == stage2
                    ]
                    all_values = []
                    for summary in stage2_summaries:
                        values = summary[distribution_key]["values"] if distribution_key != "connected_components" else summary["connected_components"]["values"]
                        all_values.extend(values)
                    bins = "auto"
                    if all_values:
                        unique_count = len(set(float(value) for value in all_values))
                        bins = min(max(unique_count, 10), 40)
                    for summary in stage2_summaries:
                        values = summary[distribution_key]["values"] if distribution_key != "connected_components" else summary["connected_components"]["values"]
                        axis.hist(
                            values,
                            bins=bins,
                            density=True,
                            histtype="bar",
                            alpha=0.3,
                            color=run_color(summary, color_map),
                            edgecolor=run_color(summary, color_map),
                            linewidth=1.0,
                            label=f"s3={value_to_label(summary.get('parameters', {}).get('stage3_distance', 0))}",
                        )
                    axis.set_title(f"stage2={value_to_label(stage2)}")
                    axis.set_ylabel("Density")
                    axis.grid(True, axis="y", alpha=0.25)
                    handles, labels = axis.get_legend_handles_labels()
                    if handles:
                        dedup = dict(zip(labels, handles))
                        axis.legend(dedup.values(), dedup.keys(), fontsize=8)
                axes[-1].set_xlabel(x_label)
                fig.suptitle(
                    f"{distribution_label} Histograms | one_per_genome={str(enforce_one_per_genome).lower()} | "
                    f"op={value_to_label(operonic_distance)} | s1={value_to_label(stage1_distance)}"
                )
                filename = (
                    f"histogram_{distribution_key}_onepergenome_{str(enforce_one_per_genome).lower()}_"
                    f"operondist_{value_to_label(operonic_distance)}_stage1_{value_to_label(stage1_distance)}"
                    .replace("connected_components", "component_sizes")
                )
                save_figure(fig, output_dir / filename)


def plot_all_distribution_comparisons(summaries: Sequence[Dict], output_dir: Path):
    dist_root = ensure_dir(output_dir / "distributions")
    ecdf_dir = ensure_dir(dist_root / "ecdf")
    hist_dir = ensure_dir(dist_root / "histograms")
    distribution_specs = [
        ("proteins_per_node", "Proteins per Node", "Proteins per Node"),
        ("node_degree", "Node Degree", "Undirected Degree"),
        ("connected_components", "Connected Component Sizes", "Component Size"),
    ]
    for distribution_key, distribution_label, x_label in distribution_specs:
        plot_distribution_ecdf(
            summaries,
            distribution_key=distribution_key,
            distribution_label=distribution_label,
            x_label=x_label,
            output_dir=ecdf_dir,
        )
        plot_distribution_histograms(
            summaries,
            distribution_key=distribution_key,
            distribution_label=distribution_label,
            x_label=x_label,
            output_dir=hist_dir,
        )


def plot_hub_heatmap(summaries: Sequence[Dict], output_dir: Path, max_labels: int = 20):
    if not summaries:
        return

    output_dir = ensure_dir(output_dir)
    annotation_max_degree = {}
    for summary in summaries:
        for hub in summary.get("top_hubs", []):
            label = hub["label"]
            annotation_max_degree[label] = max(annotation_max_degree.get(label, 0), hub["degree"])

    selected_labels = [
        label for label, _ in sorted(annotation_max_degree.items(), key=lambda item: (-item[1], item[0]))[:max_labels]
    ]
    if not selected_labels:
        return

    ordered = sorted(
        summaries,
        key=lambda summary: (
            str(summary.get("parameters", {}).get("enforce_one_per_genome")),
            float(summary.get("parameters", {}).get("operonic_distance", 0)),
            float(summary.get("parameters", {}).get("stage1_distance", 0)),
            float(summary.get("parameters", {}).get("stage2_distance", 0)),
            float(summary.get("parameters", {}).get("stage3_distance", 0)),
        ),
    )
    run_labels = [format_run_label(summary) for summary in ordered]
    matrix = []
    for label in selected_labels:
        row = []
        for summary in ordered:
            degree = 0
            for hub in summary.get("top_hubs", []):
                if hub["label"] == label:
                    degree = max(degree, hub["degree"])
            row.append(degree)
        matrix.append(row)

    fig, ax = plt.subplots(figsize=(max(12, len(run_labels) * 0.7), max(8, len(selected_labels) * 0.4)))
    image = ax.imshow(matrix, aspect="auto", cmap="YlOrRd")
    ax.set_title("High-degree Node Comparison Across Runs")
    ax.set_xlabel("Parameter Setting")
    ax.set_ylabel("Node Annotation")
    ax.set_xticks(range(len(run_labels)))
    ax.set_xticklabels(run_labels, rotation=55, ha="right")
    ax.set_yticks(range(len(selected_labels)))
    ax.set_yticklabels(selected_labels)
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("Undirected Degree")
    save_figure(fig, output_dir / "hub_degree_heatmap")


def load_summary_json(summary_json: Path) -> Dict:
    return json.loads(summary_json.read_text())


def summarize_graph_file(
    graph_gml: Path,
    summary_json: Path,
    plots_dir: Path,
    run_name: Optional[str] = None,
    parameters: Optional[Dict] = None,
    top_hubs: int = 15,
) -> Dict:
    graph = nx.read_gml(graph_gml)
    summary = summarize_graph(
        graph=graph,
        run_name=run_name or graph_gml.stem,
        parameters=parameters,
        top_hubs=top_hubs,
    )
    plot_run_distributions(summary, plots_dir)
    write_summary_json(summary, summary_json)
    return summary


def collect_summaries_from_sweep_dir(sweep_dir: Path) -> List[Dict]:
    summaries = []
    for summary_json in sorted(sweep_dir.glob("run_*/summary.json")):
        summaries.append(load_summary_json(summary_json))
    return summaries


def main():
    parser = argparse.ArgumentParser(
        description="Summarize one operon graph or aggregate a parameter sweep of graph summaries."
    )
    parser.add_argument("--graph_gml", default="", help="Single run graph.gml to summarize.")
    parser.add_argument("--summary_json", default="", help="Where to write the single-run summary JSON.")
    parser.add_argument("--plots_dir", default="", help="Directory for single-run plots.")
    parser.add_argument("--run_name", default="", help="Optional human-readable run name.")
    parser.add_argument(
        "--sweep_dir",
        default="",
        help="Sweep output directory containing run_*/summary.json files for cross-run aggregation.",
    )
    parser.add_argument(
        "--master_summary_csv",
        default="",
        help="Optional output CSV for aggregate scalar metrics. Defaults inside --sweep_dir.",
    )
    parser.add_argument("--top_hubs", type=int, default=15, help="Top hub nodes to keep per run.")
    args = parser.parse_args()

    if bool(args.graph_gml) == bool(args.sweep_dir):
        raise SystemExit("Set exactly one of --graph_gml or --sweep_dir.")

    if args.graph_gml:
        graph_gml = Path(args.graph_gml)
        summary_json = Path(args.summary_json) if args.summary_json else graph_gml.with_name("summary.json")
        plots_dir = Path(args.plots_dir) if args.plots_dir else graph_gml.parent / "plots"
        summarize_graph_file(
            graph_gml=graph_gml,
            summary_json=summary_json,
            plots_dir=plots_dir,
            run_name=args.run_name or graph_gml.parent.name,
            parameters=None,
            top_hubs=args.top_hubs,
        )
        print(f"[wrote] {summary_json}")
        print(f"[wrote] {plots_dir}")
        return

    sweep_dir = Path(args.sweep_dir)
    summaries = collect_summaries_from_sweep_dir(sweep_dir)
    if not summaries:
        raise SystemExit(f"No run summaries found under {sweep_dir}")

    master_summary_csv = (
        Path(args.master_summary_csv) if args.master_summary_csv else sweep_dir / "master_summary.csv"
    )
    write_master_summary_csv(summaries, master_summary_csv)
    comparison_dir = sweep_dir / "comparison_plots"
    plot_all_heatmaps(summaries, comparison_dir)
    plot_all_distribution_comparisons(summaries, comparison_dir)
    plot_hub_heatmap(summaries, comparison_dir / "heatmaps")
    print(f"[wrote] {master_summary_csv}")
    print(f"[wrote] {comparison_dir}")


if __name__ == "__main__":
    main()
