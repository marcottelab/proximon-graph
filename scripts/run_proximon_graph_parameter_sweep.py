import argparse
import contextlib
import csv
import itertools
import json
import traceback
from pathlib import Path
from typing import Dict, List

from graph_builder_orthologs_compatibility_scaled import GraphBuilderProst
from summarize_operon_graph_sweep import (
    plot_all_distribution_comparisons,
    plot_all_heatmaps,
    plot_hub_heatmap,
    summarize_graph_file,
    summary_to_master_row,
)


DEFAULT_ONE_PER_GENOME_VALUES = [True, False]
DEFAULT_STAGE1_DISTANCE = 1000
DEFAULT_STAGE2_DISTANCES = [1000, 3000, 5000]
DEFAULT_STAGE3_DISTANCES = [1000, 3000, 5000]
DEFAULT_OPERONIC_DISTANCES = [40]


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def parse_int_list(raw_value: str) -> List[int]:
    return [int(part.strip()) for part in raw_value.split(",") if part.strip()]


def parse_bool_list(raw_value: str) -> List[bool]:
    values = []
    for part in raw_value.split(","):
        token = part.strip().lower()
        if not token:
            continue
        if token in {"true", "t", "1", "yes", "y"}:
            values.append(True)
        elif token in {"false", "f", "0", "no", "n"}:
            values.append(False)
        else:
            raise ValueError(f"Unrecognized boolean value: {part}")
    return values


def build_run_name(
    operonic_distance: int,
    stage1_distance: int,
    stage2_distance: int,
    stage3_distance: int,
    enforce_one_per_genome: bool,
) -> str:
    return (
        f"operondist_{operonic_distance}_stage1_{stage1_distance}_stage2_{stage2_distance}_"
        f"stage3_{stage3_distance}_onepergenome_{str(enforce_one_per_genome).lower()}"
    )


def enumerate_parameter_grid(
    operonic_distances: List[int],
    stage1_distances: List[int],
    one_per_genome_values: List[bool],
    stage2_distances: List[int],
    stage3_distances: List[int],
) -> List[Dict]:
    combinations = []
    for operonic_distance, stage1_distance, enforce_one_per_genome, stage2_distance, stage3_distance in itertools.product(
        operonic_distances,
        stage1_distances,
        one_per_genome_values,
        stage2_distances,
        stage3_distances,
    ):
        combinations.append(
            {
                "operonic_distance": operonic_distance,
                "stage1_distance": stage1_distance,
                "stage2_distance": stage2_distance,
                "stage3_distance": stage3_distance,
                "enforce_one_per_genome": enforce_one_per_genome,
                "run_name": build_run_name(
                    operonic_distance=operonic_distance,
                    stage1_distance=stage1_distance,
                    stage2_distance=stage2_distance,
                    stage3_distance=stage3_distance,
                    enforce_one_per_genome=enforce_one_per_genome,
                ),
            }
        )
    return combinations


def write_master_rows(rows: List[Dict], out_csv: Path):
    if not rows:
        return
    ensure_dir(out_csv.parent)
    fieldnames = []
    seen = set()
    for row in rows:
        for key in row.keys():
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_single_configuration(args, params: Dict, sweep_root: Path) -> Dict:
    run_dir = ensure_dir(sweep_root / f"run_{params['run_name']}")
    log_path = run_dir / "log.txt"
    graph_path = run_dir / "graph.gml"
    summary_json = run_dir / "summary.json"
    plots_dir = run_dir / "plots"

    result = {
        "run_name": params["run_name"],
        "run_dir": str(run_dir),
        "status": "pending",
        "error": "",
        "graph_path": str(graph_path),
        "summary_json": str(summary_json),
        "operonic_distance": params["operonic_distance"],
        "stage1_distance": params["stage1_distance"],
        "stage2_distance": params["stage2_distance"],
        "stage3_distance": params["stage3_distance"],
        "enforce_one_per_genome": params["enforce_one_per_genome"],
    }

    if not args.force and graph_path.exists() and summary_json.exists():
        summary = json.loads(summary_json.read_text())
        result.update(summary_to_master_row(summary))
        result["status"] = "skipped_existing"
        return result

    if graph_path.exists() and args.force:
        graph_path.unlink()
    if summary_json.exists() and args.force:
        summary_json.unlink()

    with log_path.open("a") as log_handle, contextlib.redirect_stdout(log_handle), contextlib.redirect_stderr(
        log_handle
    ):
        print(f"[run] {params['run_name']}")
        print(f"[params] {json.dumps(params, sort_keys=True)}")

        builder = GraphBuilderProst(
            prost_tsv=args.prost_tsv,
            gff_dir=args.gff_dir,
            operonic_distance=int(params["operonic_distance"]),
            adjacency_only=args.adjacency_only,
            cluster_distance=float(params["stage3_distance"]),
            stage1_distance=float(params["stage1_distance"]),
            stage2_distance=float(params["stage2_distance"]),
            stage3_distance=float(params["stage3_distance"]),
            cluster_evalue=args.cluster_evalue,
            reciprocal_only=not args.no_reciprocal_only,
            reference_genome=args.reference_genome,
            preferred_gene_prefixes=[x.strip() for x in args.prefer_gene_prefixes.split(",") if x.strip()],
            enable_cluster_merge=args.enable_cluster_merge,
            enforce_one_per_genome=params["enforce_one_per_genome"],
            progress_every=args.progress_every,
        )

        try:
            graph = builder.run()
            print(f"[graph] nodes={graph.number_of_nodes()} edges={graph.number_of_edges()}")
            builder.export_gml(str(graph_path))
            summary = summarize_graph_file(
                graph_gml=graph_path,
                summary_json=summary_json,
                plots_dir=plots_dir,
                run_name=params["run_name"],
                parameters=params,
                top_hubs=args.top_hubs,
            )
            print(f"[wrote] {graph_path}")
            print(f"[wrote] {summary_json}")
            result.update(summary_to_master_row(summary))
            result["status"] = "success"
        except Exception as exc:
            traceback.print_exc()
            result["status"] = "failed"
            result["error"] = str(exc)

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Run a parameter sweep over operon graph clustering settings and summarize each graph."
    )
    parser.add_argument("--prost_tsv", required=True, help="PROST results TSV (combined.tsv).")
    parser.add_argument(
        "--gff_dir",
        required=True,
        help="Directory containing genome subfolders, each with a GFF file.",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory for run folders, master_summary.csv, and comparison_plots.",
    )
    parser.add_argument(
        "--one_per_genome_values",
        default="true,false",
        help="Comma-separated sweep values for the one-per-genome constraint.",
    )
    parser.add_argument(
        "--stage1_distance",
        type=int,
        default=1000,
        help="Legacy single stage 1 RBH seed distance. Used only if --stage1_distances is not provided.",
    )
    parser.add_argument(
        "--stage1_distances",
        default="",
        help="Comma-separated PROST distance cutoffs to sweep for stage 1 RBH seeding.",
    )
    parser.add_argument(
        "--stage2_distances",
        default="1000,3000,5000",
        help="Comma-separated PROST distance cutoffs for stage 2 expansion.",
    )
    parser.add_argument(
        "--stage3_distances",
        default="1000,3000,5000",
        help="Comma-separated PROST distance cutoffs for stage 3 consolidation.",
    )
    parser.add_argument(
        "--operonic_distance",
        type=int,
        default=40,
        help="Legacy single operonic distance. Used only if --operonic_distances is not provided.",
    )
    parser.add_argument(
        "--operonic_distances",
        default="",
        help="Comma-separated operonic-distance values to sweep. Defaults to the single --operonic_distance value.",
    )
    parser.add_argument(
        "--adjacency_only",
        action="store_true",
        help="Connect adjacent genes regardless of bp gap.",
    )
    parser.add_argument(
        "--cluster_evalue",
        type=float,
        default=1e-10,
        help="Maximum e-value to keep a PROST hit in clustering.",
    )
    parser.add_argument(
        "--no_reciprocal_only",
        action="store_true",
        help="Allow seed pairing from any thresholded hit instead of RBH-only.",
    )
    parser.add_argument(
        "--reference_genome",
        default=None,
        help="Optional genome name to prioritize for node labels.",
    )
    parser.add_argument(
        "--prefer_gene_prefixes",
        default="",
        help="Comma-separated preferred gene prefixes for node naming.",
    )
    parser.add_argument(
        "--enable_cluster_merge",
        dest="enable_cluster_merge",
        action="store_true",
        help="Enable stage 3 cluster-cluster consolidation. This is the default for the sweep.",
    )
    parser.add_argument(
        "--disable_cluster_merge",
        dest="enable_cluster_merge",
        action="store_false",
        help="Skip stage 3 cluster-cluster consolidation during the sweep.",
    )
    parser.set_defaults(enable_cluster_merge=True)
    parser.add_argument(
        "--progress_every",
        type=int,
        default=1_000_000,
        help="Accepted PROST hits between stage 0 progress messages.",
    )
    parser.add_argument("--top_hubs", type=int, default=15, help="Top hubs to retain in each summary.")
    parser.add_argument("--force", action="store_true", help="Rerun and overwrite existing outputs.")
    args = parser.parse_args()

    output_dir = ensure_dir(Path(args.output_dir))
    one_per_genome_values = parse_bool_list(args.one_per_genome_values)
    operonic_distances = (
        parse_int_list(args.operonic_distances) if args.operonic_distances else [int(args.operonic_distance)]
    )
    stage1_distances = (
        parse_int_list(args.stage1_distances) if args.stage1_distances else [int(args.stage1_distance)]
    )
    stage2_distances = parse_int_list(args.stage2_distances)
    stage3_distances = parse_int_list(args.stage3_distances)

    combinations = enumerate_parameter_grid(
        operonic_distances=operonic_distances,
        stage1_distances=stage1_distances,
        one_per_genome_values=one_per_genome_values,
        stage2_distances=stage2_distances,
        stage3_distances=stage3_distances,
    )

    expected_count = (
        len(operonic_distances)
        * len(stage1_distances)
        * len(one_per_genome_values)
        * len(stage2_distances)
        * len(stage3_distances)
    )
    if len(combinations) != expected_count:
        raise SystemExit(
            f"Combination count mismatch: enumerated {len(combinations)} but expected {expected_count}."
        )

    print(f"[sweep] total_runs={len(combinations)}")
    if (
        sorted(operonic_distances) == sorted(DEFAULT_OPERONIC_DISTANCES)
        and
        sorted(one_per_genome_values) == sorted(DEFAULT_ONE_PER_GENOME_VALUES)
        and sorted(stage1_distances) == [DEFAULT_STAGE1_DISTANCE]
        and sorted(stage2_distances) == sorted(DEFAULT_STAGE2_DISTANCES)
        and sorted(stage3_distances) == sorted(DEFAULT_STAGE3_DISTANCES)
        and len(combinations) != 18
    ):
        raise SystemExit(f"Expected exactly 18 runs for the default sweep, found {len(combinations)}.")

    all_rows = []
    successful_summaries = []
    for params in combinations:
        row = run_single_configuration(args=args, params=params, sweep_root=output_dir)
        all_rows.append(row)
        if row["status"] in {"success", "skipped_existing"} and Path(row["summary_json"]).exists():
            successful_summaries.append(json.loads(Path(row["summary_json"]).read_text()))
        print(f"[run_status] {row['run_name']} status={row['status']}")
        if row["error"]:
            print(f"[run_error] {row['run_name']} error={row['error']}")

    write_master_rows(all_rows, output_dir / "master_summary.csv")

    if successful_summaries:
        comparison_dir = output_dir / "comparison_plots"
        plot_all_heatmaps(successful_summaries, comparison_dir)
        plot_all_distribution_comparisons(successful_summaries, comparison_dir)
        plot_hub_heatmap(successful_summaries, comparison_dir / "heatmaps")

    print(f"[wrote] {output_dir / 'master_summary.csv'}")
    if successful_summaries:
        print(f"[wrote] {output_dir / 'comparison_plots'}")


if __name__ == "__main__":
    main()
