# proximon-graph
ProximonGraph constructs protein networks where nodes are homolog clusters and edges encode operonic adjacency, enabling detection of conserved gene neighborhoods and functional modules across genomes and metagenomes.

---

## Overview

This pipeline takes PROST clustering results and Prokka GFF annotations and produces a graph where **nodes are protein clusters** (ortholog groups) and **edges connect clusters whose members are genomically adjacent** (operonic neighbors). The resulting graph can be filtered, queried, laid out, and explored interactively.

**Core workflow:**

```
PROST TSV + GFF files
        ↓
proximon_graph_builder.py        →   graph.gml
        ↓
filter_graph_by_weight.py        →   high_weight.gml          (optional)
        ↓
marker_components_to_cyjs.py     →   subgraph.cyjs.json        (targeted query)
        ↓
layout_to_cyjs.py                →   layout.cyjs.json          (full graph layout)
        ↓
Cytoscape / Cytoscape.js
```

**Parameter sweep workflow** (to optimize clustering settings before committing to a single graph):

```
PROST TSV + GFF files
        ↓
run_proximon_graph_parameter_sweep.py   →   sweep_dir/run_*/graph.gml
                                             sweep_dir/master_summary.csv
                                             sweep_dir/comparison_plots/
        ↓
summarize_proximon_graph_parameter_sweep.py  →  re-aggregate or summarize a single graph
```

---

## 1. Environment Setup

### Prerequisites

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda

### Create the environment

```bash
conda create -n operon_graph python=3.10
conda activate operon_graph

# Core dependencies
pip install networkx pandas

# Optional: Graphviz-based layouts for large hairball components
conda install -c conda-forge graphviz pygraphviz
```

### Verify

```bash
python -c "import networkx, pandas; print('OK')"
```

---

## 2. Build the Graph — `proximon_graph_builder.py`

Builds a directed gene-cluster adjacency graph from PROST clustering output and Prokka GFF annotation folders.

**What it does:**
- Parses GFF files to extract CDS features and genomic positions
- Groups genes into operons by strand and intergenic distance
- Clusters proteins using RBH (reciprocal best hit) seed pairs with multi-stage expansion
- Creates graph nodes (one per cluster) and edges (one per adjacent cluster pair)
- Weights edges by the number of genomes supporting each adjacency

### Input files

| File | Description |
|------|-------------|
| `--prost_tsv` | PROST results TSV (`combined.tsv`), columns: head, member, distance, e-value |
| `--gff_dir` | Directory of genome subfolders, each containing a single `.gff` file |

**GFF directory structure:**

```
gff_dir/
├── genome_A/
│   └── genome_A.gff
├── genome_B/
│   └── genome_B.gff
└── ...
```

### Run

```bash
python proximon_graph_builder.py \
  --prost_tsv ./data/prost_results/combined.tsv \
  --gff_dir ./data/prokka_results/genomes/ \
  --out_gml ./graphs/my_graph.gml \
  --operonic_distance 100 \
  --stage1_distance 3000 \
  --stage2_distance 3000 \
  --stage3_distance 3000 \
  --enable_cluster_merge
```

### Example call (from a real run)

```bash
python proximon_graph_builder.py \
  --prost_tsv ./data/prost_results/nylon_combined_reference.tsv \
  --gff_dir ./data/prokka_results/nyl_genomes \
  --out_gml ./graphs/16MAGS_nylon_3k_3k_3k_100bp.gml \
  --operonic_distance 100 \
  --stage1_distance 3000 \
  --stage2_distance 3000 \
  --stage3_distance 3000 \
  --enable_cluster_merge \
  --no_enforce_one_per_genome
```

### Key arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--operonic_distance` | `40` | Max intergenic bp to connect adjacent genes |
| `--adjacency_only` | off | Connect all adjacent genes regardless of gap |
| `--stage1_distance` | `cluster_distance` | PROST distance cutoff for RBH seed selection |
| `--stage2_distance` | `cluster_distance` | PROST distance cutoff for cluster expansion |
| `--stage3_distance` | `cluster_distance` | PROST distance cutoff for cluster-cluster merge |
| `--cluster_evalue` | `1e-10` | Max e-value for clustering |
| `--enable_cluster_merge` | off | Enable the final cluster-cluster consolidation pass |
| `--enforce_one_per_genome` | on | Keep at most one protein per genome per cluster |
| `--reference_genome` | — | Genome to use for reference annotations on nodes |
| `--prefer_gene_prefixes` | — | Comma-separated gene name prefixes to prioritize for node labels |

### Differential expression overlay (optional)

You can overlay per-gene DE values onto graph nodes:

```bash
  --de_csv results.csv \
  --de_gene_col Geneid \
  --de_log2fc_col log2FoldChange \
  --de_padj_col padj \
  --de_padj_threshold 0.05 \
  --de_abs_log2fc_threshold 1.0
```

Node attributes added: `node_log2FC`, `node_min_padj`, `node_num_DE`, `node_frac_DE`, `node_mean_expr`.

### Outputs

| File | Description |
|------|-------------|
| `*.gml` | Main graph file for all downstream steps |
| `*_nodes.tsv` | Node table (with `--out_nodes_tsv`) |
| `*_edges.tsv` | Edge table (with `--out_edges_tsv`) |

---

## 3. Filter for Core Operons (Optional) — `filter_graph_by_weight.py`

Retains only high-confidence edges, isolates core operons, and reduces noise before visualization.

```bash
python filter_graph_by_weight.py \
  --input_gml full_graph.gml \
  --output_gml high_weight.gml \
  --min_weight 3 \
  --weight_key abs_weight
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--min_weight` | `3.0` | Minimum edge weight to keep |
| `--weight_key` | `abs_weight` | Use `abs_weight` (raw count) or `norm_weight` (cluster-size normalized) |

Isolated nodes (no remaining edges) are removed automatically.

---

## 4. Query Targeted Subgraphs — `marker_components_to_cyjs.py`

Tags graph nodes from a marker TSV (a list of genes of interest) and extracts a focused subgraph for visualization. Useful for exploring the genomic neighborhood of specific genes or pathways.

### Marker TSV format

The TSV must have at least these three columns (see `examples/example_markers.tsv`):

```
raw_id          label                       role
lacZ            beta-galactosidase          enzyme
trpA            tryptophan synthase subunit A   enzyme
nylB            nylon oligomer hydrolase    enzyme
```

`raw_id` must match the GFF `ID` field (or a member protein ID in the graph).

### Extraction modes

**Component mode** (default): extract the entire connected component(s) containing each matched node.

```bash
python marker_components_to_cyjs.py \
  --gml your_graph.gml \
  --markers_tsv markers.tsv \
  --out_prefix out/my_markers \
  --selection_mode component
```

**Hop mode**: extract an n-hop neighborhood around matched nodes.

```bash
python marker_components_to_cyjs.py \
  --gml your_graph.gml \
  --markers_tsv markers.tsv \
  --out_prefix out/markers_2hop \
  --selection_mode hops \
  --hops 2
```

### Key arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--selection_mode` | `component` | `component` or `hops` |
| `--hops` | `1` | Neighborhood radius when using `--selection_mode hops` |
| `--directed` | `either` | Hop traversal direction: `either`, `out`, or `in` |
| `--id_column` | `raw_id` | Column in marker TSV to match against node/member IDs |
| `--write_gml` | off | Also write the extracted subgraph as GML |
| `--write_tagged_full_gml` | off | Also write the full graph with marker annotations as GML |
| `--write_unmatched_ids` | off | Write marker IDs that did not match any graph node |

### Node attributes added

| Attribute | Description |
|-----------|-------------|
| `marker_hit` | `1` if this node matched a marker ID |
| `marker_seed` | `1` for matched seed nodes |
| `marker_raw_ids` | Semicolon-separated matched marker IDs |
| `marker_label` / `marker_role` | From marker TSV |
| `marker_class` | `marker_hit` or `context` |
| `marker_seed_distance` | Hop distance to the nearest seed node |

### Outputs

| File | Description |
|------|-------------|
| `*.cyjs.json` | Cytoscape.js JSON ready for visualization |
| `*.positions.tsv` | Node coordinates for Cytoscape Desktop import |

---

## 5. Full-Graph Layout and Export — `layout_to_cyjs.py`

Computes optimized layouts for the full graph (or any filtered GML), handles large "hairball" components separately, and exports Cytoscape.js JSON.

```bash
python layout_to_cyjs.py \
  --gml graph.gml \
  --out_prefix out/my_layout \
  --scale 500 \
  --padding 180 \
  --target_aspect 1.2
```

### Hairball handling

If the largest connected component exceeds `--hairball_min_nodes` (default: 300), the script splits the graph:

- **hairball** — the large component, laid out with hub-relaxed spring or SFDP
- **cruise** — all remaining smaller components, laid out with spring per component
- **allinone** — both packed together with the hairball at top-left

When no hairball is detected, a single layout file is written.

### Key layout arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--scale` | `450.0` | Layout scale for small components |
| `--padding` | `180.0` | Padding between packed components |
| `--target_aspect` | `1.2` | Desired width/height ratio |
| `--hairball_min_nodes` | `300` | Component size threshold for hairball treatment |
| `--hairball_layout` | `hub_relaxed_spring` | `hub_relaxed_spring` or `sfdp` (requires pygraphviz) |
| `--hairball_scale` | `1400.0` | Scale for hairball component |
| `--max_k` | `0` | Limit edges per node to top-K (0 = no limit) |
| `--drop_singletons` | off | Remove isolated nodes before layout |
| `--drop_doublets` | off | Remove 2-node components before layout |
| `--compact_passes` | `25` | Post-layout edge compaction passes |

### Outputs

| File | Description |
|------|-------------|
| `*_allinone.cyjs.json` | Full layout combining all components |
| `*_hairball.cyjs.json` | Large component only |
| `*_cruise.cyjs.json` | All smaller components |
| `*.positions.tsv` | Node coordinate TSV for each of the above |

---

## 6. Load into Cytoscape Desktop

1. **Download Cytoscape** at [https://cytoscape.org](https://cytoscape.org)

2. **Import network**
   `File → Import → Network from File`
   Select `*_allinone.cyjs.json` (or any `.cyjs.json` output)

3. **Import pre-computed layout** (optional, recommended for large graphs)
   `File → Import → Table from File`
   Select `*.positions.tsv`
   Map: `x → position X`, `y → position Y`

4. **Suggested visual mappings**

   | Visual property | Node attribute |
   |-----------------|---------------|
   | Node size | `cluster_size` |
   | Node label | `label` |
   | Node color | `marker_class` (highlight hits vs. context) |
   | Edge width | `abs_weight` |
   | Edge opacity | `norm_weight` |

5. **Example session**
   A pre-built Cytoscape session file is included at `examples/example_session.cys`. Open it with `File → Open` to explore a sample network with styles already applied.

---

## Examples

The `examples/` directory contains:

| File | Description |
|------|-------------|
| `example_markers.tsv` | Sample marker TSV showing the required format |
| `example_session.cys` | Pre-built Cytoscape session to explore |

---

## Output File Reference

| File | Script | Description |
|------|--------|-------------|
| `*.gml` | `proximon_graph_builder.py` | Full annotated graph |
| `*_filtered.gml` | `filter_graph_by_weight.py` | High-confidence edge subgraph |
| `*.cyjs.json` | `marker_components_to_cyjs.py` / `layout_to_cyjs.py` | Cytoscape.js JSON |
| `*.positions.tsv` | All visualization scripts | Node x/y coordinates |
| `*_nodes.tsv` | `proximon_graph_builder.py` | Node attribute table |
| `*_edges.tsv` | `proximon_graph_builder.py` | Edge attribute table |
| `master_summary.csv` | `run_proximon_graph_parameter_sweep.py` | Scalar metrics for all sweep runs |
| `summary.json` | `summarize_proximon_graph_parameter_sweep.py` | Per-run graph metrics JSON |
| `comparison_plots/` | Both sweep scripts | Heatmaps, delta heatmaps, ECDFs, histograms |

---

## 7. Parameter Sweep — `run_proximon_graph_parameter_sweep.py`

Before committing to a single set of clustering parameters, you can systematically sweep over combinations of operonic distance, PROST distance thresholds, and the one-per-genome constraint. Each combination produces its own graph and summary, and results are aggregated into a single CSV and a set of comparison plots.

### What it sweeps

| Axis | Argument | Default sweep values |
|------|----------|---------------------|
| Operonic distance | `--operonic_distances` | `40` |
| Stage 1 RBH seed distance | `--stage1_distances` | `1000` |
| Stage 2 expansion distance | `--stage2_distances` | `1000, 3000, 5000` |
| Stage 3 merge distance | `--stage3_distances` | `1000, 3000, 5000` |
| One protein per genome | `--one_per_genome_values` | `true, false` |

The default configuration runs **18 parameter combinations** (2 × 1 × 3 × 3).

### Run

```bash
python run_proximon_graph_parameter_sweep.py \
  --prost_tsv ./data/prost_results/combined.tsv \
  --gff_dir ./data/prokka_results/genomes/ \
  --output_dir ./sweep_results/ \
  --stage2_distances 1000,3000,5000 \
  --stage3_distances 1000,3000,5000 \
  --one_per_genome_values true,false \
  --enable_cluster_merge
```

To sweep additional stage 1 distances or operonic distances:

```bash
  --stage1_distances 1000,3000 \
  --operonic_distances 40,100
```

### Key arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--output_dir` | required | Root directory; each run lands in `run_<name>/` |
| `--stage2_distances` | `1000,3000,5000` | Comma-separated stage 2 distance values to sweep |
| `--stage3_distances` | `1000,3000,5000` | Comma-separated stage 3 distance values to sweep |
| `--stage1_distances` | `1000` | Comma-separated stage 1 distance values to sweep |
| `--operonic_distances` | `40` | Comma-separated operonic distance values to sweep |
| `--one_per_genome_values` | `true,false` | Whether to enforce one protein per genome per cluster |
| `--enable_cluster_merge` | on | Enable stage 3 cluster-cluster consolidation (default for sweep) |
| `--force` | off | Rerun and overwrite existing outputs |
| `--top_hubs` | `15` | Top hub nodes to record per run |

### Outputs

```
sweep_results/
├── master_summary.csv
├── comparison_plots/
│   ├── heatmaps/
│   ├── delta_heatmaps/
│   └── distributions/
│       ├── ecdf/
│       └── histograms/
└── run_operondist_40_stage1_1000_stage2_3000_stage3_3000_onepergenome_true/
    ├── graph.gml
    ├── summary.json
    ├── log.txt
    └── plots/
        ├── distributions.png
        └── top_hubs.png
```

`master_summary.csv` contains one row per run with columns for all parameter settings and graph metrics: `total_nodes`, `total_edges`, `compression_ratio`, `largest_component_size`, `largest_component_fraction`, `singleton_node_fraction`, `number_of_components`, `mean_degree`, `degree_p95`, `degree_p99`, `proteins_per_node_mean`, `top_hub_label`, `top_hub_degree`, and all parameter columns.

Comparison plots are produced for every combination of `enforce_one_per_genome`, `operonic_distance`, and `stage1_distance`, so each heatmap shows a clean 2D grid of stage2 vs stage3 values. Runs that already have both `graph.gml` and `summary.json` are skipped automatically — use `--force` to rerun.

---

## 8. Summarize a Single Graph or Re-aggregate a Sweep — `summarize_proximon_graph_parameter_sweep.py`

This script has two modes: summarize a single graph GML in isolation, or re-aggregate all `summary.json` files from an existing sweep directory. The second mode is useful for regenerating plots after a sweep without re-running the graph builder.

### Summarize a single graph

```bash
python summarize_proximon_graph_parameter_sweep.py \
  --graph_gml ./graphs/my_graph.gml \
  --summary_json ./graphs/my_graph_summary.json \
  --plots_dir ./graphs/plots/ \
  --run_name my_run
```

### Re-aggregate an existing sweep

```bash
python summarize_proximon_graph_parameter_sweep.py \
  --sweep_dir ./sweep_results/ \
  --master_summary_csv ./sweep_results/master_summary_v2.csv
```

Reads all `run_*/summary.json` files under `--sweep_dir` and regenerates `master_summary.csv` plus all comparison plots, without rebuilding any graphs.

### Key arguments

| Argument | Description |
|----------|-------------|
| `--graph_gml` | Single graph to summarize (mutually exclusive with `--sweep_dir`) |
| `--summary_json` | Output path for the single-run JSON summary |
| `--plots_dir` | Directory for distribution and hub plots |
| `--sweep_dir` | Sweep root to re-aggregate (mutually exclusive with `--graph_gml`) |
| `--master_summary_csv` | Output CSV path (defaults to `sweep_dir/master_summary.csv`) |
| `--top_hubs` | Top hub nodes to keep per run (default: `15`) |

### Per-run summary metrics

Each `summary.json` (and each row of `master_summary.csv`) captures graph metrics (total nodes, edges, compression ratio, density, clustering coefficient), component structure (number of components, largest component size and fraction, singleton fraction, size histogram), proteins-per-node distribution (mean, median, std, max), node degree distribution (mean, median, max, p95, p99), and a ranked table of top hub nodes with their labels, degrees, and product annotations.

---

## Tips

- **Start with `--operonic_distance 100`** for a well-connected graph; decrease to `40` for strict operon-only edges.
- **For very large runs**, leave `--enable_cluster_merge` off (it is expensive) and use `--no_enforce_one_per_genome` only if paralogs are relevant to your analysis.
- **`--prefer_gene_prefixes`** is useful when well-characterized gene names (e.g., `nyl`, `lac`, `trp`) are present in only a subset of genomes — it biases node labels toward those names for readability.
- **`--hairball_layout sfdp`** produces better separation in very dense components but requires pygraphviz (`conda install -c conda-forge graphviz pygraphviz`).
- **Run the parameter sweep first** when working with a new dataset. The `master_summary.csv` and heatmaps make it easy to identify which stage2/stage3 distance combination gives the best balance of compression ratio, component structure, and singleton fraction before running the full visualization pipeline.
