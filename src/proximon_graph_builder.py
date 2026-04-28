import glob
import math
import os
import re
from collections import Counter
from time import perf_counter
from typing import Dict, List, Set, Tuple

import networkx as nx
import pandas as pd

from rbh_seed_utils import (
    compute_rbh_seed_state,
    passes_margin,
    support_score_key,
)


class GraphBuilderProst:
    """
    Builds a gene-cluster adjacency graph from:
    - PROST clustering results
    - GFF-derived operons

    GFF ID is assumed to match PROST FASTA headers exactly.
    """

    def __init__(
        self,
        prost_tsv: str,
        gff_dir: str,
        operonic_distance: int = 40,
        adjacency_only: bool = False,
        cluster_distance: float = 5000.0,
        stage1_distance: float = None,
        stage2_distance: float = None,
        stage3_distance: float = None,
        stage12_cluster_distance: float = None,
        stage3_cluster_distance: float = None,
        cluster_evalue: float = 1e-10,
        reciprocal_only: bool = True,
        reference_genome: str = None,
        preferred_gene_prefixes: List[str] = None,
        de_csv: str = None,
        de_gene_col: str = "Geneid",
        de_log2fc_col: str = None,
        de_padj_col: str = None,
        de_expr_col: str = None,
        de_padj_threshold: float = 0.05,
        de_abs_log2fc_threshold: float = 1.0,
        de_expression_agg: str = "mean",
        enable_cluster_merge: bool = False,
        enforce_one_per_genome: bool = True,
        progress_every: int = 1_000_000,
    ):
        self.prost_df = self._load_and_clean_prost(prost_tsv)
        self.gff_dir = gff_dir
        self.operonic_distance = operonic_distance
        self.adjacency_only = adjacency_only

        self.cluster_distance = cluster_distance
        resolved_stage1 = stage1_distance
        if resolved_stage1 is None:
            resolved_stage1 = stage12_cluster_distance
        if resolved_stage1 is None:
            resolved_stage1 = cluster_distance

        resolved_stage2 = stage2_distance
        if resolved_stage2 is None:
            resolved_stage2 = stage12_cluster_distance
        if resolved_stage2 is None:
            resolved_stage2 = cluster_distance

        resolved_stage3 = stage3_distance
        if resolved_stage3 is None:
            resolved_stage3 = stage3_cluster_distance
        if resolved_stage3 is None:
            resolved_stage3 = cluster_distance

        self.stage1_distance = resolved_stage1
        self.stage2_distance = resolved_stage2
        self.stage3_distance = resolved_stage3
        self.cluster_evalue = cluster_evalue
        self.reciprocal_only = reciprocal_only

        self.genes: List[Dict] = []
        self.operons: List[List[Dict]] = []

        self.head_to_members: Dict[str, Set[str]] = {}
        self.member_to_cluster: Dict[str, str] = {}
        self.gff_to_prost_map: Dict[str, str] = {}

        self.nodes: Set[str] = set()
        self.gene_pairs: List[Tuple[str, str, str, int]] = []
        self.edges: List[Tuple[str, str, str, str, str, int]] = []
        self.edge_properties: Dict[Tuple[str, str], Dict] = {}
        self.weighted_edges: List[Tuple[str, str, Dict]] = []

        self.reference_genome = reference_genome
        self.preferred_gene_prefixes = [p.lower() for p in (preferred_gene_prefixes or []) if p]
        self.de_csv = de_csv
        self.de_gene_col = de_gene_col
        self.de_log2fc_col = de_log2fc_col
        self.de_padj_col = de_padj_col
        self.de_expr_col = de_expr_col
        self.de_padj_threshold = de_padj_threshold
        self.de_abs_log2fc_threshold = de_abs_log2fc_threshold
        self.de_expression_agg = de_expression_agg
        self.de_by_gene = self._load_de_table(de_csv) if de_csv else {}
        self.enable_cluster_merge = enable_cluster_merge
        self.enforce_one_per_genome = enforce_one_per_genome
        self.progress_every = max(1, int(progress_every))

        self.G = nx.DiGraph()

    # ---------------- PROST ----------------
    def _load_and_clean_prost(self, prost_tsv: str) -> pd.DataFrame:
        df = pd.read_csv(
            prost_tsv,
            sep="\t",
            header=None,
            usecols=[0, 1, 5, 6],
            names=["head", "member", "distance", "e_value"],
        )

        df["head"] = df["head"].astype(str).str.split().str[0]
        df["member"] = df["member"].astype(str).str.split().str[0]
        df["distance"] = pd.to_numeric(df["distance"], errors="coerce")
        df["e_value"] = pd.to_numeric(df["e_value"], errors="coerce")

        df = df.dropna(subset=["head", "member", "distance", "e_value"])
        return df

    @staticmethod
    def _first_existing_column(df: pd.DataFrame, candidates: List[str]) -> str:
        for col in candidates:
            if col in df.columns:
                return col
        return ""

    @staticmethod
    def _clean_number(value):
        if value is None or pd.isna(value):
            return None
        try:
            number = float(value)
        except (TypeError, ValueError):
            return None
        if not math.isfinite(number):
            return None
        return number

    def _resolve_de_columns(self, df: pd.DataFrame):
        gene_col = self.de_gene_col or "Geneid"
        if gene_col not in df.columns:
            raise ValueError(
                f"DE CSV missing gene column '{gene_col}'. Found columns: {list(df.columns)}"
            )

        log2fc_col = self.de_log2fc_col or self._first_existing_column(
            df,
            [
                "log2FoldChange",
                "log2FC",
                "logFC",
                "Nylon-6,6",
            ],
        )
        padj_col = self.de_padj_col or self._first_existing_column(
            df,
            [
                "padj",
                "FDR",
                "adj.P.Val",
                "qvalue",
                "q_value",
            ],
        )
        expr_col = self.de_expr_col or self._first_existing_column(
            df,
            [
                "AveExpr",
                "baseMean",
                "mean_expr",
                "expression",
            ],
        )

        missing = []
        if not log2fc_col:
            missing.append("log2FC")
        if not padj_col:
            missing.append("adjusted p-value")
        if not expr_col:
            missing.append("expression")
        if missing:
            raise ValueError(
                "Could not infer DE CSV column(s): "
                + ", ".join(missing)
                + ". Set --de_log2fc_col, --de_padj_col, or --de_expr_col explicitly."
            )
        return gene_col, log2fc_col, padj_col, expr_col

    def _load_de_table(self, de_csv: str) -> Dict[str, Dict[str, float]]:
        """Load per-gene differential-expression values keyed by gene id."""
        df = pd.read_csv(de_csv)
        gene_col, log2fc_col, padj_col, expr_col = self._resolve_de_columns(df)

        de_by_gene = {}
        duplicate_genes = 0
        for _, row in df.iterrows():
            gene_id = str(row[gene_col]).strip()
            if not gene_id or gene_id.lower() == "nan":
                continue

            if gene_id in de_by_gene:
                duplicate_genes += 1
                continue

            de_by_gene[gene_id] = {
                "log2fc": self._clean_number(row[log2fc_col]),
                "padj": self._clean_number(row[padj_col]),
                "expr": self._clean_number(row[expr_col]),
            }

        print(
            f"[load_de_table] genes={len(de_by_gene)} "
            f"duplicate_genes_skipped={duplicate_genes} "
            f"log2fc_col={log2fc_col} padj_col={padj_col} expr_col={expr_col}"
        )
        return de_by_gene

    def _de_attributes_for_members(self, members: Set[str]) -> Dict[str, object]:
        matched = [self.de_by_gene[m] for m in members if m in self.de_by_gene]

        log2fcs = [row["log2fc"] for row in matched if row["log2fc"] is not None]
        padjs = [row["padj"] for row in matched if row["padj"] is not None]
        exprs = [row["expr"] for row in matched if row["expr"] is not None]

        de_count = 0
        for row in matched:
            padj = row["padj"]
            log2fc = row["log2fc"]
            if padj is None or log2fc is None:
                continue
            if padj <= self.de_padj_threshold and abs(log2fc) >= self.de_abs_log2fc_threshold:
                de_count += 1

        attrs = {
            "node_log2FC": "",
            "node_min_padj": "",
            "node_num_DE": de_count,
            "node_frac_DE": de_count / len(members) if members else 0.0,
            "node_mean_expr": "",
            "node_mean_expr_gene_count": len(exprs),
        }
        if log2fcs:
            attrs["node_log2FC"] = float(pd.Series(log2fcs).median())
        if padjs:
            attrs["node_min_padj"] = min(padjs)
        if exprs:
            if self.de_expression_agg == "sum":
                attrs["node_mean_expr"] = sum(exprs)
            else:
                attrs["node_mean_expr"] = sum(exprs) / len(exprs)
        return attrs

    # ---------------- reset before running graph ----------------
    def reset(self):
        self.genes = []
        self.operons = []

        self.head_to_members = {}
        self.member_to_cluster = {}
        self.gff_to_prost_map = {}

        self.nodes = set()
        self.gene_pairs = []
        self.edges = []
        self.edge_properties = {}
        self.weighted_edges = []

        self.G = nx.DiGraph()

    # ---------------- GFF parsing ----------------
    @staticmethod
    def _derive_genome_name(folder_genome_name: str, raw_id: str) -> str:
        """Infer the genome label from GFF IDs that already encode genome identity."""
        if "|" in raw_id:
            return raw_id.split("|", 1)[0]

        numeric_suffix_match = re.match(r"(.+)_\d+$", raw_id)
        if numeric_suffix_match:
            return numeric_suffix_match.group(1)

        return folder_genome_name

    def parse_gffs(self) -> List[Dict]:
        """
        Parse GFFs and extract CDS features.
        Uses GFF 'ID' as the PROST-compatible identifier.
        """
        all_genes = []

        genome_folders = sorted(glob.glob(os.path.join(self.gff_dir, "*")))
        for genome_path in genome_folders:
            genome_name = os.path.basename(genome_path)
            gff_files = sorted(glob.glob(os.path.join(genome_path, "*.gff")))

            if not gff_files:
                continue

            gff_file = gff_files[0]

            with open(gff_file) as handle:
                for line in handle:
                    if line.startswith("#"):
                        continue

                    parts = line.rstrip().split("\t")
                    if len(parts) < 9 or parts[2] != "CDS":
                        continue

                    seqid = parts[0]
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = parts[6]

                    attr = {}
                    for field in parts[8].split(";"):
                        if "=" in field:
                            key, value = field.split("=", 1)
                            attr[key] = value

                    raw_id = attr.get("ID")
                    if raw_id is None:
                        continue

                    product = attr.get("product", "hypothetical protein")
                    gene_name = attr.get("gene") or attr.get("Name") or None
                    derived_genome_name = self._derive_genome_name(genome_name, raw_id)
                    orf_id = raw_id

                    gene = {
                        "raw_id": raw_id,
                        "orf_id": orf_id,
                        "contig": seqid,
                        "start": min(start, end),
                        "end": max(start, end),
                        "strand": strand,
                        "product": product,
                        "gene_name": gene_name,
                        "genome": derived_genome_name,
                    }

                    all_genes.append(gene)
                    self.gff_to_prost_map[raw_id] = orf_id
                    self.gff_to_prost_map[orf_id] = orf_id

        self.genes = all_genes
        self.orf_lookup = {g["orf_id"]: g for g in self.genes}
        print(f"[parse_gffs] Parsed {len(all_genes)} CDS features")
        return all_genes

    # ---------------- Operons ----------------
    @staticmethod
    def gene_distance(g1: Dict, g2: Dict):
        if g1["contig"] != g2["contig"] or g1["strand"] != g2["strand"]:
            return None
        if g1["strand"] == "+":
            return g2["start"] - g1["end"]
        return g1["start"] - g2["end"]

    def group_operons(self):
        genes_sorted = sorted(
            self.genes,
            key=lambda g: (
                g["genome"],
                g["contig"],
                g["strand"],
                g["start"] if g["strand"] == "+" else -g["start"],
            ),
        )
        if not genes_sorted:
            return "No genes present"

        current = [genes_sorted[0]]
        for gene in genes_sorted[1:]:
            dist = self.gene_distance(current[-1], gene)
            if dist is not None and (self.adjacency_only or dist <= self.operonic_distance):
                current.append(gene)
            else:
                self.operons.append(current)
                current = [gene]
        self.operons.append(current)

    def load_operons(self):
        self.parse_gffs()
        self.group_operons()

    # ---------------- Clusters ---------------
    def resolve_clusters(self):
        """
        Scalable clustering for large PROST runs.

        Key differences from the original compatibility builder:
        - vectorized PROST filtering before Python loops
        - `itertuples()` instead of `iterrows()`
        - cached cluster genome sets during expansion
        - precomputed reciprocal-neighbor maps for fast candidate screening
        - optional disabling of the expensive final cluster merge stage
        - stage timing/progress logs for long runs
        """
        t0 = perf_counter()

        proteins = set(g["raw_id"] for g in self.genes)
        proteins.update(self.prost_df["head"].tolist())
        proteins.update(self.prost_df["member"].tolist())

        if not proteins:
            self.head_to_members = {}
            self.member_to_cluster = {}
            return

        protein_to_genome = {g["raw_id"]: g["genome"] for g in self.genes}
        stage1_distance_cutoff = self.stage1_distance
        stage2_distance_cutoff = self.stage2_distance
        stage3_distance_cutoff = self.stage3_distance
        initial_distance_cutoff = max(stage1_distance_cutoff, stage2_distance_cutoff, stage3_distance_cutoff)
        expansion_distance_cutoff = stage2_distance_cutoff
        expansion_evalue_cutoff = self.cluster_evalue
        merge_distance_cutoff = stage3_distance_cutoff
        merge_evalue_cutoff = self.cluster_evalue

        def log_stage(stage_name: str, start_time: float):
            elapsed = perf_counter() - start_time
            total = perf_counter() - t0
            print(f"[resolve_clusters:{stage_name}] elapsed={elapsed:.1f}s total={total:.1f}s")

        # Stages 0-1: thresholded hit collection and RBH seed selection
        s = perf_counter()
        seed_state = compute_rbh_seed_state(
            prost_df=self.prost_df,
            filter_distance_threshold=initial_distance_cutoff,
            seed_distance_threshold=stage1_distance_cutoff,
            evalue_threshold=self.cluster_evalue,
            proteins=proteins,
            protein_to_genome=protein_to_genome,
            reciprocal_only=self.reciprocal_only,
            enforce_one_per_genome=self.enforce_one_per_genome,
            progress_every=self.progress_every,
            progress_callback=lambda processed, unique_queries: print(
                f"[resolve_clusters:stage0] processed_hits={processed} "
                f"unique_queries={unique_queries}"
            ),
        )
        filtered = seed_state["filtered_df"]
        directed_hits = seed_state["directed_hits"]
        accepted_hits = seed_state["accepted_hits"]
        best_hit = seed_state["best_hit"]
        second_hit = seed_state["second_hit"]
        seed_pairs = seed_state["seed_pairs"]
        skipped_self = seed_state["filter_stats"]["skipped_self"]
        skipped_threshold = seed_state["filter_stats"]["skipped_threshold"]
        skipped_same_genome = seed_state["filter_stats"]["skipped_same_genome"]
        kept_rows = seed_state["filter_stats"]["accepted_directed_hits"]
        log_stage("stage0_stage1_seed_setup", s)

        def support_score(protein: str):
            return support_score_key(best_hit, protein)

        def passes_seed_margin(protein: str) -> bool:
            return passes_margin(protein, best_hit, second_hit)

        clusters = [set(pair) for pair in seed_pairs]
        cluster_genome_sets = [
            {protein_to_genome.get(member, f"UNKNOWN::{member}") for member in cluster}
            for cluster in clusters
        ]
        assigned = {member for pair in seed_pairs for member in pair}

        # Stage 2a: compatibility neighbor maps
        s = perf_counter()
        expansion_neighbors = {}
        merge_neighbors = {}
        for (a, b), hit_ab in accepted_hits.items():
            hit_ba = accepted_hits.get((b, a))
            if hit_ba is None:
                continue

            if (
                hit_ab["distance"] <= expansion_distance_cutoff
                and hit_ba["distance"] <= expansion_distance_cutoff
                and hit_ab["e_value"] <= expansion_evalue_cutoff
                and hit_ba["e_value"] <= expansion_evalue_cutoff
            ):
                expansion_neighbors.setdefault(a, set()).add(b)

            if self.enable_cluster_merge and (
                hit_ab["distance"] <= merge_distance_cutoff
                and hit_ba["distance"] <= merge_distance_cutoff
                and hit_ab["e_value"] <= merge_evalue_cutoff
                and hit_ba["e_value"] <= merge_evalue_cutoff
            ):
                merge_neighbors.setdefault(a, set()).add(b)

        def candidate_cluster_score(candidate: str, cluster: Set[str]):
            distances = []
            evalues = []
            for member in cluster:
                hit_1 = accepted_hits[(candidate, member)]
                hit_2 = accepted_hits[(member, candidate)]
                distances.append(max(hit_1["distance"], hit_2["distance"]))
                evalues.append(max(hit_1["e_value"], hit_2["e_value"]))
            worst_distance = max(distances) if distances else float("inf")
            worst_evalue = max(evalues) if evalues else float("inf")
            return (worst_distance, worst_evalue, candidate)

        def compatible_cluster_ids_for_candidate(candidate: str):
            candidate_genome = protein_to_genome.get(candidate, f"UNKNOWN::{candidate}")
            neighbors = expansion_neighbors.get(candidate, set())
            if not neighbors:
                return []

            compatible_ids = []
            for idx, cluster in enumerate(clusters):
                if self.enforce_one_per_genome and candidate_genome in cluster_genome_sets[idx]:
                    continue
                if not (cluster & neighbors):
                    continue
                if cluster.issubset(neighbors):
                    compatible_ids.append(idx)
            return compatible_ids

        def expand_clusters(unassigned: List[str], round_name: str):
            expansions_local = 0
            changed = True
            pass_num = 0

            while changed:
                pass_num += 1
                changed = False
                next_unassigned = []

                for j, protein in enumerate(unassigned, start=1):
                    compatible_ids = compatible_cluster_ids_for_candidate(protein)

                    if len(compatible_ids) == 1:
                        idx = compatible_ids[0]
                    elif len(compatible_ids) > 1:
                        ranked = sorted(
                            (candidate_cluster_score(protein, clusters[idx]), idx)
                            for idx in compatible_ids
                        )
                        idx = ranked[0][1]
                    else:
                        next_unassigned.append(protein)
                        continue

                    clusters[idx].add(protein)
                    cluster_genome_sets[idx].add(
                        protein_to_genome.get(protein, f"UNKNOWN::{protein}")
                    )
                    assigned.add(protein)
                    expansions_local += 1
                    changed = True

                    if expansions_local and expansions_local % 10000 == 0:
                        print(
                            f"[resolve_clusters:{round_name}] expansions={expansions_local} "
                            f"remaining={len(next_unassigned) + (len(unassigned) - j)} "
                            f"pass={pass_num}"
                        )

                unassigned = next_unassigned
                print(
                    f"[resolve_clusters:{round_name}] pass={pass_num} "
                    f"remaining_unassigned={len(unassigned)} expansions={expansions_local}"
                )

            return expansions_local, unassigned

        log_stage("stage2a_neighbors", s)

        # Stage 2: expansion
        s = perf_counter()
        expansions = 0
        initial_unassigned = sorted((p for p in proteins if p not in assigned), key=support_score)
        exp_1, unassigned = expand_clusters(initial_unassigned, "stage2_expand")
        expansions += exp_1
        log_stage("stage2_expand", s)

        # Stage 2b: leftover seeds + second expansion
        s = perf_counter()
        leftover_pairs = []
        leftover_seen = set()

        if self.reciprocal_only:
            for a in unassigned:
                if a not in best_hit:
                    continue
                b = best_hit[a]["target"]
                if b not in unassigned or b not in best_hit:
                    continue
                if best_hit[b]["target"] != a:
                    continue
                if not passes_seed_margin(a) or not passes_seed_margin(b):
                    continue
                pair = tuple(sorted((a, b)))
                if pair not in leftover_seen:
                    leftover_seen.add(pair)
                    leftover_pairs.append(pair)
        else:
            unused = set(unassigned)
            for a in sorted(unused, key=support_score):
                if a not in unused or a not in best_hit:
                    continue
                b = best_hit[a]["target"]
                if b not in unused or a == b:
                    continue
                pair = tuple(sorted((a, b)))
                if pair not in leftover_seen:
                    leftover_seen.add(pair)
                    leftover_pairs.append(pair)
                    unused.discard(a)
                    unused.discard(b)

        for a, b in leftover_pairs:
            if a in assigned or b in assigned:
                continue
            clusters.append({a, b})
            cluster_genome_sets.append(
                {
                    protein_to_genome.get(a, f"UNKNOWN::{a}"),
                    protein_to_genome.get(b, f"UNKNOWN::{b}"),
                }
            )
            assigned.add(a)
            assigned.add(b)

        second_unassigned = sorted((p for p in proteins if p not in assigned), key=support_score)
        exp_2, unassigned = expand_clusters(second_unassigned, "stage2b_expand")
        expansions += exp_2
        log_stage("stage2b_expand", s)

        # Stage 3: optional cluster-cluster merge
        s = perf_counter()
        cluster_merges = 0
        if self.enable_cluster_merge and clusters:
            changed = True
            while changed:
                changed = False
                member_to_cluster_ids = {}
                for idx, cluster in enumerate(clusters):
                    for member in cluster:
                        member_to_cluster_ids.setdefault(member, set()).add(idx)

                used = [False] * len(clusters)
                new_clusters = []
                new_cluster_genome_sets = []

                for i in range(len(clusters)):
                    if used[i]:
                        continue

                    current = set(clusters[i])
                    current_genomes = set(cluster_genome_sets[i])
                    used[i] = True

                    grew = True
                    while grew:
                        grew = False
                        merge_candidates = []
                        common_candidates = None

                        for member in current:
                            neighbors = merge_neighbors.get(member, set())
                            candidate_ids = set()
                            for neighbor in neighbors:
                                candidate_ids.update(member_to_cluster_ids.get(neighbor, set()))
                            candidate_ids = {cid for cid in candidate_ids if not used[cid]}
                            if common_candidates is None:
                                common_candidates = candidate_ids
                            else:
                                common_candidates &= candidate_ids
                            if not common_candidates:
                                break

                        for j in sorted(common_candidates or []):
                            other = clusters[j]
                            if (
                                self.enforce_one_per_genome
                                and current_genomes & cluster_genome_sets[j]
                            ):
                                continue

                            ok = True
                            for member in current:
                                if not other.issubset(merge_neighbors.get(member, set())):
                                    ok = False
                                    break
                            if not ok:
                                continue

                            distances = []
                            evalues = []
                            for a in current:
                                for b in other:
                                    hit_1 = accepted_hits[(a, b)]
                                    hit_2 = accepted_hits[(b, a)]
                                    distances.append(max(hit_1["distance"], hit_2["distance"]))
                                    evalues.append(max(hit_1["e_value"], hit_2["e_value"]))
                            merge_candidates.append(((max(distances), max(evalues)), j))

                        if merge_candidates:
                            merge_candidates.sort()
                            _, best_j = merge_candidates[0]
                            current |= clusters[best_j]
                            current_genomes |= cluster_genome_sets[best_j]
                            used[best_j] = True
                            grew = True
                            changed = True
                            cluster_merges += 1

                    new_clusters.append(current)
                    new_cluster_genome_sets.append(current_genomes)

                clusters = new_clusters
                cluster_genome_sets = new_cluster_genome_sets

        log_stage("stage3_merge" if self.enable_cluster_merge else "stage3_merge_skipped", s)

        # Stage 4: enforce one protein per genome, then recover singletons
        s = perf_counter()
        cleaned_clusters = []
        for cluster in clusters:
            if not self.enforce_one_per_genome:
                cleaned_clusters.append(set(cluster))
                continue

            by_genome = {}
            for member in cluster:
                genome = protein_to_genome.get(member, f"UNKNOWN::{member}")
                by_genome.setdefault(genome, []).append(member)

            kept = set()
            for genome, members in by_genome.items():
                if len(members) == 1:
                    kept.add(members[0])
                else:
                    winner = sorted(members, key=support_score)[0]
                    kept.add(winner)
            cleaned_clusters.append(kept)

        unique_clusters = []
        seen = set()
        for cluster in cleaned_clusters:
            frozen = frozenset(cluster)
            if cluster and frozen not in seen:
                seen.add(frozen)
                unique_clusters.append(cluster)

        final_assigned = set()
        for cluster in unique_clusters:
            final_assigned.update(cluster)

        final_clusters = list(unique_clusters)
        for protein in sorted(proteins):
            if protein not in final_assigned:
                final_clusters.append({protein})

        self.head_to_members = {}
        self.member_to_cluster = {}

        multi_member_clusters = 0
        singleton_clusters = 0
        for members in final_clusters:
            canonical = sorted(members)[0]
            self.head_to_members[canonical] = set(members)
            for member in members:
                self.member_to_cluster[member] = canonical

            if len(members) > 1:
                multi_member_clusters += 1
            else:
                singleton_clusters += 1

        log_stage("stage4_finalize", s)

        print(
            "[resolve_clusters] "
            f"proteins={len(proteins)} "
            f"accepted_directed_hits={kept_rows} "
            f"rbh_seed_pairs={len(seed_pairs)} "
            f"expansions={expansions} "
            f"cluster_merges={cluster_merges} "
            f"multi_member_clusters={multi_member_clusters} "
            f"singleton_clusters={singleton_clusters} "
            f"clusters={len(self.head_to_members)} "
            f"skipped_self={skipped_self} "
            f"skipped_threshold={skipped_threshold} "
            f"skipped_same_genome={skipped_same_genome} "
            f"merge_stage={'on' if self.enable_cluster_merge else 'off'} "
            f"enforce_one_per_genome={'on' if self.enforce_one_per_genome else 'off'} "
            f"stage1_distance={self.stage1_distance} "
            f"stage2_distance={self.stage2_distance} "
            f"stage3_distance={self.stage3_distance} "
            f"total_time={perf_counter() - t0:.1f}s"
        )

    # ------------- Cluster Summary, PROST QC --------------
    def print_cluster_summary(self, max_clusters: int = 20):
        print("\n[cluster summary]")
        items = sorted(
            self.head_to_members.items(),
            key=lambda kv: len(kv[1]),
            reverse=True,
        )

        for cluster_id, members in items[:max_clusters]:
            products = []
            for member in members:
                orf = self.gff_to_prost_map.get(member)
                if not orf:
                    continue
                gene = self.orf_lookup.get(orf)
                if gene and gene.get("product"):
                    products.append(gene["product"])

            top_products = Counter(products).most_common(5)
            print(
                f"cluster={cluster_id} "
                f"size={len(members)} "
                f"top_products={top_products}"
            )

    # ---------------- Nodes ----------------
    def _choose_preferred_gene_name(self, gene_names: List[str]) -> str:
        """Choose a representative gene name, optionally preferring configured prefixes."""
        cleaned = [g for g in gene_names if g and str(g).strip()]
        if not cleaned:
            return "Unknown Gene"

        gene_counts = Counter(cleaned)
        for prefix in self.preferred_gene_prefixes:
            matches = [(g, c) for g, c in gene_counts.items() if g.lower().startswith(prefix)]
            if matches:
                matches.sort(key=lambda x: (-x[1], x[0]))
                return matches[0][0]

        return gene_counts.most_common(1)[0][0]

    def create_nodes(self):
        self.nodes = set(self.head_to_members.keys())
        self.node_attributes = {}

        for head, members in self.head_to_members.items():
            genomes = set()
            products = []
            gene_names = []

            reference_ids = []
            reference_products = []
            reference_gene_names = []

            for member in members:
                orf = self.gff_to_prost_map.get(member)
                if not orf:
                    continue

                gene = self.orf_lookup.get(orf)
                if gene is None:
                    continue

                genome = gene.get("genome")
                if genome:
                    genomes.add(genome)

                product = gene.get("product")
                if product:
                    products.append(product)

                gene_name = gene.get("gene_name")
                if gene_name:
                    gene_names.append(gene_name)

                if self.reference_genome and genome == self.reference_genome:
                    reference_ids.append(member)
                    if product:
                        reference_products.append(product)
                    if gene_name:
                        reference_gene_names.append(gene_name)

            if products:
                product_counts = Counter(products)
                sorted_products = product_counts.most_common()

                chosen_product = None
                for product, _ in sorted_products:
                    if "hypothetical protein" not in product.lower():
                        chosen_product = product
                        break
                if chosen_product is None:
                    chosen_product = sorted_products[0][0]
                most_common_product = chosen_product
            else:
                most_common_product = "Unknown Product"

            most_common_gene = self._choose_preferred_gene_name(gene_names)

            if reference_products:
                ref_product_counts = Counter(reference_products)
                sorted_ref_products = ref_product_counts.most_common()

                chosen_ref_product = None
                for product, _ in sorted_ref_products:
                    if "hypothetical protein" not in product.lower():
                        chosen_ref_product = product
                        break
                if chosen_ref_product is None:
                    chosen_ref_product = sorted_ref_products[0][0]
            else:
                chosen_ref_product = ""

            if reference_gene_names:
                chosen_ref_gene = self._choose_preferred_gene_name(reference_gene_names)
            else:
                chosen_ref_gene = ""

            all_gene_names = sorted(set(gene_names))
            all_products = sorted(set(products))
            chosen_ref_id = sorted(reference_ids)[0] if reference_ids else ""

            self.node_attributes[head] = {
                "cluster_size": len(members),
                "genomes": ",".join(sorted(genomes)),
                "product": most_common_product,
                "gene_name": most_common_gene,
                "members": ";".join(sorted(members)),
                "all_gene_names": ";".join(all_gene_names),
                "all_products": ";".join(all_products),
                "reference_id": chosen_ref_id,
                "reference_members": ";".join(sorted(reference_ids)) if reference_ids else "",
                "reference_product": chosen_ref_product,
                "reference_gene_name": chosen_ref_gene,
            }
            if self.de_by_gene:
                self.node_attributes[head].update(self._de_attributes_for_members(members))

    # ---------------- Edges ----------------
    def create_gene_pairs(self):
        """
        Create adjacent gene pairs within each operon list.

        Pairs are always emitted in 5'->3' transcription order.
        If self.adjacency_only is True: connect all adjacent genes
        regardless of gap.
        Else: connect adjacent genes only if gap <= self.operonic_distance.
        """
        self.gene_pairs = []
        skipped_gap = 0
        skipped_mismatch = 0
        skipped_empty = 0

        for operon in self.operons:
            if not operon:
                skipped_empty += 1
                continue

            contigs = {g["contig"] for g in operon}
            strands = {g["strand"] for g in operon}
            if len(contigs) != 1 or len(strands) != 1:
                skipped_mismatch += max(len(operon) - 1, 1)
                continue

            strand = next(iter(strands))
            if strand == "+":
                ordered = sorted(operon, key=lambda g: (g["start"], g["end"], g["raw_id"]))
            else:
                ordered = sorted(operon, key=lambda g: (-g["start"], -g["end"], g["raw_id"]))

            for i in range(len(ordered) - 1):
                a = ordered[i]
                b = ordered[i + 1]

                if a["contig"] != b["contig"] or a["strand"] != b["strand"]:
                    skipped_mismatch += 1
                    continue

                if strand == "+":
                    gap = b["start"] - a["end"] - 1
                else:
                    gap = a["start"] - b["end"] - 1

                if gap < 0:
                    gap = 0

                if getattr(self, "adjacency_only", False):
                    self.gene_pairs.append((a["raw_id"], b["raw_id"], a["genome"], gap))
                else:
                    if gap <= self.operonic_distance:
                        self.gene_pairs.append((a["raw_id"], b["raw_id"], a["genome"], gap))
                    else:
                        skipped_gap += 1

        print(
            f"[create_gene_pairs] pairs={len(self.gene_pairs)} "
            f"skipped_by_gap={skipped_gap} "
            f"skipped_mismatch={skipped_mismatch} "
            f"skipped_empty={skipped_empty}"
        )

    def create_edges(self):
        self.edges = []
        skipped_missing = 0
        self_loops = 0

        for g1, g2, genome, gap in self.gene_pairs:
            c1 = self.member_to_cluster.get(g1)
            c2 = self.member_to_cluster.get(g2)

            if c1 is None or c2 is None:
                skipped_missing += 1
                continue

            if c1 == c2:
                self_loops += 1
                if self_loops <= 20:
                    print(f"[self_loop] {g1} {g2} -> {c1}")

            self.edges.append((c1, c2, genome, g1, g2, gap))

        print(
            f"[create_edges] edges={len(self.edges)} "
            f"skipped_missing={skipped_missing} self_loops={self_loops}"
        )

    # ---------------- Weights ----------------
    def calculate_edge_properties(self):
        for u, v, genome, g1, g2, gap in self.edges:
            n1 = len(self.head_to_members[u])
            n2 = len(self.head_to_members[v])
            weight = 1 / (n1 * n2)

            edge_key = (u, v)
            if edge_key not in self.edge_properties:
                self.edge_properties[edge_key] = {
                    "abs_weight": 0,
                    "norm_weight": 0,
                    "genome_counts": {},
                    "member_gene_pairs": [],
                    "intergenic_distances": [],
                }

            self.edge_properties[edge_key]["abs_weight"] += 1
            self.edge_properties[edge_key]["norm_weight"] += weight

            genome_counts = self.edge_properties[edge_key]["genome_counts"]
            genome_counts[genome] = genome_counts.get(genome, 0) + 1
            self.edge_properties[edge_key]["member_gene_pairs"].append((genome, g1, g2))
            self.edge_properties[edge_key]["intergenic_distances"].append(int(gap))

    def make_weighted_edges(self):
        self.weighted_edges = []
        for (u, v), data in self.edge_properties.items():
            genome_counts_str = ";".join(
                f"{genome}:{count}" for genome, count in data["genome_counts"].items()
            )
            member_gene_pairs_str = ";".join(
                f"{genome}|{g1}|{g2}"
                for genome, g1, g2 in sorted(data["member_gene_pairs"])
            )
            distances = list(data["intergenic_distances"])
            intergenic_distances_str = ";".join(str(distance) for distance in distances)
            avg_intergenic_distance = (
                sum(distances) / len(distances) if distances else 0.0
            )
            edge_data = {
                "abs_weight": data["abs_weight"],
                "norm_weight": data["norm_weight"],
                "genome_counts": genome_counts_str,
                "num_genomes": len(data["genome_counts"]),
                "member_gene_pairs": member_gene_pairs_str,
                "intergenic_distances": intergenic_distances_str,
                "avg_intergenic_distance": avg_intergenic_distance,
            }
            self.weighted_edges.append((u, v, edge_data))

    # ---------------- Run ----------------
    def run(self):
        self.load_operons()
        self.resolve_clusters()
        self.print_cluster_summary(max_clusters=20)
        self.create_nodes()
        self.create_gene_pairs()
        self.create_edges()
        self.calculate_edge_properties()
        self.make_weighted_edges()

        self.G.add_nodes_from(self.nodes)
        nx.set_node_attributes(self.G, self.node_attributes)
        self.G.add_edges_from(self.weighted_edges)
        return self.G

    # ---------------- Export helpers ----------------
    def export_gml(self, out_gml: str):
        """Write the built graph to a .gml file for Cytoscape Desktop."""
        if self.G is None or self.G.number_of_nodes() == 0:
            raise ValueError("Graph is empty. Run builder.run() before exporting.")
        nx.write_gml(self.G, out_gml)
        return out_gml

    def export_tsv(self, nodes_tsv: str, edges_tsv: str):
        """Export node and edge tables as TSVs for downstream use."""
        if self.G is None or self.G.number_of_nodes() == 0:
            raise ValueError("Graph is empty. Run builder.run() before exporting.")

        node_rows = []
        for node in self.G.nodes():
            attrs = self.G.nodes[node]
            node_rows.append(
                {
                    "id": node,
                    "cluster_size": attrs.get("cluster_size", None),
                    "genomes": attrs.get("genomes", None),
                    "product": attrs.get("product", None),
                    "gene_name": attrs.get("gene_name", None),
                    "all_gene_names": attrs.get("all_gene_names", None),
                    "node_log2FC": attrs.get("node_log2FC", None),
                    "node_min_padj": attrs.get("node_min_padj", None),
                    "node_num_DE": attrs.get("node_num_DE", None),
                    "node_frac_DE": attrs.get("node_frac_DE", None),
                    "node_mean_expr": attrs.get("node_mean_expr", None),
                    "node_mean_expr_gene_count": attrs.get("node_mean_expr_gene_count", None),
                }
            )
        pd.DataFrame(node_rows).to_csv(nodes_tsv, sep="\t", index=False)

        edge_rows = []
        for u, v, attrs in self.G.edges(data=True):
            edge_rows.append(
                {
                    "source": u,
                    "target": v,
                    "abs_weight": attrs.get("abs_weight", None),
                    "norm_weight": attrs.get("norm_weight", None),
                    "member_gene_pairs": attrs.get("member_gene_pairs", None),
                    "intergenic_distances": attrs.get("intergenic_distances", None),
                    "avg_intergenic_distance": attrs.get("avg_intergenic_distance", None),
                }
            )
        pd.DataFrame(edge_rows).to_csv(edges_tsv, sep="\t", index=False)
        return nodes_tsv, edges_tsv


def main():
    import argparse
    from pathlib import Path

    ap = argparse.ArgumentParser(
        description="Build OperonGraph adjacency network from PROST TSV + Prokka GFF folders"
    )
    ap.add_argument("--prost_tsv", required=True, help="PROST results TSV (combined.tsv)")
    ap.add_argument(
        "--gff_dir",
        required=True,
        help="Directory containing genome subfolders each with a .gff",
    )
    ap.add_argument("--out_gml", required=True, help="Output .gml path")
    ap.add_argument(
        "--operonic_distance",
        type=int,
        default=40,
        help="Max bp distance to connect adjacent genes when not using --adjacency_only",
    )
    ap.add_argument("--out_nodes_tsv", default="", help="Optional nodes TSV path")
    ap.add_argument("--out_edges_tsv", default="", help="Optional edges TSV path")
    ap.add_argument("--force", action="store_true", help="Overwrite outputs if they exist")
    ap.add_argument(
        "--adjacency_only",
        action="store_true",
        help="If set, connect adjacent genes regardless of bp gap.",
    )
    ap.add_argument(
        "--cluster_distance",
        type=float,
        default=5000.0,
        help=(
            "Legacy shared PROST distance cutoff used for stages 1, 2, and 3 "
            "unless stage-specific distances are provided."
        ),
    )
    ap.add_argument(
        "--stage1_distance",
        type=float,
        default=None,
        help="Maximum PROST distance for stage 1 RBH seed selection.",
    )
    ap.add_argument(
        "--stage2_distance",
        type=float,
        default=None,
        help="Maximum PROST distance for stage 2 singleton expansion and vacuuming.",
    )
    ap.add_argument(
        "--stage3_distance",
        type=float,
        default=None,
        help="Maximum PROST distance for stage 3 cluster-cluster consolidation.",
    )
    ap.add_argument(
        "--cluster_evalue",
        type=float,
        default=1e-10,
        help="Maximum e-value to merge proteins into the same cluster.",
    )
    ap.add_argument(
        "--no_reciprocal_only",
        action="store_true",
        help="If set, cluster on any thresholded hit rather than requiring reciprocal hits.",
    )
    ap.add_argument(
        "--reference_genome",
        default=None,
        help="Optional genome name to use for node reference annotations.",
    )
    ap.add_argument(
        "--prefer_gene_prefixes",
        default="",
        help="Comma-separated gene prefixes to preferentially use for node naming.",
    )
    ap.add_argument(
        "--de_csv",
        default="",
        help="Optional differential-expression CSV to aggregate onto graph nodes.",
    )
    ap.add_argument(
        "--de_gene_col",
        default="Geneid",
        help="DE CSV gene-id column. Default: Geneid.",
    )
    ap.add_argument(
        "--de_log2fc_col",
        default="",
        help="DE CSV log2 fold-change column. Auto-detected if omitted.",
    )
    ap.add_argument(
        "--de_padj_col",
        default="",
        help="DE CSV adjusted p-value column. Auto-detected if omitted.",
    )
    ap.add_argument(
        "--de_expr_col",
        default="",
        help="DE CSV expression column for node_mean_expr. Auto-detected if omitted.",
    )
    ap.add_argument(
        "--de_padj_threshold",
        type=float,
        default=0.05,
        help="Adjusted p-value threshold for node_num_DE/node_frac_DE. Default: 0.05.",
    )
    ap.add_argument(
        "--de_abs_log2fc_threshold",
        type=float,
        default=1.0,
        help="Absolute log2FC threshold for node_num_DE/node_frac_DE. Default: 1.0.",
    )
    ap.add_argument(
        "--de_expression_agg",
        choices=["mean", "sum"],
        default="mean",
        help="Aggregate expression across matched node genes as mean or sum. Default: mean.",
    )
    ap.add_argument(
        "--enable_cluster_merge",
        action="store_true",
        help="Enable the expensive final cluster-cluster merge pass. Off by default for large runs.",
    )
    ap.add_argument(
        "--enforce_one_per_genome",
        dest="enforce_one_per_genome",
        action="store_true",
        help=(
            "Keep at most one protein per genome in each cluster. "
            "This is the default behavior."
        ),
    )
    ap.add_argument(
        "--no_enforce_one_per_genome",
        dest="enforce_one_per_genome",
        action="store_false",
        help=(
            "Allow multiple proteins from the same genome to remain in the same cluster. "
            "Default behavior keeps at most one protein per genome per cluster."
        ),
    )
    ap.set_defaults(enforce_one_per_genome=True)
    ap.add_argument(
        "--progress_every",
        type=int,
        default=1_000_000,
        help="Accepted PROST hits between progress messages during stage 0 filtering.",
    )

    args = ap.parse_args()

    out_gml = Path(args.out_gml)
    if out_gml.exists() and not args.force:
        raise SystemExit(f"ERROR: {out_gml} exists. Use --force to overwrite.")

    if args.out_nodes_tsv and not args.out_edges_tsv:
        raise SystemExit("ERROR: If you set --out_nodes_tsv you must also set --out_edges_tsv.")
    if args.out_edges_tsv and not args.out_nodes_tsv:
        raise SystemExit("ERROR: If you set --out_edges_tsv you must also set --out_nodes_tsv.")

    builder = GraphBuilderProst(
        prost_tsv=args.prost_tsv,
        gff_dir=args.gff_dir,
        operonic_distance=args.operonic_distance,
        adjacency_only=args.adjacency_only,
        cluster_distance=args.cluster_distance,
        stage1_distance=args.stage1_distance,
        stage2_distance=args.stage2_distance,
        stage3_distance=args.stage3_distance,
        cluster_evalue=args.cluster_evalue,
        reciprocal_only=not args.no_reciprocal_only,
        reference_genome=args.reference_genome,
        preferred_gene_prefixes=[x.strip() for x in args.prefer_gene_prefixes.split(",") if x.strip()],
        de_csv=args.de_csv or None,
        de_gene_col=args.de_gene_col,
        de_log2fc_col=args.de_log2fc_col or None,
        de_padj_col=args.de_padj_col or None,
        de_expr_col=args.de_expr_col or None,
        de_padj_threshold=args.de_padj_threshold,
        de_abs_log2fc_threshold=args.de_abs_log2fc_threshold,
        de_expression_agg=args.de_expression_agg,
        enable_cluster_merge=args.enable_cluster_merge,
        enforce_one_per_genome=args.enforce_one_per_genome,
        progress_every=args.progress_every,
    )

    print(
        "[mode] adjacency_rule="
        + ("adjacent_same_strand_any_gap" if args.adjacency_only else f"adjacent_same_strand_gap_le_{args.operonic_distance}bp")
    )
    print(
        "[mode] clustering "
        f"stage1_distance={builder.stage1_distance} "
        f"stage2_distance={builder.stage2_distance} "
        f"stage3_distance={builder.stage3_distance} "
        f"enforce_one_per_genome={builder.enforce_one_per_genome}"
    )

    graph = builder.run()
    print(f"[done] nodes={graph.number_of_nodes()} edges={graph.number_of_edges()}")

    builder.export_gml(str(out_gml))
    print(f"[wrote] {out_gml}")

    if args.out_nodes_tsv and args.out_edges_tsv:
        if Path(args.out_nodes_tsv).exists() and not args.force:
            raise SystemExit(f"ERROR: {args.out_nodes_tsv} exists. Use --force to overwrite.")
        if Path(args.out_edges_tsv).exists() and not args.force:
            raise SystemExit(f"ERROR: {args.out_edges_tsv} exists. Use --force to overwrite.")
        builder.export_tsv(args.out_nodes_tsv, args.out_edges_tsv)
        print(f"[wrote] {args.out_nodes_tsv}")
        print(f"[wrote] {args.out_edges_tsv}")


if __name__ == "__main__":
    main()
