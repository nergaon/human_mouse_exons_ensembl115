#!/usr/bin/env python3
"""Merge Leafcutter analysis results into summary tables matching clusters_sum_table_HN6.txt format.

Output: one row per junction from AS_clusters_value_HN6.txt.
Columns: h_junction, m_junction, symbol_h, ensembl_h, rank_h, rank_m, cluster, genes,
         {cell_type}_abs_deltapsi, {cell_type}_p.adjust,
         {dataset_a}_avg_{cell_type}, {dataset_b}_avg_{cell_type}, ...
"""

import csv
import re
from pathlib import Path
from typing import Dict, List, Tuple

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None

SHARED_JUNCTIONS_DIR = Path("/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026")
METADATA_FILE = Path("/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/unique_points_HN6.txt")
VALUE_FILE = SHARED_JUNCTIONS_DIR / "AS_clusters_value_HN6.txt"


def load_metadata(metadata_file: Path) -> Dict[str, Dict]:
    """Load unique_points_HN6.txt and index by position_h for coordinate lookup."""
    metadata = {}
    with metadata_file.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pos_h = row.get("position_h", "").strip()
            if pos_h:
                metadata[pos_h] = row
    return metadata


def parse_dataset(sample_name: str) -> str:
    match = re.search(r"_(GSE\d+)$", sample_name)
    return match.group(1) if match else ""


def infer_celltype_candidates(sample_name: str) -> List[str]:
    base = re.sub(r"_GSE\d+$", "", sample_name)

    if base.startswith("NveB"):
        return ["NveB", "BCell"]
    if base.startswith(("BCell", "Bcells", "B_fo")):
        return ["BCell", "NveB"]

    if base.startswith(("CD4T", "EffCD4T", "MemCD4T", "NveCD4T", "CD4_", "T_4_")):
        return ["CD4T"]

    if base.startswith(("CD8T", "MemCD8T", "NveCd8T", "CD8_", "T_8_")):
        return ["CD8T", "Cd8T"]

    if base.startswith(("Mono", "InfMono", "Monocytes", "Mo_6C+")):
        return ["Mono"]

    if base.startswith(("Neut", "Neutrophils")):
        return ["Neut"]

    if base.startswith("NK"):
        return ["NK"]

    return []


def load_junctions_with_as_averages(
        value_file: Path,
        dataset_a: str,
        dataset_b: str,
        cell_types: List[str],
) -> List[Dict]:
    """Load junction rows and compute avg expression per (dataset, cell_type) from AS value file."""
    junctions = []

    with value_file.open("r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header:
            return junctions

        # Build column mapping once: AS column idx -> (dataset, cell_type)
        sample_mappings: List[Tuple[int, str, str]] = []
        valid_cell_types = set(cell_types)
        for col_idx, sample_name in enumerate(header[1:], start=1):
            dataset = parse_dataset(sample_name)
            if dataset not in {dataset_a, dataset_b}:
                continue

            candidates = infer_celltype_candidates(sample_name)
            mapped_cell = next((ct for ct in candidates if ct in valid_cell_types), "")
            if mapped_cell:
                sample_mappings.append((col_idx, dataset, mapped_cell))

        for row in reader:
            if not row:
                continue
            junction_id = row[0].strip()
            if not junction_id or ":" not in junction_id:
                continue

            parts = junction_id.rsplit(":", 1)
            if len(parts) != 2:
                continue

            junction = parts[0]       # chr:start:end
            cluster_suffix = parts[1] # clu_XXXXX
            chrom = junction.split(":")[0]
            cluster_id = f"{chrom}:{cluster_suffix}"

            sums: Dict[Tuple[str, str], float] = {}
            counts: Dict[Tuple[str, str], int] = {}
            for col_idx, dataset, cell_type in sample_mappings:
                if col_idx >= len(row):
                    continue
                try:
                    value = float(row[col_idx])
                except (ValueError, TypeError):
                    continue

                key = (dataset, cell_type)
                sums[key] = sums.get(key, 0.0) + value
                counts[key] = counts.get(key, 0) + 1

            avgs: Dict[Tuple[str, str], str] = {}
            for key, total in sums.items():
                n = counts.get(key, 0)
                avgs[key] = f"{(total / n):.10g}" if n > 0 else ""

            junctions.append({
                "h_junction": junction,
                "cluster_id": cluster_id,
                "as_avgs": avgs,
            })

    return junctions


def lookup_junction_metadata(junction: str, metadata: Dict[str, Dict]) -> Dict[str, str]:
    """Build metadata for a junction by looking up BOTH start and end positions.

    m_junction = chr_m:start_m:end_m
    rank_h     = startRank_endRank  (e.g. 'E6_E5')
    rank_m     = startRank_m_endRank_m
    """
    result = {"m_junction": "", "ensembl_h": "", "symbol_h": "", "rank_h": "", "rank_m": ""}

    parts = junction.split(":")
    if len(parts) < 3:
        return result

    chrom, start_pos, end_pos = parts[0], parts[1], parts[2]
    start_key = f"{chrom}:{start_pos}"
    end_key   = f"{chrom}:{end_pos}"

    start_meta = metadata.get(start_key, {})
    end_meta   = metadata.get(end_key, {})

    # --- m_junction: combine mouse start and mouse end ---
    m_start = start_meta.get("position_m", "")
    m_end   = end_meta.get("position_m", "")
    if m_start and m_end:
        m_start_parts = m_start.split(":")
        m_end_parts   = m_end.split(":")
        if len(m_start_parts) >= 2 and len(m_end_parts) >= 2:
            m_chrom = m_start_parts[0]
            result["m_junction"] = f"{m_chrom}:{m_start_parts[1]}:{m_end_parts[1]}"
    elif m_start:
        result["m_junction"] = m_start
    elif m_end:
        result["m_junction"] = m_end

    # --- rank_h: start_rank_end_rank ---
    s_rank_h = start_meta.get("rank_h", "")
    e_rank_h = end_meta.get("rank_h", "")
    result["rank_h"] = f"{s_rank_h}_{e_rank_h}" if s_rank_h and e_rank_h else s_rank_h or e_rank_h

    # --- rank_m: start_rank_m_end_rank_m ---
    s_rank_m = start_meta.get("rank_m", "")
    e_rank_m = end_meta.get("rank_m", "")
    result["rank_m"] = f"{s_rank_m}_{e_rank_m}" if s_rank_m and e_rank_m else s_rank_m or e_rank_m

    # --- gene info: prefer start, fall back to end ---
    result["ensembl_h"] = start_meta.get("ensembl", "") or end_meta.get("ensembl", "")
    result["symbol_h"]  = start_meta.get("symbol",  "") or end_meta.get("symbol",  "")

    return result


def load_cluster_significance(sig_file: Path) -> Dict[str, Dict]:
    """Load leafcutter_ds_cluster_significance.txt by cluster ID."""
    clusters = {}
    with sig_file.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cluster_id = row.get("cluster", "").strip()
            if cluster_id:
                clusters[cluster_id] = row
    return clusters


def load_effect_sizes_by_junction(effect_file: Path) -> Dict[str, Dict]:
    """Load leafcutter_ds_effect_sizes.txt indexed by junction coordinates (chr:start:end).

    The intron column format is 'chr:start:end:clu_XXXXX'; we strip the cluster suffix
    so the key is 'chr:start:end' for per-junction lookup.
    """
    effects = {}
    with effect_file.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            intron = row.get("intron", "").strip()
            if not intron or ":" not in intron:
                continue
            parts = intron.rsplit(":", 1)
            if len(parts) == 2:
                junction_key = parts[0]   # chr:start:end
                effects[junction_key] = row
    return effects


def _safe_float(value: str) -> float:
    try:
        return float(value)
    except (ValueError, TypeError):
        return float("nan")


def _cluster_id_from_intron(intron: str) -> str:
    """Convert intron format chr:start:end:clu_X -> chr:clu_X."""
    if not intron or ":" not in intron:
        return ""
    parts = intron.rsplit(":", 1)
    if len(parts) != 2:
        return ""
    chrom = parts[0].split(":", 1)[0]
    return f"{chrom}:{parts[1]}"


def summarize_leafcutter_celltype(cell_dir: Path) -> Dict[str, int]:
    """Summarize Leafcutter counts for one cell type."""
    sig_file = cell_dir / "leafcutter_ds_cluster_significance.txt"
    effect_file = cell_dir / "leafcutter_ds_effect_sizes.txt"

    if not sig_file.exists() or not effect_file.exists():
        return {
            "success": 0,
            "sig_p05": 0,
            "sig_dpsi_02": 0,
            "sig_dpsi_01": 0,
            "sig_dpsi_005": 0,
        }

    success_clusters = set()
    sig_clusters = set()

    with sig_file.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cluster = row.get("cluster", "").strip()
            if not cluster:
                continue

            status = row.get("status", "").strip()
            if status == "Success":
                success_clusters.add(cluster)

                p_adjust = _safe_float(row.get("p.adjust", ""))
                if p_adjust == p_adjust and p_adjust <= 0.05:
                    sig_clusters.add(cluster)

    # Cluster-level max abs(deltapsi)
    cluster_max_abs_dpsi: Dict[str, float] = {}
    with effect_file.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            intron = row.get("intron", "").strip()
            cluster = _cluster_id_from_intron(intron)
            if not cluster:
                continue

            dpsi = _safe_float(row.get("deltapsi", ""))
            if dpsi != dpsi:
                continue
            abs_dpsi = abs(dpsi)

            prev = cluster_max_abs_dpsi.get(cluster)
            if prev is None or abs_dpsi > prev:
                cluster_max_abs_dpsi[cluster] = abs_dpsi

    sig_dpsi_02 = sum(1 for c in sig_clusters if cluster_max_abs_dpsi.get(c, 0.0) >= 0.2)
    sig_dpsi_01 = sum(1 for c in sig_clusters if cluster_max_abs_dpsi.get(c, 0.0) >= 0.1)
    sig_dpsi_005 = sum(1 for c in sig_clusters if cluster_max_abs_dpsi.get(c, 0.0) >= 0.05)

    return {
        "success": len(success_clusters),
        "sig_p05": len(sig_clusters),
        "sig_dpsi_02": sig_dpsi_02,
        "sig_dpsi_01": sig_dpsi_01,
        "sig_dpsi_005": sig_dpsi_005,
    }


def write_sum_table_and_stacked_bar(leafcutter_dir: Path, dataset_a: str, dataset_b: str) -> None:
    """Write sum_table_HN6.txt (last 5 rows) and a stacked success/sig bar plot."""
    cell_dirs = sorted([
        d for d in leafcutter_dir.iterdir()
        if d.is_dir()
        and d.name != "genes_figs"
        and (d / "leafcutter_ds_cluster_significance.txt").exists()
        and (d / "leafcutter_ds_effect_sizes.txt").exists()
    ])

    cell_types: List[str] = []
    summary_by_cell: Dict[str, Dict[str, int]] = {}
    for cell_dir in cell_dirs:
        cell_type = cell_dir.name
        cell_types.append(cell_type)
        summary_by_cell[cell_type] = summarize_leafcutter_celltype(cell_dir)

    # Prefer dataset-specific order for readability in outputs.
    dataset_pair = {dataset_a, dataset_b}
    if dataset_pair == {"GSE116177", "GSE180020"}:
        preferred_order = ["CD4T", "Cd8T", "BCell", "Mono"]
    else:
        preferred_order = ["CD4T", "CD8T", "NveB", "NK", "Mono", "Neut"]

    display_cell_types = [ct for ct in preferred_order if ct in summary_by_cell] + [
        ct for ct in cell_types if ct not in preferred_order
    ]

    # Write table matching the "last 5 rows" format
    sum_table = leafcutter_dir / "sum_table_HN6.txt"
    with sum_table.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["", "All"] + display_cell_types)

        rows = [
            ("Leafcutter success clusters", "success"),
            ("Leafcutter sig clusters (p<0.05)", "sig_p05"),
            ("Leafcutter sig and deltapsi>=0.2", "sig_dpsi_02"),
            ("Leafcutter sig and deltapsi>=0.1", "sig_dpsi_01"),
            ("Leafcutter sig and deltapsi>=0.05", "sig_dpsi_005"),
        ]

        for label, key in rows:
            row = [label, ""]
            for ct in display_cell_types:
                row.append(f"{float(summary_by_cell[ct][key]):.1f}")
            writer.writerow(row)

    print(f"Wrote {sum_table}")

    if plt is None:
        print(
            "Skipping stacked bar plot: matplotlib is not installed in the active environment. "
            "Install it to enable plot output."
        )
        return

    # Stacked bar: not-sig (bottom), sig-but-no-dpsi02 (middle), sig+dpsi>=0.2 (top).
    # Use exact count decomposition:
    #   bottom = success - sig_p05
    #   middle = sig_p05 - sig_dpsi_02
    #   top    = sig_dpsi_02
    success_vals      = [summary_by_cell[ct]["success"] for ct in display_cell_types]
    sig_p05_vals      = [summary_by_cell[ct]["sig_p05"] for ct in display_cell_types]
    sig02_vals_raw    = [summary_by_cell[ct]["sig_dpsi_02"] for ct in display_cell_types]
    not_sig_vals      = [s - p for s, p in zip(success_vals, sig_p05_vals)]
    sig_not_dpsi_vals = [p - d for p, d in zip(sig_p05_vals, sig02_vals_raw)]
    sig02_vals        = sig02_vals_raw
    bottom2           = [a + b for a, b in zip(not_sig_vals, sig_not_dpsi_vals)]

    fig, ax = plt.subplots(figsize=(8, 5))
    x = list(range(len(display_cell_types)))
    ax.bar(x, not_sig_vals,      label="Success but not sig (p>0.05)", color="#deebf7")
    ax.bar(x, sig_not_dpsi_vals, bottom=not_sig_vals, label="sig (p<0.05) but abs(deltapsi)<0.2", color="#9ecae1")
    ax.bar(x, sig02_vals,        bottom=bottom2,      label="sig (p<0.05) + abs(deltapsi)>=0.2",  color="#08519c")

    ax.set_xticks(x)
    ax.set_xticklabels(display_cell_types, rotation=0)
    ax.set_ylabel("Number of clusters")
    ax.set_title(f"{dataset_a} vs {dataset_b}: Leafcutter success and significant clusters")
    ax.legend(fontsize=8)

    for i, total in enumerate(success_vals):
        ax.text(i, total + 0.3, str(int(total)), ha="center", va="bottom", fontsize=8)

    plt.tight_layout()
    plot_png = leafcutter_dir / f"{dataset_a}_{dataset_b}_stacked_success_sig_deltapsi02.png"
    plot_svg = leafcutter_dir / f"{dataset_a}_{dataset_b}_stacked_success_sig_deltapsi02.svg"
    fig.savefig(plot_png, dpi=200)
    fig.savefig(plot_svg)
    plt.close(fig)

    print(f"Wrote {plot_png}")
    print(f"Wrote {plot_svg}")


def merge_leafcutter_results(
        leafcutter_dir: Path,
        dataset_a: str,
        dataset_b: str,
        metadata: Dict[str, Dict],
) -> None:
    """Write clusters_sum_table_HN6.txt: one row per junction, per-cell-type stats as columns."""

    # Load data for each cell type
    cell_types: List[str] = []
    cell_sig: Dict[str, Dict[str, Dict]] = {}   # cell_type -> cluster_id -> row
    cell_eff: Dict[str, Dict[str, Dict]] = {}   # cell_type -> junction_key -> row

    for cell_type_dir in sorted(d for d in leafcutter_dir.iterdir()
                                if d.is_dir() and d.name != "genes_figs"):
        sig_file    = cell_type_dir / "leafcutter_ds_cluster_significance.txt"
        effect_file = cell_type_dir / "leafcutter_ds_effect_sizes.txt"
        if not sig_file.exists() or not effect_file.exists():
            continue
        ct = cell_type_dir.name
        cell_types.append(ct)
        cell_sig[ct] = load_cluster_significance(sig_file)
        cell_eff[ct] = load_effect_sizes_by_junction(effect_file)

    all_junctions = load_junctions_with_as_averages(VALUE_FILE, dataset_a, dataset_b, cell_types)

    # Build output: one row per junction
    output_rows = []
    for junc_info in all_junctions:
        h_junction = junc_info["h_junction"]
        cluster_id = junc_info["cluster_id"]

        junc_meta = lookup_junction_metadata(h_junction, metadata)

        row = {
            "h_junction": h_junction,
            "m_junction": junc_meta["m_junction"],
            "symbol_h":   junc_meta["symbol_h"],
            "ensembl_h":  junc_meta["ensembl_h"],
            "rank_h":     junc_meta["rank_h"],
            "rank_m":     junc_meta["rank_m"],
            "cluster":    cluster_id,
            "genes":      "",
        }

        for ct in cell_types:
            eff_row = cell_eff[ct].get(h_junction, {})
            if eff_row:
                # Recover cluster ID from the full intron field (chr:start:end:clu_XXXXX)
                intron_full = eff_row.get("intron", "")
                intron_parts = intron_full.rsplit(":", 1)
                ct_cluster = ""
                if len(intron_parts) == 2:
                    junc_chrom = intron_parts[0].split(":")[0]
                    ct_cluster = f"{junc_chrom}:{intron_parts[1]}"

                sig_row = cell_sig[ct].get(ct_cluster, {}) if ct_cluster else {}

                try:
                    deltapsi_abs = f"{abs(float(eff_row.get('deltapsi', 0))):.10g}"
                except (ValueError, TypeError):
                    deltapsi_abs = ""

                p_adjust    = sig_row.get("p.adjust", "")
                gene_name   = sig_row.get("genes", "")
                if gene_name and not row["genes"]:
                    row["genes"] = gene_name
                    if not row["symbol_h"]:
                        row["symbol_h"] = gene_name
            else:
                deltapsi_abs = ""
                p_adjust     = ""

            as_avgs = junc_info.get("as_avgs", {})
            dataset_a_v = as_avgs.get((dataset_a, ct), "")
            dataset_b_v = as_avgs.get((dataset_b, ct), "")

            row[f"{ct}_abs_deltapsi"]      = deltapsi_abs
            row[f"{ct}_p.adjust"]          = p_adjust
            row[f"{dataset_a}_avg_{ct}"]   = dataset_a_v
            row[f"{dataset_b}_avg_{ct}"]   = dataset_b_v

        # Last resort: fill genes from significance using cluster from value file
        if not row["genes"]:
            for ct in cell_types:
                sig_row = cell_sig[ct].get(cluster_id, {})
                if sig_row.get("genes"):
                    row["genes"] = sig_row["genes"]
                    if not row["symbol_h"]:
                        row["symbol_h"] = sig_row["genes"]
                    break

        output_rows.append(row)

    # Write TSV
    output_file = leafcutter_dir / "clusters_sum_table_HN6.txt"
    base_cols = ["h_junction", "m_junction", "symbol_h", "ensembl_h",
                 "rank_h", "rank_m", "cluster", "genes"]
    stat_cols = []
    for ct in cell_types:
        stat_cols += [
            f"{ct}_abs_deltapsi",
            f"{ct}_p.adjust",
            f"{dataset_a}_avg_{ct}",
            f"{dataset_b}_avg_{ct}",
        ]
    final_cols = base_cols + stat_cols

    with output_file.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=final_cols, delimiter="\t", lineterminator="\n")
        writer.writerow({c: c for c in final_cols})
        for row in output_rows:
            writer.writerow({c: row.get(c, "") for c in final_cols})

    print(f"Wrote {output_file}: {len(output_rows)} rows, {len(cell_types)} cell types")


def main() -> None:
    metadata      = load_metadata(METADATA_FILE)

    leafcutter_dirs = sorted(d for d in SHARED_JUNCTIONS_DIR.iterdir()
                             if d.is_dir() and d.name.startswith("leafcutter_GSE"))

    for leaf_dir in leafcutter_dirs:
        match = re.match(r"leafcutter_(GSE\d+)_(GSE\d+)$", leaf_dir.name)
        if not match:
            continue
        dataset_a, dataset_b = match.group(1), match.group(2)
        print(f"Processing {leaf_dir.name} ({dataset_a} vs {dataset_b})...")
        merge_leafcutter_results(leaf_dir, dataset_a, dataset_b, metadata)
        write_sum_table_and_stacked_bar(leaf_dir, dataset_a, dataset_b)


if __name__ == "__main__":
    main()
