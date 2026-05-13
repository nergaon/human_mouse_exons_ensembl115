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
from typing import Dict, List, Optional, Tuple

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None

SHARED_JUNCTIONS_DIR = Path("/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026")
METADATA_FILE = Path("/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/unique_points_HN6.txt")
VALUE_FILE = SHARED_JUNCTIONS_DIR / "AS_clusters_value_HN6.txt"
FIBRO_VALUE_FILE = SHARED_JUNCTIONS_DIR / "AS_clusters_value_fibroblast_HN6.txt"

COMPARISON_LABELS: Dict[Tuple[str, str], str] = {
    ("GSE115736", "GSE116177"): "Human vs Mouse",
    ("GSE115736", "GSE60424"): "Human vs Human",
    ("GSE116177", "GSE180020"): "Mouse vs Mouse",
}

POSSIBLE_AS_CLUSTERS = 1823


def _normalize_pair(dataset_a: str, dataset_b: str) -> Tuple[str, str]:
    return tuple(sorted((dataset_a, dataset_b)))


def _comparison_title(dataset_a: str, dataset_b: str) -> str:
    key = _normalize_pair(dataset_a, dataset_b)
    return COMPARISON_LABELS.get(key, f"{dataset_a} vs {dataset_b}")


def _display_cell_name(cell_type: str, dataset_a: str, dataset_b: str) -> str:
    """Normalize display names for specific comparisons."""
    pair = _normalize_pair(dataset_a, dataset_b)
    if cell_type == "Cd8T":
        return "CD8T"
    if pair == _normalize_pair("GSE116177", "GSE180020") and cell_type == "BCell":
        return "NveB"
    return cell_type


def _value_file_for_pair(dataset_a: str, dataset_b: str) -> Path:
    if _normalize_pair(dataset_a, dataset_b) == _normalize_pair("EMTAB5919H", "EMTAB5919M"):
        return FIBRO_VALUE_FILE
    return VALUE_FILE


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
    is_fibro_pair = _normalize_pair(dataset_a, dataset_b) == _normalize_pair("EMTAB5919H", "EMTAB5919M")

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
            candidates = infer_celltype_candidates(sample_name)
            mapped_cell = next((ct for ct in candidates if ct in valid_cell_types), "")

            if dataset in {dataset_a, dataset_b} and mapped_cell:
                sample_mappings.append((col_idx, dataset, mapped_cell))
                continue

            # Fibroblast merged matrix uses HS*/MM* sample names, not *_GSE* names.
            if is_fibro_pair and "Fibroblast" in valid_cell_types:
                upper_name = sample_name.upper()
                if upper_name.startswith("HS"):
                    sample_mappings.append((col_idx, "EMTAB5919H", "Fibroblast"))
                elif upper_name.startswith("MM"):
                    sample_mappings.append((col_idx, "EMTAB5919M", "Fibroblast"))


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
            "unchanged_dpsi_lt_005": 0,
            "not_informative": 0,
        }

    success_clusters = set()
    sig_clusters = set()
    p_adjust_by_cluster: Dict[str, float] = {}

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
                p_adjust_by_cluster[cluster] = p_adjust
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
    unchanged_dpsi_lt_005 = sum(
        1
        for c in success_clusters
        if (cluster_max_abs_dpsi.get(c, 0.0) < 0.05) or (p_adjust_by_cluster.get(c, float("nan")) > 0.05)
    )
    not_informative = max(0, len(success_clusters) - sig_dpsi_01 - unchanged_dpsi_lt_005)

    return {
        "success": len(success_clusters),
        "sig_p05": len(sig_clusters),
        "sig_dpsi_02": sig_dpsi_02,
        "sig_dpsi_01": sig_dpsi_01,
        "sig_dpsi_005": sig_dpsi_005,
        "unchanged_dpsi_lt_005": unchanged_dpsi_lt_005,
        "not_informative": not_informative,
    }


def load_summary_from_sum_table(sum_table_file: Path) -> Optional[Dict[str, int]]:
    """Load summary counts from an existing sum_table_HN6.txt file."""
    if not sum_table_file.exists():
        return None

    key_map = {
        "Leafcutter success clusters": "success",
        "Leafcutter sig clusters (p<=0.05, deltapsi>=0.1)": "sig_dpsi_01",
        "Leafcutter sig clusters (p≤0.05 and deltapsi≥0.1)": "sig_dpsi_01",
        "Unchanged (deltapsi<0.05 or p>0.05)": "unchanged_dpsi_lt_005",
        "Not informative": "not_informative",
    }
    out: Dict[str, int] = {}

    with sum_table_file.open("r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header:
            return None

        # Prefer Fibroblast column when present; otherwise use All.
        value_idx = 1 if len(header) > 1 else None
        if "Fibroblast" in header:
            value_idx = header.index("Fibroblast")
        elif "All" in header:
            value_idx = header.index("All")

        if value_idx is None:
            return None

        for row in reader:
            if not row:
                continue
            label = row[0]
            target = key_map.get(label)
            if not target:
                continue
            raw = row[value_idx].strip() if value_idx < len(row) else ""
            try:
                out[target] = int(round(float(raw))) if raw else 0
            except ValueError:
                out[target] = 0

    if not out:
        return None
    return {
        "success": out.get("success", 0),
        "sig_p05": 0,
        "sig_dpsi_02": 0,
        "sig_dpsi_01": out.get("sig_dpsi_01", 0),
        "sig_dpsi_005": 0,
        "unchanged_dpsi_lt_005": out.get("unchanged_dpsi_lt_005", 0),
        "not_informative": out.get("not_informative", 0),
    }


def collect_leafcutter_units(leafcutter_dir: Path) -> List[Tuple[str, Path, Path]]:
    """Collect processable units as (cell_type, sig_file, effect_file).

    Standard layout has one subfolder per cell type. Some folders (e.g. fibroblast)
    provide files directly at leafcutter root; these are treated as one pseudo cell
    type named "Fibroblast".
    """
    units: List[Tuple[str, Path, Path]] = []

    for d in sorted(
        p for p in leafcutter_dir.iterdir()
        if p.is_dir() and p.name != "genes_figs"
    ):
        sig_file = d / "leafcutter_ds_cluster_significance.txt"
        effect_file = d / "leafcutter_ds_effect_sizes.txt"
        if sig_file.exists() and effect_file.exists():
            units.append((d.name, sig_file, effect_file))

    root_sig = leafcutter_dir / "leafcutter_ds_cluster_significance.txt"
    root_eff = leafcutter_dir / "leafcutter_ds_effect_sizes.txt"
    if root_sig.exists() and root_eff.exists():
        units.append(("Fibroblast", root_sig, root_eff))

    return units


def write_sum_table_and_stacked_bar(
    leafcutter_dir: Path,
    dataset_a: str,
    dataset_b: str,
    extra_cells: Optional[Dict[str, Dict[str, int]]] = None,
) -> Dict[str, Dict[str, int]]:
    """Write sum_table_HN6.txt and return per-cell summary counts for plotting."""
    units = collect_leafcutter_units(leafcutter_dir)

    cell_types: List[str] = []
    summary_by_cell: Dict[str, Dict[str, int]] = {}
    for cell_type, sig_file, effect_file in units:
        cell_types.append(cell_type)
        summary_by_cell[cell_type] = summarize_leafcutter_celltype(sig_file.parent)

    if extra_cells:
        for ct, summary in extra_cells.items():
            if ct not in summary_by_cell:
                cell_types.append(ct)
            summary_by_cell[ct] = summary

    # Prefer dataset-specific order for readability in outputs.
    dataset_pair = {dataset_a, dataset_b}
    if dataset_pair == {"GSE116177", "GSE180020"}:
        preferred_order = ["CD4T", "Cd8T", "BCell", "NK", "Mono"]
    elif dataset_pair == {"GSE115736", "GSE116177"}:
        preferred_order = ["CD4T", "CD8T", "NveB", "NK", "Mono", "Neut", "Fibroblast"]
    else:
        preferred_order = ["CD4T", "CD8T", "NveB", "NK", "Mono", "Neut"]

    display_cell_types = [ct for ct in preferred_order if ct in summary_by_cell] + [
        ct for ct in cell_types if ct not in preferred_order
    ]
    display_labels = [_display_cell_name(ct, dataset_a, dataset_b) for ct in display_cell_types]

    # Write table with requested summary rows.
    sum_table = leafcutter_dir / "sum_table_HN6.txt"
    with sum_table.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["", "All"] + display_labels)

        rows = [
            ("Leafcutter success clusters", "success"),
            ("Leafcutter sig clusters (p≤0.05 and deltapsi≥0.1)", "sig_dpsi_01"),
            ("Unchanged (deltapsi<0.05 or p>0.05)", "unchanged_dpsi_lt_005"),
            ("Not informative", "not_informative"),
        ]

        for label, key in rows:
            row = [label, ""]
            for ct in display_cell_types:
                row.append(f"{float(summary_by_cell[ct][key]):.1f}")
            writer.writerow(row)

    print(f"Wrote {sum_table}")

    ordered_summary = {
        _display_cell_name(ct, dataset_a, dataset_b): summary_by_cell[ct]
        for ct in display_cell_types
    }
    return ordered_summary


def write_combined_stacked_bar(
    summaries_by_pair: Dict[Tuple[str, str], Dict[str, Dict[str, int]]],
) -> None:
    """Write one 3-panel stacked bar figure for HH, HM, and MM comparisons.

    Bars are percentages where each cell type's success clusters are 100%.
    """

    if plt is None:
        print(
            "Skipping stacked bar plot: matplotlib is not installed in the active environment. "
            "Install it to enable plot output."
        )
        return

    panel_order = [
        ("GSE115736", "GSE60424"),
        ("GSE115736", "GSE116177"),
        ("GSE116177", "GSE180020"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    colors = {
        "not_informative": "#deebf7",
        "unchanged": "#9ecae1",
        "sig_dpsi_01": "#74c476",
    }

    for panel_idx, (dataset_a, dataset_b) in enumerate(panel_order):
        ax = axes[panel_idx]
        pair_key = _normalize_pair(dataset_a, dataset_b)
        summary_by_cell = summaries_by_pair.get(pair_key, {})

        display_cell_types = list(summary_by_cell.keys())
        success_vals = [summary_by_cell[ct]["success"] for ct in display_cell_types]
        sig_dpsi_01_vals = [summary_by_cell[ct]["sig_dpsi_01"] for ct in display_cell_types]
        unchanged_vals = [summary_by_cell[ct]["unchanged_dpsi_lt_005"] for ct in display_cell_types]
        not_informative_vals = [summary_by_cell[ct]["not_informative"] for ct in display_cell_types]

        not_informative_pct = [
            ((v / s) * 100.0) if s > 0 else 0.0
            for v, s in zip(not_informative_vals, success_vals)
        ]
        unchanged_pct = [
            ((v / s) * 100.0) if s > 0 else 0.0
            for v, s in zip(unchanged_vals, success_vals)
        ]
        sig_dpsi_01_pct = [
            ((v / s) * 100.0) if s > 0 else 0.0
            for v, s in zip(sig_dpsi_01_vals, success_vals)
        ]
        bottom2 = [a + b for a, b in zip(not_informative_pct, unchanged_pct)]

        x = list(range(len(display_cell_types)))
        ax.bar(x, not_informative_pct, label="Not informative", color=colors["not_informative"])
        ax.bar(
            x,
            unchanged_pct,
            bottom=not_informative_pct,
            label="Unchanged (deltapsi<0.05 or p>0.05)",
            color=colors["unchanged"],
        )
        ax.bar(
            x,
            sig_dpsi_01_pct,
            bottom=bottom2,
            label="Leafcutter sig clusters (p≤0.05 and deltapsi≥0.1)",
            color=colors["sig_dpsi_01"],
        )

        ax.set_xticks(x)
        ax.set_xticklabels(display_cell_types, rotation=25, ha="right")
        ax.set_ylim(0, 110)
        ax.set_title(_comparison_title(dataset_a, dataset_b), pad=10)

        for i, total in enumerate(success_vals):
            ax.text(i, 102, str(int(total)), ha="center", va="bottom", fontsize=8)

    axes[0].set_ylabel("Percent of success clusters (%)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=max(1, len(labels)),
        fontsize=9,
        frameon=False,
    )
    fig.suptitle("Leafcutter cluster outcomes by comparison (success = 100%)", y=0.98)

    plt.tight_layout(rect=(0, 0.08, 1, 0.93))
    plot_png = SHARED_JUNCTIONS_DIR / "combined_stacked_success_sig_deltapsi01_HH_HM_MM.png"
    plot_svg = SHARED_JUNCTIONS_DIR / "combined_stacked_success_sig_deltapsi01_HH_HM_MM.svg"
    fig.savefig(plot_png, dpi=200)
    fig.savefig(plot_svg)
    plt.close(fig)

    print(f"Wrote {plot_png}")
    print(f"Wrote {plot_svg}")


def write_combined_stacked_bar_possible_as(
    summaries_by_pair: Dict[Tuple[str, str], Dict[str, Dict[str, int]]],
    possible_as_clusters: int = POSSIBLE_AS_CLUSTERS,
) -> None:
    """Write one 3-panel stacked bar figure where 100% equals possible AS clusters."""

    if plt is None:
        print(
            "Skipping stacked bar plot: matplotlib is not installed in the active environment. "
            "Install it to enable plot output."
        )
        return

    if possible_as_clusters <= 0:
        print("Skipping possible-AS stacked bar plot: denominator must be positive.")
        return

    panel_order = [
        ("GSE115736", "GSE60424"),
        ("GSE115736", "GSE116177"),
        ("GSE116177", "GSE180020"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    colors = {
        "not_success": "#f0f0f0",
        "not_informative": "#deebf7",
        "unchanged": "#9ecae1",
        "sig_dpsi_01": "#74c476",
    }

    for panel_idx, (dataset_a, dataset_b) in enumerate(panel_order):
        ax = axes[panel_idx]
        pair_key = _normalize_pair(dataset_a, dataset_b)
        summary_by_cell = summaries_by_pair.get(pair_key, {})

        display_cell_types = list(summary_by_cell.keys())
        success_vals = [summary_by_cell[ct]["success"] for ct in display_cell_types]
        not_success_vals = [max(0, possible_as_clusters - s) for s in success_vals]
        sig_dpsi_01_vals = [summary_by_cell[ct]["sig_dpsi_01"] for ct in display_cell_types]
        unchanged_vals = [summary_by_cell[ct]["unchanged_dpsi_lt_005"] for ct in display_cell_types]
        not_informative_vals = [summary_by_cell[ct]["not_informative"] for ct in display_cell_types]

        denom = float(possible_as_clusters)
        not_success_pct = [(v / denom) * 100.0 for v in not_success_vals]
        not_informative_pct = [(v / denom) * 100.0 for v in not_informative_vals]
        unchanged_pct = [(v / denom) * 100.0 for v in unchanged_vals]
        sig_dpsi_01_pct = [(v / denom) * 100.0 for v in sig_dpsi_01_vals]

        bottom_not_informative = not_success_pct
        bottom_unchanged = [a + b for a, b in zip(not_success_pct, not_informative_pct)]
        bottom_sig = [a + b for a, b in zip(bottom_unchanged, unchanged_pct)]

        x = list(range(len(display_cell_types)))
        ax.bar(x, not_success_pct, label="Not success", color=colors["not_success"])
        ax.bar(
            x,
            not_informative_pct,
            bottom=bottom_not_informative,
            label="Not informative",
            color=colors["not_informative"],
        )
        ax.bar(
            x,
            unchanged_pct,
            bottom=bottom_unchanged,
            label="Unchanged (deltapsi<0.05 or p>0.05)",
            color=colors["unchanged"],
        )
        ax.bar(
            x,
            sig_dpsi_01_pct,
            bottom=bottom_sig,
            label="Leafcutter sig clusters (p≤0.05 and deltapsi≥0.1)",
            color=colors["sig_dpsi_01"],
        )

        ax.set_xticks(x)
        ax.set_xticklabels(display_cell_types, rotation=25, ha="right")
        ax.set_title(_comparison_title(dataset_a, dataset_b), pad=10)

        for i, total in enumerate(success_vals):
            ax.text(i, 102, str(int(total)), ha="center", va="bottom", fontsize=8)

    for ax in axes:
        ax.set_ylim(0, 110)

    axes[0].set_ylabel(f"Percent of possible AS clusters (n={possible_as_clusters})")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=max(1, len(labels)),
        fontsize=9,
        frameon=False,
    )
    fig.suptitle(
        f"Leafcutter cluster outcomes by comparison (possible AS clusters = {possible_as_clusters})",
        y=0.98,
    )

    plt.tight_layout(rect=(0, 0.08, 1, 0.93))
    plot_png = SHARED_JUNCTIONS_DIR / "combined_stacked_possibleAS1823_sig_deltapsi01_HH_HM_MM.png"
    plot_svg = SHARED_JUNCTIONS_DIR / "combined_stacked_possibleAS1823_sig_deltapsi01_HH_HM_MM.svg"
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

    for ct, sig_file, effect_file in collect_leafcutter_units(leafcutter_dir):
        cell_types.append(ct)
        cell_sig[ct] = load_cluster_significance(sig_file)
        cell_eff[ct] = load_effect_sizes_by_junction(effect_file)

    value_file = _value_file_for_pair(dataset_a, dataset_b)
    all_junctions = load_junctions_with_as_averages(value_file, dataset_a, dataset_b, cell_types)

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
    fibro_dir = SHARED_JUNCTIONS_DIR / "leafcutter_EMTAB5919H_EMTAB5919M"
    fibro_summary = None
    if fibro_dir.exists():
        fibro_summary = summarize_leafcutter_celltype(fibro_dir)
        if fibro_summary.get("success", 0) == 0:
            fibro_summary = load_summary_from_sum_table(fibro_dir / "sum_table_HN6.txt") or fibro_summary

    leafcutter_dirs = sorted(d for d in SHARED_JUNCTIONS_DIR.iterdir()
                             if d.is_dir() and d.name.startswith("leafcutter_"))

    summaries_by_pair: Dict[Tuple[str, str], Dict[str, Dict[str, int]]] = {}

    for leaf_dir in leafcutter_dirs:
        match = re.match(r"leafcutter_([^_]+)_([^_]+)$", leaf_dir.name)
        if not match:
            continue
        dataset_a, dataset_b = match.group(1), match.group(2)

        print(f"Processing {leaf_dir.name} ({dataset_a} vs {dataset_b})...")
        merge_leafcutter_results(leaf_dir, dataset_a, dataset_b, metadata)

        extra_cells = None
        if _normalize_pair(dataset_a, dataset_b) == _normalize_pair("GSE115736", "GSE116177") and fibro_summary:
            extra_cells = {"Fibroblast": fibro_summary}

        summary = write_sum_table_and_stacked_bar(leaf_dir, dataset_a, dataset_b, extra_cells=extra_cells)
        summaries_by_pair[_normalize_pair(dataset_a, dataset_b)] = summary

    write_combined_stacked_bar(summaries_by_pair)
    write_combined_stacked_bar_possible_as(summaries_by_pair)


if __name__ == "__main__":
    main()
