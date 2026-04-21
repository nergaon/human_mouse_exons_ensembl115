#!/usr/bin/env python3
"""Plot stacked sample-level PSI/count bars for 4 datasets.

Style matches plot_bar_success_clusters_HN6.py:
- x-axis uses samples
- upper panel: stacked PSI
- lower panel: stacked counts
- samples from the same dataset are grouped together

Outputs:
- Per-cell plots under: <base>/4dbPlots/<cell>/
    Filename: <gene>_<cluster>_<cell>.png
- Combined all-cell plots under: <base>/4dbPlots/all_cell_types/
    Filename: <gene>_<cluster>_all_cells.png
"""

import argparse
import re
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import transforms
from matplotlib.patches import Patch


DEFAULT_BASE_DIR = Path("/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026")
DEFAULT_OUTPUT_DIR = DEFAULT_BASE_DIR / "4dbPlots"

DATASET_ORDER = ["GSE180020", "GSE116177", "GSE115736", "GSE60424"]
DATASET_LABELS = {
    "GSE180020": "M2",
    "GSE116177": "M1",
    "GSE115736": "H1",
    "GSE60424": "H2",
}
COMPARISON_LABELS = {
    "leafcutter_GSE115736_GSE116177": "H1_M1",
    "leafcutter_GSE115736_GSE60424": "H1_H2",
    "leafcutter_GSE116177_GSE180020": "M1_M2",
}
COMPARISON_DATASETS = [
    ("H1_M1", ("GSE115736", "GSE116177")),
    ("H1_H2", ("GSE115736", "GSE60424")),
    ("M1_M2", ("GSE116177", "GSE180020")),
]
CELL_ORDER = ["Mono", "NveB", "NK", "CD8", "CD4"]

VALUE_FILE = DEFAULT_BASE_DIR / "AS_clusters_value_HN6.txt"
GENE_TABLE_FILES = [
    DEFAULT_BASE_DIR / "leafcutter_GSE115736_GSE116177" / "clusters_sum_table_HN6.txt",
    DEFAULT_BASE_DIR / "leafcutter_GSE115736_GSE60424" / "clusters_sum_table_HN6.txt",
    DEFAULT_BASE_DIR / "leafcutter_GSE116177_GSE180020" / "clusters_sum_table_HN6.txt",
]


def sanitize_filename(text: str) -> str:
    text = str(text).strip().replace(" ", "_")
    text = text.replace(":", "_")
    text = re.sub(r"[^A-Za-z0-9._+-]+", "_", text)
    return re.sub(r"_+", "_", text).strip("_") or "NA"


def first_gene_name(genes_value: object) -> str:
    genes = str(genes_value).strip()
    if not genes or genes.lower() == "nan":
        return "NA"
    # If multiple genes are comma-separated, keep only the first.
    return genes.split(",", 1)[0].strip() or "NA"


def parse_dataset(sample_name: str) -> str:
    m = re.search(r"_(GSE\d+)$", sample_name)
    return m.group(1) if m else ""


def normalize_cell_name(cell_name: str) -> str:
    value = str(cell_name)
    if value in {"BCell", "Bcells", "NveB"}:
        return "NveB"
    if value in {"CD4T", "CD4"}:
        return "CD4"
    if value in {"Cd8T", "CD8T", "CD8"}:
        return "CD8"
    if value in {"Mono", "NK", "Neut"}:
        return value
    return value


def natural_sort_key(text: str) -> List[object]:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", text)]


def infer_celltype_candidates(sample_name: str) -> List[str]:
    base = re.sub(r"_GSE\d+$", "", sample_name)

    if base.startswith("NveB"):
        return ["NveB"]
    if base.startswith(("BCell", "Bcells", "B_fo")):
        return ["NveB"]

    if base.startswith(("CD4T", "EffCD4T", "MemCD4T", "NveCD4T", "CD4_", "T_4_")):
        return ["CD4"]

    if base.startswith(("CD8T", "MemCD8T", "NveCd8T", "Cd8T", "CD8_", "T_8_")):
        return ["CD8"]

    if base.startswith(("Mono", "InfMono", "Monocytes", "Mo_6C+")):
        return ["Mono"]

    if base.startswith("NK"):
        return ["NK"]

    return []


def sample_short_label(sample_name: str) -> str:
    return re.sub(r"_GSE\d+$", "", sample_name)


def dataset_display_label(dataset_id: str) -> str:
    return DATASET_LABELS.get(dataset_id, dataset_id)


def intron_legend_label(intron_id: str) -> str:
    # Remove trailing cluster suffix from row ids like chr:start:end:clu_123.
    return re.sub(r":clu_[^:]+$", "", str(intron_id))


def extract_cluster_id(row_id: str) -> str:
    if ":" not in row_id:
        return ""
    junc, clu = row_id.rsplit(":", 1)
    chrom = junc.split(":", 1)[0]
    return f"{chrom}:{clu}"


def load_gene_map(gene_table_files: List[Path]) -> Dict[str, str]:
    gene_map: Dict[str, str] = {}
    for path in gene_table_files:
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t", usecols=["cluster", "genes"])
        for _, row in df.dropna(subset=["cluster"]).iterrows():
            cl = str(row["cluster"])
            if cl in gene_map:
                continue
            gene_map[cl] = first_gene_name(row.get("genes", ""))
    return gene_map


def load_p_adjust_map(gene_table_files: List[Path]) -> Dict[str, Dict[str, Dict[str, float]]]:
    p_adjust_map: Dict[str, Dict[str, Dict[str, float]]] = {}

    for path in gene_table_files:
        if not path.exists():
            continue

        comparison_label = COMPARISON_LABELS.get(path.parent.name, path.parent.name)
        df = pd.read_csv(path, sep="\t")
        p_adjust_columns = [col for col in df.columns if col.endswith("_p.adjust")]
        if "cluster" not in df.columns or not p_adjust_columns:
            continue

        for _, row in df.dropna(subset=["cluster"]).iterrows():
            cluster = str(row["cluster"])
            cluster_map = p_adjust_map.setdefault(cluster, {})
            comparison_map = cluster_map.setdefault(comparison_label, {})

            for col in p_adjust_columns:
                cell_name = normalize_cell_name(col[: -len("_p.adjust")])
                p_value = pd.to_numeric(row.get(col), errors="coerce")
                if pd.isna(p_value):
                    continue
                comparison_map[cell_name] = float(p_value)

    return p_adjust_map


def load_abs_deltapsi_map(gene_table_files: List[Path]) -> Dict[str, Dict[str, Dict[str, float]]]:
    abs_deltapsi_map: Dict[str, Dict[str, Dict[str, float]]] = {}

    for path in gene_table_files:
        if not path.exists():
            continue

        comparison_label = COMPARISON_LABELS.get(path.parent.name, path.parent.name)
        df = pd.read_csv(path, sep="\t")
        abs_deltapsi_columns = [col for col in df.columns if col.endswith("_abs_deltapsi")]
        if "cluster" not in df.columns or not abs_deltapsi_columns:
            continue

        for _, row in df.dropna(subset=["cluster"]).iterrows():
            cluster = str(row["cluster"])
            cluster_map = abs_deltapsi_map.setdefault(cluster, {})
            comparison_map = cluster_map.setdefault(comparison_label, {})

            for col in abs_deltapsi_columns:
                cell_name = normalize_cell_name(col[: -len("_abs_deltapsi")])
                delta_value = pd.to_numeric(row.get(col), errors="coerce")
                if pd.isna(delta_value):
                    continue
                existing_value = comparison_map.get(cell_name)
                if existing_value is None or float(delta_value) > existing_value:
                    comparison_map[cell_name] = float(delta_value)

    return abs_deltapsi_map


def format_p_adjust(value: float | None) -> str:
    if value is None or pd.isna(value):
        return "NA"
    return f"{value:.2f}"


def lookup_p_adjust(
    cluster: str,
    cell: str,
    comparison_label: str,
    p_adjust_map: Dict[str, Dict[str, Dict[str, float]]],
) -> float | None:
    return p_adjust_map.get(cluster, {}).get(comparison_label, {}).get(cell)


def lookup_abs_deltapsi(
    cluster: str,
    cell: str,
    comparison_label: str,
    abs_deltapsi_map: Dict[str, Dict[str, Dict[str, float]]],
) -> float | None:
    return abs_deltapsi_map.get(cluster, {}).get(comparison_label, {}).get(cell)


def build_stat_label(
    cluster: str,
    cell: str,
    comparison_label: str,
    p_adjust_map: Dict[str, Dict[str, Dict[str, float]]],
    abs_deltapsi_map: Dict[str, Dict[str, Dict[str, float]]],
) -> str:
    p_adjust = format_p_adjust(lookup_p_adjust(cluster, cell, comparison_label, p_adjust_map))
    abs_deltapsi = format_p_adjust(lookup_abs_deltapsi(cluster, cell, comparison_label, abs_deltapsi_map))
    return f"{comparison_label}: p={p_adjust} | deltapsi={abs_deltapsi}"


def load_counts_table(value_file: Path) -> pd.DataFrame:
    df = pd.read_csv(value_file, sep="\t", index_col=0)
    idx = df.index.astype(str).tolist()
    df["cluster"] = [extract_cluster_id(x) for x in idx]
    return df


def collect_cell_samples(all_columns: List[str]) -> Dict[str, List[str]]:
    samples_by_cell: Dict[str, List[str]] = {c: [] for c in CELL_ORDER}

    for sample in all_columns:
        ds = parse_dataset(sample)
        if ds not in DATASET_ORDER:
            continue
        candidates = infer_celltype_candidates(sample)
        if not candidates:
            continue
        cell = candidates[0]
        if cell in samples_by_cell:
            samples_by_cell[cell].append(sample)

    for cell in CELL_ORDER:
        samples = samples_by_cell[cell]
        samples.sort(key=lambda s: (DATASET_ORDER.index(parse_dataset(s)), natural_sort_key(s)))
        samples_by_cell[cell] = samples

    return samples_by_cell


def create_color_map(intron_legend: List[str], legend_color_map: Dict[str, tuple]) -> Dict[str, tuple]:
    cmap = plt.get_cmap("tab20")
    color_index = 0
    used_colors = [legend_color_map[i] for i in intron_legend if i in legend_color_map]

    for intron in intron_legend:
        if intron in legend_color_map:
            continue
        while cmap(color_index % cmap.N) in used_colors:
            color_index += 1
        legend_color_map[intron] = cmap(color_index % cmap.N)
        used_colors.append(legend_color_map[intron])
        color_index += 1

    return legend_color_map


def build_positions_with_gaps(samples: List[str], cell_break: bool = False) -> np.ndarray:
    if not samples:
        return np.array([])

    within_step = 0.23
    dataset_gap = 0.5
    cell_gap = 0.95

    positions = np.zeros(len(samples), dtype=float)
    prev_ds = parse_dataset(samples[0])
    prev_cell = infer_celltype_candidates(samples[0])[0] if infer_celltype_candidates(samples[0]) else ""

    for i in range(1, len(samples)):
        positions[i] = positions[i - 1] + within_step
        cur_ds = parse_dataset(samples[i])
        cur_cell = infer_celltype_candidates(samples[i])[0] if infer_celltype_candidates(samples[i]) else ""
        if cur_ds != prev_ds:
            positions[i] += dataset_gap
        if cell_break and cur_cell != prev_cell:
            positions[i] += cell_gap
        prev_ds = cur_ds
        prev_cell = cur_cell

    return positions


def build_dataset_start_labels(samples: List[str]) -> List[str]:
    labels = ["" for _ in samples]
    prev_ds = ""
    for i, s in enumerate(samples):
        ds = parse_dataset(s)
        if ds != prev_ds:
            labels[i] = ds
            prev_ds = ds
    return labels


def contiguous_blocks(samples: List[str], positions: np.ndarray, key_fn) -> List[tuple[str, float, float, float]]:
    """Return (block_key, start_x, end_x, center_x) for each contiguous block."""
    blocks: List[tuple[str, float, float, float]] = []
    if not samples:
        return blocks

    start = 0
    prev_key = key_fn(samples[0])
    for i in range(1, len(samples) + 1):
        end_block = i == len(samples) or key_fn(samples[i]) != prev_key
        if end_block:
            start_x = float(positions[start])
            end_x = float(positions[i - 1])
            center_x = float((start_x + end_x) / 2.0)
            blocks.append((prev_key, start_x, end_x, center_x))
            if i < len(samples):
                start = i
                prev_key = key_fn(samples[i])
    return blocks


def build_upper_axis_labels_for_cell(
    samples: List[str],
    positions: np.ndarray,
    cluster: str,
    cell: str,
    p_adjust_map: Dict[str, Dict[str, Dict[str, float]]],
    abs_deltapsi_map: Dict[str, Dict[str, Dict[str, float]]],
) -> tuple[List[float], List[str]]:
    dataset_centers = {dataset: center for dataset, _, _, center in contiguous_blocks(samples, positions, parse_dataset)}
    ticks: List[float] = []
    labels: List[str] = []

    for comparison_label, datasets in COMPARISON_DATASETS:
        centers = [dataset_centers[dataset] for dataset in datasets if dataset in dataset_centers]
        if not centers:
            continue
        ticks.append(float(sum(centers) / len(centers)))
        labels.append(build_stat_label(cluster, cell, comparison_label, p_adjust_map, abs_deltapsi_map))

    return ticks, labels


def build_upper_axis_labels_for_combined(
    samples: List[str],
    positions: np.ndarray,
    cluster: str,
    p_adjust_map: Dict[str, Dict[str, Dict[str, float]]],
    abs_deltapsi_map: Dict[str, Dict[str, Dict[str, float]]],
) -> tuple[List[float], List[str]]:
    cell_blocks = contiguous_blocks(samples, positions, lambda s: infer_celltype_candidates(s)[0] if infer_celltype_candidates(s) else "")
    ticks: List[float] = []
    labels: List[str] = []

    for cell, _, _, center in cell_blocks:
        ticks.append(center)
        labels.append(
            "\n".join(
                [cell] + [
                    build_stat_label(cluster, cell, comparison_label, p_adjust_map, abs_deltapsi_map)
                    for comparison_label, _ in COMPARISON_DATASETS
                ]
            )
        )

    return ticks, labels


def plot_stacked_sample_panels(
    psi_cluster: pd.DataFrame,
    counts_cluster: pd.DataFrame,
    fig_title: str,
    output_path: Path,
    legend_color_map: Dict[str, tuple],
    upper_xticks: List[float],
    upper_xticklabels: List[str],
    combined_mode: bool = False,
) -> Dict[str, tuple]:
    psi_cluster = psi_cluster.sort_index()
    counts_cluster = counts_cluster.loc[psi_cluster.index, psi_cluster.columns]

    intron_legend = psi_cluster.index.astype(str).tolist()
    legend_labels = [intron_legend_label(x) for x in intron_legend]
    legend_color_map = create_color_map(intron_legend, legend_color_map)

    samples = psi_cluster.columns.tolist()
    if not samples:
        return legend_color_map

    positions = build_positions_with_gaps(samples, cell_break=combined_mode)
    n = len(samples)
    bar_width = 0.19

    introns_psi = psi_cluster.values
    introns_counts = counts_cluster.values

    bottom_psi = np.zeros(n)
    bottom_counts = np.zeros(n)

    plt.close("all")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 6), gridspec_kw={"height_ratios": [1, 1]})
    legend_handles = []

    for i, intron in enumerate(intron_legend):
        color = legend_color_map[intron]
        ax1.bar(positions, introns_psi[i], width=bar_width, bottom=bottom_psi, edgecolor="none", color=color)
        bottom_psi = np.add(bottom_psi, introns_psi[i])

        ax2.bar(positions, introns_counts[i], width=bar_width, bottom=bottom_counts, edgecolor="none", color=color)
        bottom_counts = np.add(bottom_counts, introns_counts[i])
        legend_handles.append(Patch(facecolor=color, edgecolor="none", label=legend_labels[i]))

    # Keep bars close to panel borders (minimal left/right padding).
    x_min = float(np.min(positions) - bar_width * 0.55)
    x_max = float(np.max(positions) + bar_width * 0.55)
    ax1.set_xlim(x_min, x_max)
    ax2.set_xlim(x_min, x_max)

    ax1.set_xticks(upper_xticks)
    ax1.set_xticklabels(upper_xticklabels, fontsize=7, rotation=0, ha="center")
    ax1.tick_params(axis="x", length=0)
    ax1.set_ylim(0, 1)
    ax1.set_yticks([0, 1])
    ax1.set_ylabel("PSI")

    # Bottom panel: in combined mode, x-axis is dataset labels with angle.
    if combined_mode:
        cell_blocks = contiguous_blocks(samples, positions, lambda s: infer_celltype_candidates(s)[0] if infer_celltype_candidates(s) else "")
        ds_blocks = contiguous_blocks(samples, positions, parse_dataset)
        ax2.set_xticks([center for _, _, _, center in ds_blocks])
        ax2.set_xticklabels([dataset_display_label(k) for k, _, _, _ in ds_blocks], fontsize=8, rotation=0, ha="center")
        ax2.tick_params(axis="x", length=0)

        # Make cell boundaries clearer.
        for i in range(1, len(cell_blocks)):
            prev_end = cell_blocks[i - 1][2]
            cur_start = cell_blocks[i][1]
            boundary_x = (prev_end + cur_start) / 2.0
            ax1.axvline(boundary_x, color="#4d4d4d", linewidth=1.2, alpha=0.95)
            ax2.axvline(boundary_x, color="#4d4d4d", linewidth=1.2, alpha=0.95)
    else:
        ax2.set_xticks(positions)
        ax2.set_xticklabels(["" for _ in samples])
        ax2.tick_params(axis="x", length=0)
    ax2.set_ylabel("Counts")

    # Per-cell plots: keep dataset labels centered below axis (horizontal).
    if not combined_mode:
        trans = transforms.blended_transform_factory(ax2.transData, ax2.transAxes)
        for ds, _, _, x_center in contiguous_blocks(samples, positions, parse_dataset):
            ax2.text(x_center, -0.12, dataset_display_label(ds), transform=trans, ha="center", va="top", fontsize=8)

    fig.legend(legend_handles, legend_labels, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.03), ncol=4)
    fig.suptitle(fig_title, fontsize=10)
    plt.subplots_adjust(hspace=0.35)

    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close("all")
    return legend_color_map


def build_cluster_tables(counts_df: pd.DataFrame, cluster: str, sample_cols: List[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not sample_cols:
        return pd.DataFrame(), pd.DataFrame()

    sub = counts_df.loc[counts_df["cluster"] == cluster, sample_cols].copy()
    if sub.empty:
        return pd.DataFrame(), pd.DataFrame()

    # PSI per sample = junction count / total cluster count per sample.
    totals = sub.sum(axis=0)
    psi = sub.div(totals.where(totals != 0, np.nan), axis=1).fillna(0)

    return psi, sub


def run(base_dir: Path, output_dir: Path) -> None:
    value_file = base_dir / "AS_clusters_value_HN6.txt"
    if not value_file.exists():
        raise FileNotFoundError(f"Missing input file: {value_file}")

    counts_df = load_counts_table(value_file)
    gene_map = load_gene_map(GENE_TABLE_FILES)
    p_adjust_map = load_p_adjust_map(GENE_TABLE_FILES)
    abs_deltapsi_map = load_abs_deltapsi_map(GENE_TABLE_FILES)

    all_samples = [c for c in counts_df.columns if c != "cluster"]
    samples_by_cell = collect_cell_samples(all_samples)

    output_dir.mkdir(parents=True, exist_ok=True)
    for cell in CELL_ORDER:
        (output_dir / cell).mkdir(parents=True, exist_ok=True)
    all_cells_dir = output_dir / "all_cell_types"
    all_cells_dir.mkdir(parents=True, exist_ok=True)

    clusters = sorted([c for c in counts_df["cluster"].dropna().astype(str).unique().tolist() if c])

    n_cell_plots = 0
    n_all_cells_plots = 0

    for cluster in clusters:
        gene = gene_map.get(cluster, "NA")
        safe_gene = sanitize_filename(gene)
        safe_cluster = sanitize_filename(cluster)
        cluster_color_map: Dict[str, tuple] = {}

        cluster_has_any = False
        for cell in CELL_ORDER:
            sample_cols = samples_by_cell.get(cell, [])
            if not sample_cols:
                continue
            psi_cluster, counts_cluster = build_cluster_tables(counts_df, cluster, sample_cols)
            if psi_cluster.empty or counts_cluster.empty:
                continue
            if float(counts_cluster.to_numpy().sum()) <= 0:
                continue

            positions = build_positions_with_gaps(sample_cols)
            upper_xticks, upper_xticklabels = build_upper_axis_labels_for_cell(
                sample_cols,
                positions,
                cluster,
                cell,
                p_adjust_map,
                abs_deltapsi_map,
            )
            title = f"{gene}_{cluster}_{cell}"
            out_name = f"{safe_gene}_{safe_cluster}_{cell}.png"
            out_path = output_dir / cell / out_name

            cluster_color_map = plot_stacked_sample_panels(
                psi_cluster=psi_cluster,
                counts_cluster=counts_cluster,
                fig_title=title,
                output_path=out_path,
                legend_color_map=cluster_color_map,
                upper_xticks=upper_xticks,
                upper_xticklabels=upper_xticklabels,
                combined_mode=False,
            )
            n_cell_plots += 1
            cluster_has_any = True

        if not cluster_has_any:
            continue

        combined_samples: List[str] = []
        for cell in CELL_ORDER:
            combined_samples.extend(samples_by_cell.get(cell, []))

        psi_all, counts_all = build_cluster_tables(counts_df, cluster, combined_samples)
        if not psi_all.empty and not counts_all.empty and float(counts_all.to_numpy().sum()) > 0:
            combined_positions = build_positions_with_gaps(combined_samples, cell_break=True)
            upper_xticks_all, upper_xticklabels_all = build_upper_axis_labels_for_combined(
                combined_samples,
                combined_positions,
                cluster,
                p_adjust_map,
                abs_deltapsi_map,
            )
            title_all = f"{gene}_{cluster}_all_cells"
            out_all = all_cells_dir / f"{safe_gene}_{safe_cluster}_all_cells.png"
            cluster_color_map = plot_stacked_sample_panels(
                psi_cluster=psi_all,
                counts_cluster=counts_all,
                fig_title=title_all,
                output_path=out_all,
                legend_color_map=cluster_color_map,
                upper_xticks=upper_xticks_all,
                upper_xticklabels=upper_xticklabels_all,
                combined_mode=True,
            )
            n_all_cells_plots += 1

    print(f"Done. Per-cell plots: {n_cell_plots}")
    print(f"Done. All-cell plots: {n_all_cells_plots}")
    print(f"Output directory: {output_dir}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot sample-level stacked PSI/count bars by cluster and cell type.")
    parser.add_argument("--base-dir", type=Path, default=DEFAULT_BASE_DIR, help="Base sharedJunctions directory")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help="Output folder for plots")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    run(args.base_dir, args.output_dir)


if __name__ == "__main__":
    main()
