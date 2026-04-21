#!/usr/bin/env python3
"""Create bar plots for each Success cluster in each cell type and comparison.

Input sources:
- AS value table with junction rows: AS_clusters_value_HN6.txt
- Leafcutter comparison folders: leafcutter_GSE*_* / <cell_type> /
  - groups_file.txt
  - leafcutter_ds_cluster_significance.txt

For every (comparison, cell_type, success_cluster):
- Build stacked PSI bar plot (top panel) and stacked counts bar plot (bottom panel)
- Save SVG and PNG under: <comparison>/genes_figs/<cell_type>/
"""

import argparse
import os
import re
import shutil
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_BASE_DIR = Path("/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026")
#DEFAULT_VALUE_FILE = DEFAULT_BASE_DIR / "AS_clusters_value_HN6.txt"
DEFAULT_VALUE_FILE = DEFAULT_BASE_DIR / "AS_clusters_value_fibroblast_HN6.txt"

def create_color_map(intron_legend: List[str], legend_color_map: Dict[str, tuple]) -> Dict[str, tuple]:
    """Create/update a stable color map for introns across plots."""
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


def extract_dataset_id(sample_name: str) -> str:
    upper_name = sample_name.upper()
    if upper_name.startswith("HS"):
        return "HS"
    if upper_name.startswith("MM"):
        return "MM"

    match = re.search(r"(GSE\d+)$", sample_name)
    return match.group(1) if match else ""


def build_positions_by_dataset(
    sample_names: List[str],
    within_step: float = 0.25,
    between_gap: float = 0.55,
) -> np.ndarray:
    """Return x positions with tight spacing within dataset and extra gap between datasets."""
    if not sample_names:
        return np.array([])

    positions = np.zeros(len(sample_names), dtype=float)
    datasets = [extract_dataset_id(s) for s in sample_names]

    for i in range(1, len(sample_names)):
        positions[i] = positions[i - 1] + within_step
        if datasets[i] != datasets[i - 1]:
            positions[i] += between_gap

    return positions


def build_xtick_labels(sample_names: List[str]) -> List[str]:
    """Show one label per contiguous sample group block (e.g. HS ... MM)."""
    out = [" " for _ in sample_names]
    prev_group = None

    for i, sample_name in enumerate(sample_names):
        group = extract_dataset_id(sample_name)
        if group != prev_group:
            out[i] = group if group else " "
            prev_group = group

    return out


def sanitize_filename(text: str) -> str:
    text = text.replace("$cl$", "_")
    text = re.sub(r"[^A-Za-z0-9._:=+-]+", "_", text)
    return re.sub(r"_+", "_", text)


def plot_cluster_bars(
    psi_cluster: pd.DataFrame,
    value_cluster: pd.DataFrame,
    fig_title: str,
    output_dir: Path,
    legend_color_map: Dict[str, tuple],
) -> Dict[str, tuple]:
    """Plot stacked PSI + counts bar chart for one cluster."""
    output_svg = output_dir / f"{sanitize_filename(fig_title)}.svg"
    output_png = output_dir / f"{sanitize_filename(fig_title)}.png"

    psi_cluster = psi_cluster.sort_index()
    value_cluster = value_cluster.sort_index()

    intron_legend = psi_cluster.index.tolist()
    legend_color_map = create_color_map(intron_legend, legend_color_map)

    original_samples = psi_cluster.columns.tolist()

    # Keep count and PSI column order identical and use shared x positions.
    value_cluster = value_cluster[original_samples]

    labels = build_xtick_labels(original_samples)
    r = build_positions_by_dataset(original_samples)
    n = len(original_samples)
    bar_width = 0.2

    introns_psi = psi_cluster.values
    introns_counts = value_cluster.values

    bottom_psi = np.zeros(n)
    bottom_counts = np.zeros(n)

    plt.close("all")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 4), gridspec_kw={"height_ratios": [1, 1]})

    for i in range(len(introns_psi)):
        color = legend_color_map[intron_legend[i]]
        ax1.bar(r, introns_psi[i], width=bar_width, bottom=bottom_psi, edgecolor="none", color=color)
        bottom_psi = np.add(bottom_psi, introns_psi[i])

        ax2.bar(r, introns_counts[i], width=bar_width, bottom=bottom_counts, edgecolor="none", color=color)
        bottom_counts = np.add(bottom_counts, introns_counts[i])

    ax1.set_xticks(r)
    ax1.set_xticklabels(labels, fontsize=0)
    ax1.tick_params(axis="x", which="both", length=0)
    ax1.set_ylim(0, 1)
    ax1.set_yticks([0, 1])
    ax1.set_ylabel("PSI")

    ax2.set_xticks(r)
    ax2.set_xticklabels(labels, fontsize=8, rotation=45, ha="right")
    ax2.tick_params(axis="x", which="both", length=0)
    ax2.set_ylabel("Counts")

    fig.legend(intron_legend, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1)
    fig.suptitle(fig_title, fontsize=10)
    plt.subplots_adjust(hspace=0.35)

    fig.savefig(output_svg, bbox_inches="tight")
    fig.savefig(output_png, bbox_inches="tight")
    plt.close("all")

    return legend_color_map


def load_as_value_table(value_file: Path) -> pd.DataFrame:
    """Load AS value table and append cluster ID column from row index."""
    df = pd.read_csv(value_file, sep="\t", index_col=0)

    idx = df.index.astype(str).tolist()
    clusters = []
    for junction_id in idx:
        if ":" not in junction_id:
            clusters.append("")
            continue
        junc_part, clu_suffix = junction_id.rsplit(":", 1)
        chrom = junc_part.split(":", 1)[0]
        clusters.append(f"{chrom}:{clu_suffix}")

    df["cluster"] = clusters
    return df


def iter_cell_units(comp_dir: Path) -> List[tuple[str, Path, Path, Path]]:
    """Return processing units as (cell_name, output_dir, groups_file, sig_file)."""
    units: List[tuple[str, Path, Path, Path]] = []

    # Standard layout: one folder per cell type.
    for cell_dir in sorted([d for d in comp_dir.iterdir() if d.is_dir() and d.name != "genes_figs"]):
        groups_file = cell_dir / "groups_file.txt"
        sig_file = cell_dir / "leafcutter_ds_cluster_significance.txt"
        if groups_file.exists() and sig_file.exists():
            units.append((cell_dir.name, cell_dir, groups_file, sig_file))

    # EMTAB layout: files are directly in comparison root (no cell subfolders).
    root_groups = comp_dir / "groups_file.txt"
    root_sig = comp_dir / "leafcutter_ds_cluster_significance.txt"
    if root_groups.exists() and root_sig.exists():
        units.append(("Fibroblast", comp_dir / "Fibroblast", root_groups, root_sig))

    return units


def run(base_dir: Path, value_file: Path, clean: bool, max_clusters: int) -> None:
    counts_df = load_as_value_table(value_file)

    comparison_dirs = sorted([p for p in base_dir.iterdir() if p.is_dir() and p.name.startswith("leafcutter_")])
    if not comparison_dirs:
        raise FileNotFoundError(f"No leafcutter_* folders found in {base_dir}")

    for comp_dir in comparison_dirs:
        print(f"Processing comparison: {comp_dir.name}")
        genes_figs_root = comp_dir / "genes_figs"
        if clean and genes_figs_root.exists():
            shutil.rmtree(genes_figs_root)
        genes_figs_root.mkdir(parents=True, exist_ok=True)

        for cell_name, out_dir, groups_file, sig_file in iter_cell_units(comp_dir):
            print(f"  Cell type: {cell_name}")

            # groups_file can be tab or space separated; use regex whitespace.
            group_df = pd.read_csv(groups_file, sep=r"\s+", header=None, names=["sample", "condition"], engine="python")
            sample_cols = [s for s in group_df["sample"].astype(str).tolist() if s in counts_df.columns]
            if not sample_cols:
                print(f"    Skipping: no sample columns matched in AS value table")
                continue

            value_table = counts_df[sample_cols + ["cluster"]].copy()

            # PSI per sample = junction count / total cluster count per sample
            percent_table = value_table.copy()
            totals = value_table.groupby("cluster")[sample_cols].transform("sum")
            percent_table[sample_cols] = value_table[sample_cols].div(totals).fillna(0)

            sig_df = pd.read_csv(sig_file, sep="\t")
            success_clusters = sig_df.loc[sig_df["status"] == "Success", "cluster"].astype(str).tolist()
            if max_clusters > 0:
                success_clusters = success_clusters[:max_clusters]

            out_dir = genes_figs_root / out_dir.name
            out_dir.mkdir(parents=True, exist_ok=True)

            legend_color_map: Dict[str, tuple] = {}
            plotted = 0
            for cluster in success_clusters:
                psi_cluster = percent_table.loc[percent_table["cluster"] == cluster].copy()
                val_cluster = value_table.loc[value_table["cluster"] == cluster].copy()

                if psi_cluster.empty or val_cluster.empty:
                    continue

                psi_cluster = psi_cluster.drop(columns=["cluster"])
                val_cluster = val_cluster.drop(columns=["cluster"])

                sig_row = sig_df.loc[sig_df["cluster"] == cluster]
                gene_name = ""
                p_adjust = "NA"
                if not sig_row.empty:
                    row = sig_row.iloc[0]
                    gene_name = str(row.get("genes", ""))
                    p_raw = row.get("p.adjust", np.nan)
                    p_adjust = f"{float(p_raw):.3g}" if pd.notna(p_raw) else "NA"

                title = f"{cell_name}_{gene_name}_{cluster}_p={p_adjust}" if gene_name else f"{cell_name}_{cluster}_p={p_adjust}"
                legend_color_map = plot_cluster_bars(psi_cluster, val_cluster, title, out_dir, legend_color_map)
                plotted += 1

            print(f"    Plotted {plotted} success clusters")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot stacked PSI and counts bars for Success clusters.")
    parser.add_argument("--base-dir", type=Path, default=DEFAULT_BASE_DIR, help="Directory containing leafcutter_* folders")
    parser.add_argument("--value-file", type=Path, default=DEFAULT_VALUE_FILE, help="AS cluster value file")
    parser.add_argument("--clean", action="store_true", help="Remove existing genes_figs output before plotting")
    parser.add_argument("--max-clusters", type=int, default=0, help="Limit plots per cell type (0 = all)")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    run(args.base_dir, args.value_file, args.clean, args.max_clusters)


if __name__ == "__main__":
    main()
