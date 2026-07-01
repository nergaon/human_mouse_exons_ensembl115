#!/usr/bin/env python3
"""Build deltapsi heatmaps for immune high-confidence cluster sets.

Outputs 6 heatmaps:
1) CD4T high-confidence clusters
2) CD8T high-confidence clusters
3) NveB high-confidence clusters
4) NK high-confidence clusters
5) Mono high-confidence clusters
6) Mixed high-confidence clusters (>=1 high differential and >=1 high unchanged,
   in different immune cell types)

Data sources:
- unique_sig_clusters_HN6.xlsx, sheet: all_cluster_status
- clusters_sum_table_HN6.xlsx, sheet: deltapsi
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

BASE = Path("/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026")
LEAFCUTTER_DIR = BASE / "leafcutter_GSE115736_GSE116177"
STATUS_XLSX = LEAFCUTTER_DIR / "unique_sig_clusters_HN6.xlsx"
DELTAPSI_XLSX = LEAFCUTTER_DIR / "clusters_sum_table_HN6.xlsx"
OUT_DIR = LEAFCUTTER_DIR / "detapsi_heatmaps"

IMMUNE_CELLS = ["CD4T", "CD8T", "NveB", "NK", "Mono"]

STATUS_HIGH_DIFF = "high confidance differentially spliced"
STATUS_HIGH_UNCH = "high confidance splicing unchanged"


def _load_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    if not STATUS_XLSX.exists():
        raise FileNotFoundError(f"Missing status workbook: {STATUS_XLSX}")
    if not DELTAPSI_XLSX.exists():
        raise FileNotFoundError(f"Missing deltapsi workbook: {DELTAPSI_XLSX}")

    status_df = pd.read_excel(STATUS_XLSX, sheet_name="all_cluster_status")
    deltapsi_df = pd.read_excel(DELTAPSI_XLSX, sheet_name="deltapsi")

    required_status_cols = ["cluster", "genes"] + IMMUNE_CELLS
    missing_status = [c for c in required_status_cols if c not in status_df.columns]
    if missing_status:
        raise ValueError(f"Missing columns in all_cluster_status: {missing_status}")

    required_dpsi_cols = ["cluster", "h_junction", "genes"] + [f"{ct}_abs_deltapsi" for ct in IMMUNE_CELLS]
    missing_dpsi = [c for c in required_dpsi_cols if c not in deltapsi_df.columns]
    if missing_dpsi:
        raise ValueError(f"Missing columns in deltapsi sheet: {missing_dpsi}")

    status_df = status_df.copy()
    status_df["cluster"] = status_df["cluster"].astype(str)

    deltapsi_df = deltapsi_df.copy()
    deltapsi_df["cluster"] = deltapsi_df["cluster"].astype(str)
    for ct in IMMUNE_CELLS:
        col = f"{ct}_abs_deltapsi"
        deltapsi_df[col] = pd.to_numeric(deltapsi_df[col], errors="coerce")

    return status_df, deltapsi_df


def _high_conf_mask(status_col: pd.Series) -> pd.Series:
    return status_col.isin({STATUS_HIGH_DIFF, STATUS_HIGH_UNCH})


def _select_clusters(status_df: pd.DataFrame) -> dict[str, set[str]]:
    selected: dict[str, set[str]] = {}

    for ct in IMMUNE_CELLS:
        mask = _high_conf_mask(status_df[ct])
        selected[ct] = set(status_df.loc[mask, "cluster"].astype(str))

    has_high_diff = status_df[IMMUNE_CELLS].eq(STATUS_HIGH_DIFF)
    has_high_unch = status_df[IMMUNE_CELLS].eq(STATUS_HIGH_UNCH)

    mixed_mask = np.zeros(len(status_df), dtype=bool)
    for i in range(len(status_df)):
        diff_cells = has_high_diff.columns[has_high_diff.iloc[i]].tolist()
        unch_cells = has_high_unch.columns[has_high_unch.iloc[i]].tolist()
        mixed_mask[i] = any(d != u for d in diff_cells for u in unch_cells)

    selected["mixed_high_confidence"] = set(status_df.loc[mixed_mask, "cluster"].astype(str))
    return selected


def _prepare_heatmap_matrix(
    df: pd.DataFrame,
    include_gene_once_per_cluster: bool = False,
) -> tuple[np.ndarray, list[str], list[str], list[float]]:
    col_order = [f"{ct}_abs_deltapsi" for ct in IMMUNE_CELLS]
    mat = df[col_order].to_numpy(dtype=float)

    row_labels: list[str] = []
    cluster_breaks: list[float] = []
    prev_cluster = None
    seen_cluster_gene: set[str] = set()
    for i, (_, row) in enumerate(df.iterrows()):
        cluster = str(row.get("cluster", ""))
        rank_h = str(row.get("rank_h", ""))
        gene = str(row.get("genes", ""))
        if rank_h.lower() == "nan":
            rank_h = ""
        if gene.lower() == "nan":
            gene = ""
        if prev_cluster is not None and cluster != prev_cluster:
            cluster_breaks.append(i - 0.5)

        if include_gene_once_per_cluster:
            gene_part = gene if cluster not in seen_cluster_gene else ""
            seen_cluster_gene.add(cluster)
            row_labels.append(f"{cluster} | {rank_h} | {gene_part}")
        else:
            row_labels.append(f"{cluster} | {rank_h}")
        prev_cluster = cluster

    col_labels = IMMUNE_CELLS
    return mat, row_labels, col_labels, cluster_breaks


def _plot_heatmap(
    mat: np.ndarray,
    row_labels: list[str],
    col_labels: list[str],
    cluster_breaks: list[float],
    title: str,
    out_png: Path,
    row_fontsize: float = 6.0,
) -> None:
    if mat.size == 0:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "No rows to plot", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        fig.suptitle(title)
        fig.tight_layout()
        fig.savefig(out_png, dpi=220, bbox_inches="tight")
        plt.close(fig)
        return

    finite_vals = mat[np.isfinite(mat)]
    vmax = float(np.nanmax(np.abs(finite_vals))) if finite_vals.size else 1.0
    if vmax == 0:
        vmax = 1.0

    fig_h = max(4.0, min(30.0, 2.0 + 0.03 * len(row_labels)))
    fig, ax = plt.subplots(figsize=(6.4, fig_h))

    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", norm=norm, interpolation="nearest")

    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha="right")
    ax.set_ylabel("Cluster | rank_h")

    # Always show row labels; scale font down with row count so labels remain visible.
    if len(row_labels) > 900:
        y_font = min(row_fontsize, 1.6)
    elif len(row_labels) > 600:
        y_font = min(row_fontsize, 1.9)
    elif len(row_labels) > 350:
        y_font = min(row_fontsize, 2.2)
    elif len(row_labels) > 180:
        y_font = min(row_fontsize, 2.8)
    else:
        y_font = row_fontsize

    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=y_font)

    # Draw separators so junctions belonging to the same cluster are visually grouped.
    for y in cluster_breaks:
        ax.axhline(y=y, color="black", linewidth=0.6, alpha=0.55)

    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, shrink=0.86, fraction=0.03, pad=0.02)
    cbar.set_label("deltapsi")

    fig.tight_layout()
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def build_detapsi_heatmaps() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    status_df, deltapsi_df = _load_tables()
    selected = _select_clusters(status_df)

    categories = [
        ("CD4T", "cd4t_high_confidence_deltapsi_heatmap.png"),
        ("CD8T", "cd8t_high_confidence_deltapsi_heatmap.png"),
        ("NveB", "nveb_high_confidence_deltapsi_heatmap.png"),
        ("NK", "nk_high_confidence_deltapsi_heatmap.png"),
        ("Mono", "mono_high_confidence_deltapsi_heatmap.png"),
        ("mixed_high_confidence", "mixed_high_confidence_deltapsi_heatmap.png"),
    ]

    for key, out_name in categories:
        clusters = selected.get(key, set())
        subset = deltapsi_df[deltapsi_df["cluster"].isin(clusters)].copy()
        subset = subset.sort_values(["cluster", "h_junction"], kind="stable")

        include_gene_once = key != "mixed_high_confidence"
        row_fontsize = 4.5 if key == "mixed_high_confidence" else 6.0
        mat, row_labels, col_labels, cluster_breaks = _prepare_heatmap_matrix(
            subset,
            include_gene_once_per_cluster=include_gene_once,
        )
        n_clusters = subset["cluster"].nunique()
        title_key = key if key != "mixed_high_confidence" else "mixed high confidence"
        title = f"{title_key} | clusters={n_clusters}, junctions={len(subset)}"

        out_png = OUT_DIR / out_name
        _plot_heatmap(
            mat,
            row_labels,
            col_labels,
            cluster_breaks,
            title,
            out_png,
            row_fontsize=row_fontsize,
        )
        print(f"Wrote {out_png}")


if __name__ == "__main__":
    build_detapsi_heatmaps()
