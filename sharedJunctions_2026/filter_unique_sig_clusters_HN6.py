#!/usr/bin/env python3
"""
For each cell type (CD4T, CD8T, NveB, NK, Mono), collect rows from
leafcutter_GSE115736_GSE116177/clusters_sum_table_HN6.txt that satisfy:

    Condition 1 (unique-sig):
            sig   in GSE115736_GSE116177  (p.adjust <= 0.05 and abs_deltapsi >= 0.1)
                AND unchanged in GSE115736_GSE60424  (abs_deltapsi < 0.05 or p.adjust > 0.05)
                AND unchanged in GSE116177_GSE180020 (abs_deltapsi < 0.05 or p.adjust > 0.05)

  Condition 2 (unchanged everywhere):
      unchanged in all three comparisons (abs_deltapsi < 0.05 or p.adjust > 0.05)

Cell-type name aliases across files:
    canonical  file-A   file-B   file-C
    CD4T       CD4T     CD4T     CD4T
    CD8T       CD8T     CD8T     Cd8T
    NveB       NveB     NveB     BCell
    NK         NK       NK       NK
    Mono       Mono     Mono     Mono

Output: Excel file with
    - One sheet per cell type containing matching rows with only that cell type's
        columns from file A/B/C (plus identifiers).
    - One "summary" sheet: one row per unique cluster found in ANY cell-type sheet;
        columns = cell types; values = "sig" / "unchanged" based on file A.
"""

from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path("/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026")

FILE_A = BASE / "leafcutter_GSE115736_GSE116177" / "clusters_sum_table_HN6.txt"
FILE_B = BASE / "leafcutter_GSE115736_GSE60424"  / "clusters_sum_table_HN6.txt"
FILE_C = BASE / "leafcutter_GSE116177_GSE180020"  / "clusters_sum_table_HN6.txt"
OUTPUT = BASE / "leafcutter_GSE115736_GSE116177"  / "unique_sig_clusters_HN6.xlsx"
PLOT_STACKED = BASE / "leafcutter_GSE115736_GSE116177" / "unique_sig_clusters_HN6_counts_stacked.png"
PLOT_HEATMAP = BASE / "leafcutter_GSE115736_GSE116177" / "unique_sig_clusters_HN6_summary_heatmap.png"
PLOT_SHARED = BASE / "leafcutter_GSE115736_GSE116177" / "unique_sig_clusters_HN6_shared_sig_counts.png"
PLOT_UPSET_SIG_PNG = BASE / "leafcutter_GSE115736_GSE116177" / "unique_sig_clusters_HN6_upset_sig.png"
PLOT_UPSET_SIG_SVG = BASE / "leafcutter_GSE115736_GSE116177" / "unique_sig_clusters_HN6_upset_sig.svg"
PLOT_UPSET_UNCHANGED_PNG = BASE / "leafcutter_GSE115736_GSE116177" / "unique_sig_clusters_HN6_upset_unchanged.png"
PLOT_UPSET_UNCHANGED_SVG = BASE / "leafcutter_GSE115736_GSE116177" / "unique_sig_clusters_HN6_upset_unchanged.svg"
PLOT_A_ONLY_UPSET_SIG_PNG = BASE / "leafcutter_GSE115736_GSE116177" / "GSE115736_GSE116177_upset_sig.png"
PLOT_A_ONLY_UPSET_SIG_SVG = BASE / "leafcutter_GSE115736_GSE116177" / "GSE115736_GSE116177_upset_sig.svg"
PLOT_A_ONLY_UPSET_UNCHANGED_PNG = BASE / "leafcutter_GSE115736_GSE116177" / "GSE115736_GSE116177_upset_unchanged.png"
PLOT_A_ONLY_UPSET_UNCHANGED_SVG = BASE / "leafcutter_GSE115736_GSE116177" / "GSE115736_GSE116177_upset_unchanged.svg"

THRESHOLD = 0.05
SIG_DELTAPSI_THRESHOLD = 0.1
UNCHANGED_DELTAPSI_THRESHOLD = 0.05

# canonical CT name -> (prefix_in_A, prefix_in_B, prefix_in_C)
CELL_TYPES: dict[str, tuple[str, str, str]] = {
    "CD4T": ("CD4T", "CD4T", "CD4T"),
    "CD8T": ("CD8T", "CD8T", "Cd8T"),
    "NveB": ("NveB", "NveB", "BCell"),
    "NK":   ("NK",   "NK",   "NK"),
    "Mono": ("Mono", "Mono", "Mono"),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def padj_col(prefix: str) -> str:
    return f"{prefix}_p.adjust"

def abs_deltapsi_col(prefix: str) -> str:
    return f"{prefix}_abs_deltapsi"

def cluster_min_padj_lookup(df: pd.DataFrame, prefix: str) -> pd.Series:
    """Return per-cluster minimum p.adjust (index = cluster id, float)."""
    col = padj_col(prefix)
    if col not in df.columns:
        return pd.Series(dtype=float)
    return df.groupby("cluster", dropna=False)[col].min().astype(float)


def cluster_max_abs_deltapsi_lookup(df: pd.DataFrame, prefix: str) -> pd.Series:
    """Return per-cluster maximum absolute deltapsi (index = cluster id, float)."""
    col = abs_deltapsi_col(prefix)
    if col not in df.columns:
        return pd.Series(dtype=float)
    return df[col].astype(float).abs().groupby(df["cluster"], dropna=False).max().astype(float)


def classify_cluster_status(min_padj: float | None, max_abs_dpsi: float | None) -> str:
    """Classify one cluster using one p.adjust and one max abs_deltapsi value."""
    if min_padj is None or max_abs_dpsi is None:
        return ""
    if (min_padj <= THRESHOLD) and (max_abs_dpsi >= SIG_DELTAPSI_THRESHOLD):
        return "sig"
    if (max_abs_dpsi < UNCHANGED_DELTAPSI_THRESHOLD) or (min_padj > THRESHOLD):
        return "unchanged"
    return ""


def cell_type_cols(df: pd.DataFrame, prefix: str) -> list[str]:
    """Return only columns specific to one cell type prefix."""
    cols: list[str] = []
    abs_col = f"{prefix}_abs_deltapsi"
    p_col = f"{prefix}_p.adjust"
    if abs_col in df.columns:
        cols.append(abs_col)
    if p_col in df.columns:
        cols.append(p_col)

    # Keep dataset-specific average columns (e.g., GSE*_avg_CD4T).
    for c in df.columns:
        if c.startswith("GSE") and c.endswith(f"_{prefix}"):
            cols.append(c)
    return cols

def first_gene_label(genes_text: str) -> str:
    """Return the first gene from a delimited genes string."""
    txt = str(genes_text).strip()
    if not txt or txt.lower() == "nan":
        return ""
    return re.split(r"[;,|]", txt, maxsplit=1)[0].strip()

def plot_stacked_counts(counts_by_ct: dict[str, dict[str, int]], output_path: Path) -> None:
    """Plot stacked cluster counts: 5 categories per cell type.

    Categories (stacked bottom to top):
      1. unique-sig         : sig in A, unchanged in B and C  (solid blue)
      2. sig-not-unique     : sig in A, controls NOT all unchanged (hatched blue)
      3. unchanged-all      : unchanged in all three (solid grey)
      4. unchanged-not-all  : unchanged in A, controls NOT all unchanged (hatched grey)
      5. not-informative    : neither sig nor unchanged in A (light beige)
    """
    ct_order = list(CELL_TYPES.keys())
    v_usig  = np.array([counts_by_ct.get(ct, {}).get("unique_sig",        0) for ct in ct_order], dtype=float)
    v_snu   = np.array([counts_by_ct.get(ct, {}).get("sig_not_unique",    0) for ct in ct_order], dtype=float)
    v_uall  = np.array([counts_by_ct.get(ct, {}).get("unchanged_all",     0) for ct in ct_order], dtype=float)
    v_unall = np.array([counts_by_ct.get(ct, {}).get("unchanged_not_all", 0) for ct in ct_order], dtype=float)
    v_noinfo= np.array([counts_by_ct.get(ct, {}).get("not_informative",   0) for ct in ct_order], dtype=float)
    totals  = (v_usig + v_snu + v_uall + v_unall + v_noinfo).astype(int)

    b2 = v_usig
    b3 = b2 + v_snu
    b4 = b3 + v_uall
    b5 = b4 + v_unall

    sig_color  = "#1f77b4"
    unch_color = "#bdbdbd"
    noinfo_color = "#f5f0e8"

    x = np.arange(len(ct_order))
    plt.close("all")
    fig, ax = plt.subplots(figsize=(8.0, 4.8))

    ax.bar(x, v_usig,   color=sig_color,   label="unique-sig (sig in H-M, unchanged in H-H & M-M)")
    ax.bar(x, v_snu,    bottom=b2, color=sig_color,   hatch="///", edgecolor="white",
           label="sig in H-M (controls not all unchanged)")
    ax.bar(x, v_uall,   bottom=b3, color=unch_color,  label="unchanged-all (unchanged in all three)")
    ax.bar(x, v_unall,  bottom=b4, color=unch_color,  hatch="///", edgecolor="white",
           label="unchanged in H-M (controls not all unchanged)")
    ax.bar(x, v_noinfo, bottom=b5, color=noinfo_color, edgecolor="#aaaaaa",
           label="not informative")

    ymax = int(max(totals, default=0))
    pad  = max(5, int(0.015 * ymax))
    for i, total in enumerate(totals):
        ax.text(i, total + pad, str(total), ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(ct_order)
    ax.set_ylabel("Clusters")
    ax.set_title("Cluster counts per cell type (informative clusters in A)")
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close("all")

def plot_summary_heatmap(summary_df: pd.DataFrame, output_path: Path) -> None:
    """Plot cluster x cell-type heatmap: sig=blue, unchanged=light yellow, blank=light grey."""
    ct_order = list(CELL_TYPES.keys())
    if summary_df.empty:
        return

    heat_df = summary_df.copy()
    heat_df["genes"] = heat_df.get("genes", "").fillna("").astype(str)
    for ct in ct_order:
        heat_df[ct] = heat_df[ct].replace("", float("nan")).map({"sig": 1.0, "unchanged": 0.0})
    heat_df["sig_count"] = heat_df[ct_order].sum(axis=1, skipna=True)
    heat_df = heat_df.sort_values(["sig_count", "cluster"], ascending=[False, True]).reset_index(drop=True)

    row_labels = [first_gene_label(gn) for gn in heat_df["genes"]]
    n_rows = len(heat_df)
    n_cols = len(ct_order)

    # map: NaN -> 2 (blank), 0 -> 0 (unchanged), 1 -> 1 (sig)
    mat_raw = heat_df[ct_order].to_numpy(dtype=float)
    mat_coded = np.where(np.isnan(mat_raw), 2.0, mat_raw)

    # custom 3-color map: unchanged=#ffffcc, sig=#2171b5, blank=#d9d9d9
    from matplotlib.colors import ListedColormap, BoundaryNorm
    cmap3 = ListedColormap(["#ffffcc", "#2171b5", "#d9d9d9"])
    norm3 = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap3.N)

    fig_h = max(6.0, min(48.0, 2.0 + 0.035 * n_rows))
    y_font = 6 if n_rows <= 120 else (4 if n_rows <= 300 else 2.5)

    plt.close("all")
    fig, ax = plt.subplots(figsize=(9.4, fig_h))
    ax.pcolormesh(mat_coded, cmap=cmap3, norm=norm3)

    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_xticklabels(ct_order)
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_yticklabels(row_labels, fontsize=y_font)
    ax.set_ylabel(f"Genes (n={n_rows})")
    ax.invert_yaxis()

    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="#2171b5", label="sig"),
        Patch(facecolor="#ffffcc", label="unchanged"),
        Patch(facecolor="#d9d9d9", label="not tested / excluded"),
    ]
    ax.legend(handles=legend_handles, loc="upper right",
              bbox_to_anchor=(1.0, -0.02), frameon=False, fontsize=8)

    ax.set_title("Summary status in A (sig vs unchanged)", pad=30)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close("all")


def plot_shared_sig_counts(summary_df: pd.DataFrame, output_path: Path) -> None:
    """Plot how many clusters are shared by exactly 1..5 significant cell types."""
    if summary_df.empty or "n_sig" not in summary_df.columns:
        return

    x_vals = [1, 2, 3, 4, 5]
    counts = [int((summary_df["n_sig"] == x).sum()) for x in x_vals]

    plt.close("all")
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    bars = ax.bar(x_vals, counts, color="#4C78A8", width=0.72)

    for bar, val in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, val + max(1, int(0.01 * max(counts, default=1))),
                str(val), ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x_vals)
    ax.set_xlabel("Number of cell types with sig")
    ax.set_ylabel("Number of clusters shared")
    ax.set_title("Cluster sharing across cell types")
    fig.tight_layout()
    fig.savefig(output_path, dpi=220)
    plt.close("all")


def _exclusive_intersections_from_summary(
    summary_df: pd.DataFrame,
    status: str,
    ct_order: list[str],
) -> dict[frozenset[str], int]:
    """Build exclusive membership counts for one summary status (sig or unchanged)."""
    intersections: dict[frozenset[str], int] = {}
    if summary_df.empty:
        return intersections

    for _, row in summary_df.iterrows():
        members = frozenset(ct for ct in ct_order if str(row.get(ct, "")) == status)
        if members:
            intersections[members] = intersections.get(members, 0) + 1
    return intersections


def plot_summary_upset(
    summary_df: pd.DataFrame,
    status: str,
    title: str,
    output_png: Path,
    output_svg: Path,
    max_combinations: int = 40,
) -> None:
    """Plot UpSet-style bar + membership matrix using summary-sheet clusters only."""
    ct_order = list(CELL_TYPES.keys())
    intersections = _exclusive_intersections_from_summary(summary_df, status, ct_order)
    if not intersections:
        print(f"  no '{status}' combinations to plot")
        return

    ordered = sorted(
        intersections.items(),
        key=lambda item: (-item[1], -len(item[0]), tuple(sorted(item[0]))),
    )[:max_combinations]

    combinations = [combo for combo, _ in ordered]
    counts = [count for _, count in ordered]

    width = max(10, min(22, len(combinations) * 0.45 + 5.0))
    height = max(5.0, min(12.0, len(ct_order) * 0.6 + 3.5))

    plt.close("all")
    fig = plt.figure(figsize=(width, height))
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 2], hspace=0.04)
    ax_bar = fig.add_subplot(gs[0])
    ax_matrix = fig.add_subplot(gs[1], sharex=ax_bar)

    x = list(range(len(combinations)))
    ax_bar.bar(x, counts, color="#2b8cbe")
    ax_bar.set_ylabel("Clusters")
    ax_bar.set_title(title, pad=10)
    ax_bar.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.5)
    ax_bar.tick_params(axis="x", labelbottom=False)

    pad = max(1, int(0.01 * max(counts, default=1)))
    for idx, value in enumerate(counts):
        ax_bar.text(idx, value + pad, str(value), ha="center", va="bottom", fontsize=7)

    y_map = {cell_type: y for y, cell_type in enumerate(reversed(ct_order))}

    for idx, combo in enumerate(combinations):
        ys = sorted(y_map[cell] for cell in combo)
        for y in y_map.values():
            ax_matrix.plot(idx, y, "o", color="#d9d9d9", markersize=5)
        for y in ys:
            ax_matrix.plot(idx, y, "o", color="#08519c", markersize=6)
        if len(ys) >= 2:
            ax_matrix.plot([idx, idx], [ys[0], ys[-1]], color="#08519c", linewidth=1.5)

    ax_matrix.set_yticks([y_map[ct] for ct in reversed(ct_order)])
    ax_matrix.set_yticklabels(list(reversed(ct_order)))
    ax_matrix.set_ylim(-0.7, len(ct_order) - 0.3)
    ax_matrix.set_xlabel("Shared-cluster combinations (sorted by size)")
    ax_matrix.grid(axis="x", linestyle=":", linewidth=0.4, alpha=0.4)

    fig.tight_layout()
    fig.savefig(output_png, dpi=220)
    fig.savefig(output_svg)
    plt.close("all")


def build_a_only_summary(df_a: pd.DataFrame) -> pd.DataFrame:
    """Build cluster x cell-type status table using only GSE115736_GSE116177 results."""
    cluster_genes = (
        df_a.drop_duplicates(subset="cluster")[["cluster", "genes"]]
        .set_index("cluster")["genes"]
    )

    cluster_status_a: dict[str, dict[str, str]] = {}
    for ct, (pfx_a, _, _) in CELL_TYPES.items():
        min_padj = cluster_min_padj_lookup(df_a, pfx_a)
        max_dpsi = cluster_max_abs_deltapsi_lookup(df_a, pfx_a)
        for cl in set(min_padj.index).union(max_dpsi.index):
            status = classify_cluster_status(min_padj.get(cl), max_dpsi.get(cl))
            if status:
                cluster_status_a.setdefault(str(cl), {})[ct] = status

    rows = []
    for cl in df_a["cluster"].drop_duplicates().astype(str):
        row: dict[str, str] = {"cluster": cl, "genes": str(cluster_genes.get(cl, ""))}
        has_any_status = False
        for ct in CELL_TYPES:
            row[ct] = cluster_status_a.get(cl, {}).get(ct, "")
            has_any_status = has_any_status or bool(row[ct])
        if has_any_status:
            rows.append(row)

    a_only_df = pd.DataFrame(rows, columns=["cluster", "genes"] + list(CELL_TYPES.keys()))
    if not a_only_df.empty:
        ct_cols = list(CELL_TYPES.keys())
        a_only_df["n_sig"] = (a_only_df[ct_cols] == "sig").sum(axis=1)
    return a_only_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print("Loading input files …")
    df_a = pd.read_csv(FILE_A, sep="\t")
    df_b = pd.read_csv(FILE_B, sep="\t")
    df_c = pd.read_csv(FILE_C, sep="\t")

    print(f"  File A: {len(df_a):,} rows")
    print(f"  File B: {len(df_b):,} rows")
    print(f"  File C: {len(df_c):,} rows")

    # Pre-build cluster-level metrics for each comparison.
    min_padj_a: dict[str, pd.Series] = {}
    max_dpsi_a: dict[str, pd.Series] = {}
    min_padj_b: dict[str, pd.Series] = {}
    max_dpsi_b: dict[str, pd.Series] = {}
    min_padj_c: dict[str, pd.Series] = {}
    max_dpsi_c: dict[str, pd.Series] = {}
    for ct, (pfx_a, pfx_b, pfx_c) in CELL_TYPES.items():
        min_padj_a[ct] = cluster_min_padj_lookup(df_a, pfx_a)
        max_dpsi_a[ct] = cluster_max_abs_deltapsi_lookup(df_a, pfx_a)
        min_padj_b[ct] = cluster_min_padj_lookup(df_b, pfx_b)
        max_dpsi_b[ct] = cluster_max_abs_deltapsi_lookup(df_b, pfx_b)
        min_padj_c[ct] = cluster_min_padj_lookup(df_c, pfx_c)
        max_dpsi_c[ct] = cluster_max_abs_deltapsi_lookup(df_c, pfx_c)

    # Ordered list of clusters seen across all cell-type sheets (deduplicated)
    seen_clusters: dict[str, None] = {}  # ordered set via dict
    counts_by_ct: dict[str, dict[str, int]] = {}
    unique_sig_per_ct: dict[str, set[str]] = {ct: set() for ct in CELL_TYPES}
    unchanged_all_per_ct: dict[str, set[str]] = {ct: set() for ct in CELL_TYPES}

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(OUTPUT, engine="openpyxl") as writer:

        # ------------------------------------------------------------------ #
        # Per-cell-type sheets                                                 #
        # ------------------------------------------------------------------ #
        for ct, (pfx_a, pfx_b, pfx_c) in CELL_TYPES.items():
            col_a = padj_col(pfx_a)
            dpsi_col_a = abs_deltapsi_col(pfx_a)
            if col_a not in df_a.columns or dpsi_col_a not in df_a.columns:
                print(
                    f"  [{ct}] columns '{col_a}' and/or '{dpsi_col_a}' not found in file A — skipping"
                )
                continue

            clusters_with_all = (
                set(min_padj_a[ct].index)
                & set(max_dpsi_a[ct].index)
                & set(min_padj_b[ct].index)
                & set(max_dpsi_b[ct].index)
                & set(min_padj_c[ct].index)
                & set(max_dpsi_c[ct].index)
            )

            unique_sig_clusters: set[str] = set()
            sig_not_unique_clusters: set[str] = set()
            unchanged_all_clusters: set[str] = set()
            unchanged_not_all_clusters: set[str] = set()
            not_informative_clusters: set[str] = set()
            for cl in clusters_with_all:
                status_a = classify_cluster_status(min_padj_a[ct].get(cl), max_dpsi_a[ct].get(cl))
                status_b = classify_cluster_status(min_padj_b[ct].get(cl), max_dpsi_b[ct].get(cl))
                status_c = classify_cluster_status(min_padj_c[ct].get(cl), max_dpsi_c[ct].get(cl))

                if (status_a == "sig") and (status_b == "unchanged") and (status_c == "unchanged"):
                    unique_sig_clusters.add(str(cl))
                elif status_a == "sig":
                    sig_not_unique_clusters.add(str(cl))
                elif (status_a == "unchanged") and (status_b == "unchanged") and (status_c == "unchanged"):
                    unchanged_all_clusters.add(str(cl))
                elif status_a == "unchanged":
                    unchanged_not_all_clusters.add(str(cl))
                else:
                    not_informative_clusters.add(str(cl))

            selected_clusters = unique_sig_clusters | unchanged_all_clusters

            # Keep only this cell type columns from each source file.
            keys = [
                c
                for c in ["cluster", "h_junction", "m_junction", "genes", "rank_h", "rank_m"]
                if c in df_a.columns
            ]
            cols_a = keys + cell_type_cols(df_a, pfx_a)
            cols_b = [c for c in ["cluster", "h_junction", "m_junction"] if c in df_b.columns] + cell_type_cols(df_b, pfx_b)
            cols_c = [c for c in ["cluster", "h_junction", "m_junction"] if c in df_c.columns] + cell_type_cols(df_c, pfx_c)

            base = df_a.loc[df_a["cluster"].astype(str).isin(selected_clusters), cols_a].copy()

            merge_keys = [c for c in ["cluster", "h_junction", "m_junction"] if c in base.columns and c in df_b.columns and c in df_c.columns]
            if merge_keys:
                b_small = df_b[cols_b].copy().drop_duplicates(subset=merge_keys)
                c_small = df_c[cols_c].copy().drop_duplicates(subset=merge_keys)
                result = base.merge(b_small, on=merge_keys, how="left", suffixes=("", "_B"))
                result = result.merge(c_small, on=merge_keys, how="left", suffixes=("", "_C"))
            else:
                result = base

            # Add source tag to avoid ambiguous repeated column names.
            rename_map: dict[str, str] = {}
            for c in result.columns:
                if c in keys:
                    continue
                if c in cell_type_cols(df_a, pfx_a):
                    rename_map[c] = f"A_{c}"
                elif c in cell_type_cols(df_b, pfx_b) or c.endswith("_B"):
                    rename_map[c] = f"B_{c.removesuffix('_B')}"
                elif c in cell_type_cols(df_c, pfx_c) or c.endswith("_C"):
                    rename_map[c] = f"C_{c.removesuffix('_C')}"
            result = result.rename(columns=rename_map)

            # Final safety: keep only key identifiers + A/B/C cell-type columns.
            keep_cols = [c for c in keys if c in result.columns] + [
                c for c in result.columns if c.startswith("A_") or c.startswith("B_") or c.startswith("C_")
            ]
            result = result[keep_cols].reset_index(drop=True)
            result.to_excel(writer, sheet_name=ct, index=False)

            n_unique_sig        = len(unique_sig_clusters)
            n_sig_not_unique    = len(sig_not_unique_clusters)
            n_unchanged_all     = len(unchanged_all_clusters)
            n_unchanged_not_all = len(unchanged_not_all_clusters)
            n_not_informative   = len(not_informative_clusters)
            counts_by_ct[ct] = {
                "unique_sig":        n_unique_sig,
                "sig_not_unique":    n_sig_not_unique,
                "unchanged_all":     n_unchanged_all,
                "unchanged_not_all": n_unchanged_not_all,
                "not_informative":   n_not_informative,
            }
            unique_sig_per_ct[ct] = unique_sig_clusters
            unchanged_all_per_ct[ct] = unchanged_all_clusters

            print(f"  [{ct}] {len(result):,} rows written  "
                  f"(unique-sig: {n_unique_sig}, sig-not-unique: {n_sig_not_unique}, "
                  f"unchanged-all: {n_unchanged_all}, unchanged-not-all: {n_unchanged_not_all}, "
                  f"not-informative: {n_not_informative})")

            for cl in sorted(selected_clusters):
                seen_clusters[cl] = None

        # ------------------------------------------------------------------ #
        # Summary sheet                                                        #
        # ------------------------------------------------------------------ #
        print(f"\nBuilding summary sheet for {len(seen_clusters):,} unique clusters …")

        cluster_genes = (
            df_a.drop_duplicates(subset="cluster")[["cluster", "genes"]]
            .set_index("cluster")["genes"]
        )

        rows = []
        for cl in seen_clusters:
            row: dict[str, str] = {"cluster": cl, "genes": str(cluster_genes.get(cl, ""))}
            for ct in CELL_TYPES:
                if cl in unique_sig_per_ct[ct]:
                    row[ct] = "sig"
                elif cl in unchanged_all_per_ct[ct]:
                    row[ct] = "unchanged"
                else:
                    row[ct] = ""
            rows.append(row)

        summary_df = pd.DataFrame(rows, columns=["cluster", "genes"] + list(CELL_TYPES.keys()))
        ct_cols = list(CELL_TYPES.keys())
        summary_df["n_sig"] = (summary_df[ct_cols] == "sig").sum(axis=1)
        summary_df.to_excel(writer, sheet_name="summary", index=False)
        print(f"  summary sheet: {len(summary_df):,} rows")

    plot_stacked_counts(counts_by_ct, PLOT_STACKED)
    plot_summary_heatmap(summary_df, PLOT_HEATMAP)
    plot_shared_sig_counts(summary_df, PLOT_SHARED)
    plot_summary_upset(
        summary_df,
        status="sig",
        title="Shared significant clusters (summary sheet)",
        output_png=PLOT_UPSET_SIG_PNG,
        output_svg=PLOT_UPSET_SIG_SVG,
    )
    plot_summary_upset(
        summary_df,
        status="unchanged",
        title="Shared unchanged clusters (summary sheet)",
        output_png=PLOT_UPSET_UNCHANGED_PNG,
        output_svg=PLOT_UPSET_UNCHANGED_SVG,
    )

    a_only_summary_df = build_a_only_summary(df_a)
    plot_summary_upset(
        a_only_summary_df,
        status="sig",
        title="Shared significant clusters (GSE115736_GSE116177 only)",
        output_png=PLOT_A_ONLY_UPSET_SIG_PNG,
        output_svg=PLOT_A_ONLY_UPSET_SIG_SVG,
    )
    plot_summary_upset(
        a_only_summary_df,
        status="unchanged",
        title="Shared unchanged clusters (GSE115736_GSE116177 only)",
        output_png=PLOT_A_ONLY_UPSET_UNCHANGED_PNG,
        output_svg=PLOT_A_ONLY_UPSET_UNCHANGED_SVG,
    )

    print(f"\nDone. Output written to:\n  {OUTPUT}")
    print(
        "Plots written to:\n"
        f"  {PLOT_STACKED}\n"
        f"  {PLOT_HEATMAP}\n"
        f"  {PLOT_SHARED}\n"
        f"  {PLOT_UPSET_SIG_PNG}\n"
        f"  {PLOT_UPSET_SIG_SVG}\n"
        f"  {PLOT_UPSET_UNCHANGED_PNG}\n"
        f"  {PLOT_UPSET_UNCHANGED_SVG}\n"
        f"  {PLOT_A_ONLY_UPSET_SIG_PNG}\n"
        f"  {PLOT_A_ONLY_UPSET_SIG_SVG}\n"
        f"  {PLOT_A_ONLY_UPSET_UNCHANGED_PNG}\n"
        f"  {PLOT_A_ONLY_UPSET_UNCHANGED_SVG}"
    )


if __name__ == "__main__":
    main()
