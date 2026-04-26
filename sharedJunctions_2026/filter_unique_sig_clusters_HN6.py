#!/usr/bin/env python3
"""
For each cell type (CD4T, CD8T, NveB, NK, Mono), collect rows from
leafcutter_GSE115736_GSE116177/clusters_sum_table_HN6.txt that satisfy:

    Condition 1 (unique-sig):
            sig   in GSE115736_GSE116177  (p.adjust < 0.05 and abs_deltapsi > 0.1)
                AND unchanged in GSE115736_GSE60424  (abs_deltapsi < 0.05)
                AND unchanged in GSE116177_GSE180020 (abs_deltapsi < 0.05)

  Condition 2 (unchanged everywhere):
      unchanged in all three comparisons (abs_deltapsi < 0.05)

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

def cluster_padj_lookup(df: pd.DataFrame, prefix: str) -> pd.Series:
    """Return per-cluster p.adjust (index = cluster id, float)."""
    col = padj_col(prefix)
    if col not in df.columns:
        return pd.Series(dtype=float)
    first_per_cluster = df.drop_duplicates(subset="cluster").set_index("cluster")[col]
    return first_per_cluster.astype(float)


def cluster_abs_deltapsi_lookup(df: pd.DataFrame, prefix: str) -> pd.Series:
    """Return per-cluster abs_deltapsi (index = cluster id, float)."""
    col = abs_deltapsi_col(prefix)
    if col not in df.columns:
        return pd.Series(dtype=float)
    first_per_cluster = df.drop_duplicates(subset="cluster").set_index("cluster")[col]
    return first_per_cluster.astype(float).abs()


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
    """Plot stacked cluster counts: unique-sig vs unchanged-all for each cell type."""
    ct_order = list(CELL_TYPES.keys())
    unique_vals = [counts_by_ct.get(ct, {}).get("unique_sig", 0) for ct in ct_order]
    unchanged_vals = [counts_by_ct.get(ct, {}).get("unchanged_all", 0) for ct in ct_order]
    totals = [u + n for u, n in zip(unique_vals, unchanged_vals)]

    x = np.arange(len(ct_order))
    plt.close("all")
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.bar(x, unique_vals, color="#1f77b4", label="unique-sig")
    ax.bar(x, unchanged_vals, bottom=unique_vals, color="#bdbdbd", label="unchanged-all")

    for i, total in enumerate(totals):
        ax.text(i, total + max(5, int(0.01 * max(totals, default=1))), str(total),
                ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(ct_order)
    ax.set_ylabel("Clusters")
    ax.set_title("Selected clusters per cell type")
    ax.legend(frameon=False)
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

    # Pre-build per-cluster p.adjust and abs_deltapsi lookups from B and C
    padj_b: dict[str, pd.Series] = {}
    padj_c: dict[str, pd.Series] = {}
    dpsi_b: dict[str, pd.Series] = {}
    dpsi_c: dict[str, pd.Series] = {}
    for ct, (pfx_a, pfx_b, pfx_c) in CELL_TYPES.items():
        padj_b[ct] = cluster_padj_lookup(df_b, pfx_b)
        padj_c[ct] = cluster_padj_lookup(df_c, pfx_c)
        dpsi_b[ct] = cluster_abs_deltapsi_lookup(df_b, pfx_b)
        dpsi_c[ct] = cluster_abs_deltapsi_lookup(df_c, pfx_c)

    # Ordered list of clusters seen across all cell-type sheets (deduplicated)
    seen_clusters: dict[str, None] = {}  # ordered set via dict
    counts_by_ct: dict[str, dict[str, int]] = {}
    # tracks which clusters qualified for each cell type
    qualified_per_ct: dict[str, set[str]] = {ct: set() for ct in CELL_TYPES}

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

            pa = df_a[col_a].astype(float)
            da = df_a[dpsi_col_a].astype(float).abs()
            pb = df_a["cluster"].map(padj_b[ct])
            pc = df_a["cluster"].map(padj_c[ct])
            db = df_a["cluster"].map(dpsi_b[ct])
            dc = df_a["cluster"].map(dpsi_c[ct])

            # Exclude any cluster absent in B or C for this cell type.
            valid_in_all = pb.notna() & pc.notna() & db.notna() & dc.notna()

            sig_a = (pa < THRESHOLD) & (da > SIG_DELTAPSI_THRESHOLD)
            unchanged_a = da < UNCHANGED_DELTAPSI_THRESHOLD
            unchanged_b = db < UNCHANGED_DELTAPSI_THRESHOLD
            unchanged_c = dc < UNCHANGED_DELTAPSI_THRESHOLD

            unique_sig_mask = valid_in_all & sig_a & unchanged_b & unchanged_c
            unchanged_all_mask = valid_in_all & unchanged_a & unchanged_b & unchanged_c
            mask = unique_sig_mask | unchanged_all_mask

            # Keep only this cell type columns from each source file.
            keys = [
                c
                for c in ["cluster", "h_junction", "m_junction", "genes", "rank_h", "rank_m"]
                if c in df_a.columns
            ]
            cols_a = keys + cell_type_cols(df_a, pfx_a)
            cols_b = [c for c in ["cluster", "h_junction", "m_junction"] if c in df_b.columns] + cell_type_cols(df_b, pfx_b)
            cols_c = [c for c in ["cluster", "h_junction", "m_junction"] if c in df_c.columns] + cell_type_cols(df_c, pfx_c)

            base = df_a.loc[mask, cols_a].copy()

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

            unique_sig_clusters = set(df_a.loc[unique_sig_mask, "cluster"].astype(str))
            unchanged_all_clusters = set(df_a.loc[unchanged_all_mask, "cluster"].astype(str))

            n_unique_sig = len(unique_sig_clusters)
            n_unchanged_all = len(unchanged_all_clusters)
            counts_by_ct[ct] = {
                "unique_sig": n_unique_sig,
                "unchanged_all": n_unchanged_all,
            }
            print(f"  [{ct}] {len(result):,} rows written  "
                  f"(unique-sig: {n_unique_sig}, unchanged-all: {n_unchanged_all})")

            for cl in result["cluster"].unique():
                seen_clusters[cl] = None
                qualified_per_ct[ct].add(cl)

        # ------------------------------------------------------------------ #
        # Summary sheet                                                        #
        # ------------------------------------------------------------------ #
        print(f"\nBuilding summary sheet for {len(seen_clusters):,} unique clusters …")

        cluster_genes = (
            df_a.drop_duplicates(subset="cluster")[ ["cluster", "genes"] ]
            .set_index("cluster")["genes"]
        )

        # Build a cluster -> {ct: sig/unchanged status in A} table
        # Use first occurrence of each cluster in file A for p.adjust/abs_deltapsi lookup
        cluster_padj_a: dict[str, dict[str, float]] = {}
        cluster_dpsi_a: dict[str, dict[str, float]] = {}
        for ct, (pfx_a, _, _) in CELL_TYPES.items():
            col_a = padj_col(pfx_a)
            dpsi_col_a = abs_deltapsi_col(pfx_a)
            if col_a not in df_a.columns or dpsi_col_a not in df_a.columns:
                continue
            for cl, val in cluster_padj_lookup(df_a, pfx_a).items():
                cluster_padj_a.setdefault(cl, {})[ct] = val
            for cl, val in cluster_abs_deltapsi_lookup(df_a, pfx_a).items():
                cluster_dpsi_a.setdefault(cl, {})[ct] = val

        rows = []
        for cl in seen_clusters:
            row: dict[str, str] = {"cluster": cl, "genes": str(cluster_genes.get(cl, ""))}
            for ct in CELL_TYPES:
                if cl not in qualified_per_ct[ct]:
                    row[ct] = ""  # did not qualify for this cell type
                else:
                    padj_val = cluster_padj_a.get(cl, {}).get(ct)
                    dpsi_val = cluster_dpsi_a.get(cl, {}).get(ct)
                    if padj_val is None or dpsi_val is None:
                        row[ct] = ""
                    else:
                        if (padj_val < THRESHOLD) and (dpsi_val > SIG_DELTAPSI_THRESHOLD):
                            row[ct] = "sig"
                        elif dpsi_val < UNCHANGED_DELTAPSI_THRESHOLD:
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

    print(f"\nDone. Output written to:\n  {OUTPUT}")
    print(f"Plots written to:\n  {PLOT_STACKED}\n  {PLOT_HEATMAP}\n  {PLOT_SHARED}")


if __name__ == "__main__":
    main()
