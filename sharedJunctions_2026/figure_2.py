#!/usr/bin/env python3
"""Build Figure 2 as a publication-ready A4 PDF.

Layout:
- Left column (top to bottom):
  1) Human vs Mouse combined stacked panel (from existing SVG/PNG)
  2) HN6 stacked success-percent panel (with Neut/Fibroblast placeholders)
  3) UpSet-like panel for high confidence differentially spliced clusters
  4) UpSet-like panel for high confidence splicing unchanged clusters
- Right column (full page height):
  - Heatmap for clusters that are high confidence in at least 2 cell types

Output:
- /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/figure_2.pdf
"""

from __future__ import annotations

from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
BASE = Path("/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026")
LEAFCUTTER_DIR = BASE / "leafcutter_GSE115736_GSE116177"
HN6_XLSX = LEAFCUTTER_DIR / "unique_sig_clusters_HN6.xlsx"
OUTPUT_PDF = BASE / "figure_2.pdf"
OUTPUT_S1_PDF = BASE / "figrue_S1.pdf"
FILE_A = LEAFCUTTER_DIR / "clusters_sum_table_HN6.txt"
POSSIBLE_AS_CLUSTERS = 1823

PAIR_DIRS: list[tuple[str, str, Path]] = [
    ("Human vs Human", "GSE115736_GSE60424", BASE / "leafcutter_GSE115736_GSE60424" / "sum_table_HN6.txt"),
    ("Human vs Mouse", "GSE115736_GSE116177", BASE / "leafcutter_GSE115736_GSE116177" / "sum_table_HN6.txt"),
    ("Mouse vs Mouse", "GSE116177_GSE180020", BASE / "leafcutter_GSE116177_GSE180020" / "sum_table_HN6.txt"),
]

# -----------------------------------------------------------------------------
# Status labels from filter_unique_sig_clusters_HN6.py
# -----------------------------------------------------------------------------
STATUS_HIGH_CHANGE = "high confidance differentially spliced"
STATUS_LOW_CHANGE = "low confidance differentially spliced"
STATUS_HIGH_CONSERVED = "high confidance splicing unchanged"
STATUS_LOW_CONSERVED = "low confidance splicing unchanged"
STATUS_NOT_INFORMATIVE = "not informative"
STATUS_NOT_SUCCESS = "not success"

THRESHOLD = 0.05
SIG_DELTAPSI_THRESHOLD = 0.1
UNCHANGED_DELTAPSI_THRESHOLD = 0.05

CELL_TYPES_LEFT = ["CD4T", "CD8T", "NveB", "NK", "Mono", "Neut", "Fibroblast"]

# UpSet logic uses the five main HN6 immune cell types exactly as requested.
CELL_TYPES_UPSET = ["CD4T", "CD8T", "NveB", "NK", "Mono"]

# Palette aligned with existing script
COLOR_SIG = "#74c476"          # green
COLOR_UNCH = "#9ecae1"         # light blue
COLOR_NOINFO = "#deebf7"       # very light blue
COLOR_EDGE = "#4d4d4d"
LOW_HATCH = ".."

FONT_FAMILY = "DejaVu Sans"
FONT_SIZE = 7.0


def _configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": FONT_FAMILY,
            "font.size": FONT_SIZE,
            "axes.labelsize": FONT_SIZE,
            "xtick.labelsize": FONT_SIZE,
            "ytick.labelsize": FONT_SIZE,
            "legend.fontsize": FONT_SIZE,
            "axes.titlesize": FONT_SIZE,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 1200,
        }
    )


def _first_gene_label(genes_text: object) -> str:
    txt = str(genes_text).strip()
    if not txt or txt.lower() == "nan":
        return ""
    return re.split(r"[;,|]", txt, maxsplit=1)[0].strip()


def _read_pair_sum_table(sum_table: Path, cell_order: list[str]) -> dict[str, dict[str, float]]:
    """Read one sum_table_HN6.txt into per-cell metrics used by top stacked panels."""
    if not sum_table.exists():
        raise FileNotFoundError(f"Missing top-panel input file: {sum_table}")

    raw = pd.read_csv(sum_table, sep="\t")
    if raw.columns[0] != "Unnamed: 0":
        raw = raw.rename(columns={raw.columns[0]: "metric"})
    else:
        raw = raw.rename(columns={"Unnamed: 0": "metric"})

    metric_col = "metric"
    row_success = raw.loc[raw[metric_col] == "Leafcutter success clusters"]
    row_sig = raw.loc[raw[metric_col].str.contains("sig clusters", na=False)]
    row_unch = raw.loc[raw[metric_col].str.startswith("Unchanged", na=False)]
    row_noinfo = raw.loc[raw[metric_col] == "Not informative"]
    if row_success.empty or row_sig.empty or row_unch.empty or row_noinfo.empty:
        raise ValueError(f"Unexpected sum_table_HN6.txt format in: {sum_table}")

    use_cells = [ct for ct in cell_order if ct in raw.columns]
    out: dict[str, dict[str, float]] = {}
    for ct in use_cells:
        out[ct] = {
            "success": float(row_success.iloc[0][ct] if pd.notna(row_success.iloc[0][ct]) else 0.0),
            "sig_dpsi_01": float(row_sig.iloc[0][ct] if pd.notna(row_sig.iloc[0][ct]) else 0.0),
            "unchanged": float(row_unch.iloc[0][ct] if pd.notna(row_unch.iloc[0][ct]) else 0.0),
            "not_informative": float(row_noinfo.iloc[0][ct] if pd.notna(row_noinfo.iloc[0][ct]) else 0.0),
        }
    return out


def _read_pair_sum_table_all_cells(sum_table: Path) -> dict[str, dict[str, float]]:
    """Read one sum_table_HN6.txt keeping all known cell-type columns in stable order."""
    ordered_cells = ["CD4T", "CD8T", "NveB", "NK", "Mono", "Neut", "Fibroblast"]
    return _read_pair_sum_table(sum_table, ordered_cells)


def _draw_top_three_panels(axes: list[plt.Axes]) -> None:
    """Recreate the HH/HM/MM top panel from merge_leafcutter_results logic."""
    colors = {
        "not_informative": COLOR_NOINFO,
        "unchanged": COLOR_UNCH,
        "sig_dpsi_01": COLOR_SIG,
    }

    panel_cell_orders = [
        CELL_TYPES_UPSET,
        CELL_TYPES_UPSET + ["Neut", "Fibroblast"],
        CELL_TYPES_UPSET,
    ]

    for panel_idx, (ax, (_, _, sum_table)) in enumerate(zip(axes, PAIR_DIRS)):
        summary_by_cell = _read_pair_sum_table(sum_table, panel_cell_orders[panel_idx])
        display_cell_types = [ct for ct in panel_cell_orders[panel_idx] if ct in summary_by_cell]

        success_vals = [summary_by_cell[ct]["success"] for ct in display_cell_types]
        sig_vals = [summary_by_cell[ct]["sig_dpsi_01"] for ct in display_cell_types]
        unchanged_vals = [summary_by_cell[ct]["unchanged"] for ct in display_cell_types]
        noinfo_vals = [summary_by_cell[ct]["not_informative"] for ct in display_cell_types]

        noinfo_pct = [((v / s) * 100.0) if s > 0 else 0.0 for v, s in zip(noinfo_vals, success_vals)]
        unchanged_pct = [((v / s) * 100.0) if s > 0 else 0.0 for v, s in zip(unchanged_vals, success_vals)]
        sig_pct = [((v / s) * 100.0) if s > 0 else 0.0 for v, s in zip(sig_vals, success_vals)]
        bottom_sig = [a + b for a, b in zip(noinfo_pct, unchanged_pct)]

        x = np.arange(len(display_cell_types))
        ax.bar(x, noinfo_pct, color=colors["not_informative"])
        ax.bar(x, unchanged_pct, bottom=noinfo_pct, color=colors["unchanged"])
        ax.bar(x, sig_pct, bottom=bottom_sig, color=colors["sig_dpsi_01"])

        ax.set_xticks(x)
        rotation = 35 if panel_idx == 1 else 0
        ha = "right" if panel_idx == 1 else "center"
        ax.set_xticklabels(display_cell_types, rotation=rotation, ha=ha)
        ax.tick_params(axis="x", pad=1.5)
        ax.margins(x=0.0)
        ax.set_xlim(-0.5, len(display_cell_types) - 0.5)
        ax.set_ylim(0, 110)

    axes[0].set_ylabel("Percent of success clusters (%)")
    for idx, ax in enumerate(axes):
        if idx == 0:
            ax.tick_params(axis="y", which="both", left=True, labelleft=True)
        else:
            ax.tick_params(axis="y", which="both", left=False, labelleft=False)


def _draw_top_three_panels_possible_as(axes: list[plt.Axes]) -> None:
    """Recreate combined_stacked_possibleAS1823_sig_deltapsi01_HH_HM_MM as 3 separate panels."""
    colors = {
        "not_success": "#f0f0f0",
        "not_informative": COLOR_NOINFO,
        "unchanged": COLOR_UNCH,
        "sig_dpsi_01": COLOR_SIG,
    }

    for panel_idx, (ax, (_, _, sum_table)) in enumerate(zip(axes, PAIR_DIRS)):
        summary_by_cell = _read_pair_sum_table_all_cells(sum_table)
        display_cell_types = [ct for ct in ["CD4T", "CD8T", "NveB", "NK", "Mono", "Neut", "Fibroblast"] if ct in summary_by_cell]

        success_vals = [summary_by_cell[ct]["success"] for ct in display_cell_types]
        sig_vals = [summary_by_cell[ct]["sig_dpsi_01"] for ct in display_cell_types]
        unchanged_vals = [summary_by_cell[ct]["unchanged"] for ct in display_cell_types]
        noinfo_vals = [summary_by_cell[ct]["not_informative"] for ct in display_cell_types]
        not_success_vals = [max(0.0, float(POSSIBLE_AS_CLUSTERS) - s) for s in success_vals]

        denom = float(POSSIBLE_AS_CLUSTERS)
        not_success_pct = [(v / denom) * 100.0 for v in not_success_vals]
        noinfo_pct = [(v / denom) * 100.0 for v in noinfo_vals]
        unchanged_pct = [(v / denom) * 100.0 for v in unchanged_vals]
        sig_pct = [(v / denom) * 100.0 for v in sig_vals]

        b2 = not_success_pct
        b3 = [a + b for a, b in zip(not_success_pct, noinfo_pct)]
        b4 = [a + b for a, b in zip(b3, unchanged_pct)]

        x = np.arange(len(display_cell_types))
        ax.bar(x, not_success_pct, color=colors["not_success"])
        ax.bar(x, noinfo_pct, bottom=b2, color=colors["not_informative"])
        ax.bar(x, unchanged_pct, bottom=b3, color=colors["unchanged"])
        ax.bar(x, sig_pct, bottom=b4, color=colors["sig_dpsi_01"])

        ax.set_xticks(x)
        rotation = 35 if panel_idx == 1 else 0
        ha = "right" if panel_idx == 1 else "center"
        ax.set_xticklabels(display_cell_types, rotation=rotation, ha=ha)
        ax.tick_params(axis="x", pad=1.5)
        ax.margins(x=0.0)
        ax.set_xlim(-0.5, len(display_cell_types) - 0.5)
        ax.set_ylim(0, 110)

    axes[0].set_ylabel(f"Percent of possible AS clusters (n={POSSIBLE_AS_CLUSTERS})")
    for idx, ax in enumerate(axes):
        if idx == 0:
            ax.tick_params(axis="y", which="both", left=True, labelleft=True)
        else:
            ax.tick_params(axis="y", which="both", left=False, labelleft=False)


def _read_full_status_sheet() -> pd.DataFrame:
    if not HN6_XLSX.exists():
        raise FileNotFoundError(f"Missing input Excel file: {HN6_XLSX}")

    df = pd.read_excel(HN6_XLSX, sheet_name="all_cluster_status")
    if "cluster" not in df.columns:
        raise ValueError("Sheet 'all_cluster_status' must include a 'cluster' column")
    if "genes" not in df.columns:
        df["genes"] = ""

    df["cluster"] = df["cluster"].astype(str)
    for ct in CELL_TYPES_LEFT:
        if ct not in df.columns:
            df[ct] = ""
        df[ct] = df[ct].fillna("").astype(str)
    return df


def _compute_stacked_pct_arrays(full_status_df: pd.DataFrame) -> tuple[list[str], dict[str, np.ndarray], np.ndarray]:
    """Compute stacked success-percent bars.

    Neut and Fibroblast are kept as placeholders with invisible (zero-height) bars.
    """
    labels = ["CD4T", "CD8T", "NveB", "NK", "Mono", "Neut", "Fibroblasts"]

    n = len(CELL_TYPES_LEFT)
    v_high_unch = np.zeros(n, dtype=float)
    v_low_unch = np.zeros(n, dtype=float)
    v_high_sig = np.zeros(n, dtype=float)
    v_low_sig = np.zeros(n, dtype=float)
    v_noinfo = np.zeros(n, dtype=float)
    success_counts = np.zeros(n, dtype=float)

    for i, ct in enumerate(CELL_TYPES_LEFT):
        # Requested placeholders for Neut and Fibroblast: keep category slots with no visible bars.
        if ct in {"Neut", "Fibroblast"}:
            continue

        s = full_status_df[ct]
        success = (s != STATUS_NOT_SUCCESS).sum()
        success_counts[i] = float(success)
        if success == 0:
            continue

        high_sig = (s == STATUS_HIGH_CHANGE).sum()
        low_sig = (s == STATUS_LOW_CHANGE).sum()
        high_unch = (s == STATUS_HIGH_CONSERVED).sum()
        low_unch = (s == STATUS_LOW_CONSERVED).sum()
        noinfo = (s == STATUS_NOT_INFORMATIVE).sum()

        v_high_unch[i] = (high_unch / success) * 100.0
        v_low_unch[i] = (low_unch / success) * 100.0
        v_high_sig[i] = (high_sig / success) * 100.0
        v_low_sig[i] = (low_sig / success) * 100.0
        v_noinfo[i] = (noinfo / success) * 100.0

    arrays = {
        "high_unch": v_high_unch,
        "low_unch": v_low_unch,
        "high_sig": v_high_sig,
        "low_sig": v_low_sig,
        "noinfo": v_noinfo,
    }
    return labels, arrays, success_counts


def _exclusive_intersections(df: pd.DataFrame, status_label: str, ct_cols: list[str]) -> dict[frozenset[str], int]:
    intersections: dict[frozenset[str], int] = {}
    for _, row in df.iterrows():
        members = frozenset(ct for ct in ct_cols if str(row.get(ct, "")) == status_label)
        if members:
            intersections[members] = intersections.get(members, 0) + 1
    return intersections


def _draw_upset_on_subspec(
    fig: plt.Figure,
    outer_subspec,
    full_status_df: pd.DataFrame,
    status_label: str,
    max_combinations: int = 36,
) -> None:
    ct_cols = [ct for ct in CELL_TYPES_UPSET if ct in full_status_df.columns]
    intersections = _exclusive_intersections(full_status_df, status_label, ct_cols)

    sub = outer_subspec.subgridspec(2, 1, height_ratios=[3, 2], hspace=0.06)
    ax_bar = fig.add_subplot(sub[0])
    ax_mat = fig.add_subplot(sub[1], sharex=ax_bar)

    if not intersections:
        ax_bar.text(
            0.5,
            0.5,
            "No combinations",
            ha="center",
            va="center",
            transform=ax_bar.transAxes,
            fontsize=FONT_SIZE,
        )
        ax_bar.set_axis_off()
        ax_mat.set_axis_off()
        return

    ordered = sorted(
        intersections.items(),
        key=lambda item: (-item[1], -len(item[0]), tuple(sorted(item[0]))),
    )[:max_combinations]

    combos = [combo for combo, _ in ordered]
    counts = [count for _, count in ordered]
    x = np.arange(len(combos))

    ax_bar.bar(x, counts, color="#2b8cbe")
    ax_bar.set_ylabel("Clusters")
    ax_bar.set_ylim(0, max(counts) * 1.18)
    ax_bar.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.45)
    ax_bar.tick_params(axis="x", labelbottom=False)

    pad = max(1, int(0.01 * max(counts)))
    for i, c in enumerate(counts):
        ax_bar.text(i, c + pad, str(c), ha="center", va="bottom", fontsize=FONT_SIZE)

    y_map = {ct: y for y, ct in enumerate(reversed(ct_cols))}
    for i, combo in enumerate(combos):
        ys = sorted(y_map[ct] for ct in combo)
        for y in y_map.values():
            ax_mat.plot(i, y, "o", color="#d9d9d9", markersize=2.2)
        for y in ys:
            ax_mat.plot(i, y, "o", color=COLOR_SIG if status_label == STATUS_HIGH_CHANGE else COLOR_UNCH, markersize=2.8)
        if len(ys) >= 2:
            ax_mat.plot(
                [i, i],
                [ys[0], ys[-1]],
                color=COLOR_SIG if status_label == STATUS_HIGH_CHANGE else COLOR_UNCH,
                linewidth=1.3,
            )

    ax_mat.set_yticks([y_map[ct] for ct in reversed(ct_cols)])
    ax_mat.set_yticklabels(list(reversed(ct_cols)), fontsize=FONT_SIZE)
    ax_mat.set_ylim(-0.7, len(ct_cols) - 0.3)
    ax_mat.set_xlabel("Combinations (sorted by size)")
    ax_mat.set_xticklabels([])
    ax_mat.grid(axis="x", linestyle=":", linewidth=0.4, alpha=0.35)


def _draw_heatmap(ax: plt.Axes, full_status_df: pd.DataFrame) -> int:
    """Side heatmap for clusters high confidence in at least 2 cell types."""
    ct_cols = [ct for ct in CELL_TYPES_UPSET if ct in full_status_df.columns]

    high_mask = full_status_df[ct_cols].isin({STATUS_HIGH_CHANGE, STATUS_HIGH_CONSERVED})
    keep = high_mask.sum(axis=1) >= 2
    df_h = full_status_df.loc[keep, ["cluster", "genes"] + ct_cols].copy()

    if df_h.empty:
        ax.text(
            0.5,
            0.5,
            "No clusters with high confidence in >= 2 cell types",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=FONT_SIZE,
        )
        ax.set_axis_off()
        return 0

    # Sort by number of high-confidence statuses, then by pattern for readability.
    df_h["__n_high"] = df_h[ct_cols].isin({STATUS_HIGH_CHANGE, STATUS_HIGH_CONSERVED}).sum(axis=1)
    for ct in ct_cols:
        df_h[f"__sort_{ct}"] = df_h[ct].map(
            {
                STATUS_HIGH_CHANGE: 0,
                STATUS_HIGH_CONSERVED: 1,
            }
        ).fillna(2)

    sort_cols = ["__n_high"] + [f"__sort_{ct}" for ct in ct_cols] + ["genes", "cluster"]
    ascending = [False] + [True] * (len(sort_cols) - 1)
    df_h = df_h.sort_values(sort_cols, ascending=ascending, kind="stable").reset_index(drop=True)

    mat = np.full((len(df_h), len(ct_cols)), 2.0, dtype=float)
    for j, ct in enumerate(ct_cols):
        col = df_h[ct]
        mat[:, j] = np.where(col == STATUS_HIGH_CHANGE, 1.0, np.where(col == STATUS_HIGH_CONSERVED, 0.0, 2.0))

    cmap = ListedColormap(["#9ecae1", "#74c476", COLOR_NOINFO])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

    ax.pcolormesh(mat, cmap=cmap, norm=norm)
    ax.set_xticks(np.arange(len(ct_cols)) + 0.5)
    ax.set_xticklabels(ct_cols, fontsize=FONT_SIZE, rotation=45, ha="left")
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")

    y_font = max(2.5, FONT_SIZE - 2.5)
    row_labels_full = [_first_gene_label(v) for v in df_h["genes"]]
    row_labels = [lbl if i % 5 == 0 else "" for i, lbl in enumerate(row_labels_full)]
    ax.set_yticks(np.arange(len(df_h)) + 0.5)
    ax.set_yticklabels(row_labels, fontsize=y_font)
    ax.invert_yaxis()
    ax.set_ylabel("")

    return len(df_h)


def _draw_middle_panel(ax: plt.Axes, full_status_df: pd.DataFrame) -> None:
    labels, arr, success_counts = _compute_stacked_pct_arrays(full_status_df)
    x = np.arange(len(labels))

    # Keep ordering consistent with the top panel: not informative at the bottom.
    b2 = arr["noinfo"]
    b3 = b2 + arr["high_unch"]
    b4 = b3 + arr["low_unch"]
    b5 = b4 + arr["high_sig"]

    bar_edge = "#4d4d4d"
    ax.bar(x, arr["noinfo"], color=COLOR_NOINFO, edgecolor=bar_edge, linewidth=0.5, label="not informative")
    ax.bar(x, arr["high_unch"], bottom=b2, color=COLOR_UNCH, edgecolor=bar_edge, linewidth=0.5, label="high confidence splicing unchanged")
    ax.bar(x, arr["low_unch"], bottom=b3, color=COLOR_UNCH, hatch=LOW_HATCH, edgecolor=COLOR_EDGE, linewidth=0.5,
        label="low confidence splicing unchanged")
    ax.bar(x, arr["high_sig"], bottom=b4, color=COLOR_SIG, edgecolor=bar_edge, linewidth=0.5, label="high confidence differentially spliced")
    ax.bar(x, arr["low_sig"], bottom=b5, color=COLOR_SIG, hatch=LOW_HATCH, edgecolor=COLOR_EDGE, linewidth=0.5,
        label="low confidence differentially spliced")

    ax.set_xticks(x)
    labels_display = labels[:5] + ["", ""]
    ax.set_xticklabels(labels_display, fontsize=FONT_SIZE, rotation=35, ha="right")
    ax.tick_params(axis="x", pad=1.5)
    # Tight horizontal margins keep category centers visually aligned to the
    # middle segment above in most PDF viewers.
    ax.margins(x=0.0)
    ax.set_xlim(-0.5, len(labels) - 0.5)
    ax.set_ylim(0, 110)
    ax.set_ylabel("Percent of success clusters (%)")

    # Keep a single legend in the whole page: this panel only.
    ax.legend(frameon=False, fontsize=FONT_SIZE, loc="center left", bbox_to_anchor=(1.02, 0.5), ncol=1)


def _draw_middle_panel_all_clusters(ax: plt.Axes, full_status_df: pd.DataFrame) -> None:
    """Recreate unique_sig_clusters_HN6_counts_stacked categories with all clusters denominator."""
    cts = ["CD4T", "CD8T", "NveB", "NK", "Mono", "Neut", "Fibroblast"]
    x = np.arange(len(cts))

    v_not_success = np.zeros(len(cts), dtype=float)
    v_high_sig = np.zeros(len(cts), dtype=float)
    v_low_sig = np.zeros(len(cts), dtype=float)
    v_high_unch = np.zeros(len(cts), dtype=float)
    v_low_unch = np.zeros(len(cts), dtype=float)
    v_noinfo = np.zeros(len(cts), dtype=float)

    for i, ct in enumerate(cts):
        # Keep blank slots for Neut/Fibroblast in this panel.
        if ct in {"Neut", "Fibroblast"}:
            continue
        s = full_status_df[ct]
        v_not_success[i] = float((s == STATUS_NOT_SUCCESS).sum())
        v_high_sig[i] = float((s == STATUS_HIGH_CHANGE).sum())
        v_low_sig[i] = float((s == STATUS_LOW_CHANGE).sum())
        v_high_unch[i] = float((s == STATUS_HIGH_CONSERVED).sum())
        v_low_unch[i] = float((s == STATUS_LOW_CONSERVED).sum())
        v_noinfo[i] = float((s == STATUS_NOT_INFORMATIVE).sum())

    b2 = v_not_success
    b3 = b2 + v_high_sig
    b4 = b3 + v_low_sig
    b5 = b4 + v_high_unch
    b6 = b5 + v_low_unch

    bar_edge = "#4d4d4d"
    ax.bar(x, v_not_success, color="#f0f0f0", edgecolor=bar_edge, linewidth=0.5)
    ax.bar(x, v_high_sig, bottom=b2, color=COLOR_SIG, edgecolor=bar_edge, linewidth=0.5)
    ax.bar(x, v_low_sig, bottom=b3, color=COLOR_SIG, hatch=LOW_HATCH, edgecolor=COLOR_EDGE, linewidth=0.5)
    ax.bar(x, v_high_unch, bottom=b4, color=COLOR_UNCH, edgecolor=bar_edge, linewidth=0.5)
    ax.bar(x, v_low_unch, bottom=b5, color=COLOR_UNCH, hatch=LOW_HATCH, edgecolor=COLOR_EDGE, linewidth=0.5)
    ax.bar(x, v_noinfo, bottom=b6, color=COLOR_NOINFO, edgecolor="#7f7f7f", linewidth=0.5)

    ax.set_xticks(x)
    labels_display = cts[:5] + ["", ""]
    ax.set_xticklabels(labels_display, rotation=35, ha="right")
    ax.tick_params(axis="x", pad=1.5)
    ax.margins(x=0.0)
    ax.set_xlim(-0.5, len(cts) - 0.5)
    ax.set_ylim(0, max(POSSIBLE_AS_CLUSTERS, 1) * 1.02)
    ax.set_ylabel("Clusters")

    legend_handles = [
        Patch(facecolor="#f0f0f0", edgecolor=bar_edge, label="not success"),
        Patch(facecolor=COLOR_NOINFO, edgecolor="#7f7f7f", label="not informative"),
        Patch(facecolor=COLOR_UNCH, edgecolor=bar_edge, label="high confidence splicing unchanged"),
        Patch(facecolor=COLOR_UNCH, edgecolor=COLOR_EDGE, hatch=LOW_HATCH, label="low confidence splicing unchanged"),
        Patch(facecolor=COLOR_SIG, edgecolor=bar_edge, label="high confidence differentially spliced"),
        Patch(facecolor=COLOR_SIG, edgecolor=COLOR_EDGE, hatch=LOW_HATCH, label="low confidence differentially spliced"),
    ]
    ax.legend(handles=legend_handles, frameon=False, fontsize=FONT_SIZE, loc="center left", bbox_to_anchor=(1.02, 0.5), ncol=1)


def _cluster_min_padj_lookup(df: pd.DataFrame, prefix: str) -> pd.Series:
    col = f"{prefix}_p.adjust"
    if col not in df.columns:
        return pd.Series(dtype=float)
    out = df.groupby("cluster", dropna=False)[col].min().astype(float)
    out.index = out.index.astype(str)
    return out


def _cluster_max_abs_dpsi_lookup(df: pd.DataFrame, prefix: str) -> pd.Series:
    col = f"{prefix}_abs_deltapsi"
    if col not in df.columns:
        return pd.Series(dtype=float)
    out = df[col].astype(float).abs().groupby(df["cluster"], dropna=False).max().astype(float)
    out.index = out.index.astype(str)
    return out


def _classify_cluster_status(min_padj: float | None, max_abs_dpsi: float | None) -> str:
    if min_padj is None or max_abs_dpsi is None:
        return ""
    if (min_padj <= THRESHOLD) and (max_abs_dpsi >= SIG_DELTAPSI_THRESHOLD):
        return "sig"
    if (max_abs_dpsi < UNCHANGED_DELTAPSI_THRESHOLD) or (min_padj > THRESHOLD):
        return "unchanged"
    return ""


def _build_a_only_summary_df() -> pd.DataFrame:
    """Build cluster x cell-type status table using only GSE115736_GSE116177 results."""
    df_a = pd.read_csv(FILE_A, sep="\t")
    if "cluster" not in df_a.columns:
        raise ValueError(f"Missing 'cluster' in {FILE_A}")

    cluster_genes = (
        df_a.drop_duplicates(subset="cluster")[["cluster", "genes"]]
        .set_index("cluster")["genes"]
    )

    cluster_status: dict[str, dict[str, str]] = {}
    for ct in CELL_TYPES_UPSET:
        min_padj = _cluster_min_padj_lookup(df_a, ct)
        max_dpsi = _cluster_max_abs_dpsi_lookup(df_a, ct)
        for cl in set(min_padj.index).union(max_dpsi.index):
            status = _classify_cluster_status(min_padj.get(cl), max_dpsi.get(cl))
            if status:
                cluster_status.setdefault(str(cl), {})[ct] = status

    rows: list[dict[str, str]] = []
    for cl in df_a["cluster"].drop_duplicates().astype(str):
        row: dict[str, str] = {"cluster": cl, "genes": str(cluster_genes.get(cl, ""))}
        has_any = False
        for ct in CELL_TYPES_UPSET:
            row[ct] = cluster_status.get(cl, {}).get(ct, "")
            has_any = has_any or bool(row[ct])
        if has_any:
            rows.append(row)
    return pd.DataFrame(rows, columns=["cluster", "genes"] + CELL_TYPES_UPSET)


def _draw_upset_from_summary_status(
    fig: plt.Figure,
    outer_subspec,
    summary_df: pd.DataFrame,
    status_label: str,
    max_combinations: int = 36,
) -> None:
    ct_cols = [ct for ct in CELL_TYPES_UPSET if ct in summary_df.columns]
    intersections = _exclusive_intersections(summary_df, status_label, ct_cols)

    sub = outer_subspec.subgridspec(2, 1, height_ratios=[3, 2], hspace=0.06)
    ax_bar = fig.add_subplot(sub[0])
    ax_mat = fig.add_subplot(sub[1], sharex=ax_bar)

    if not intersections:
        ax_bar.text(0.5, 0.5, "No combinations", ha="center", va="center", transform=ax_bar.transAxes, fontsize=FONT_SIZE)
        ax_bar.set_axis_off()
        ax_mat.set_axis_off()
        return

    ordered = sorted(intersections.items(), key=lambda item: (-item[1], -len(item[0]), tuple(sorted(item[0]))))[:max_combinations]
    combos = [combo for combo, _ in ordered]
    counts = [count for _, count in ordered]
    x = np.arange(len(combos))

    ax_bar.bar(x, counts, color="#2b8cbe")
    ax_bar.set_ylabel("Clusters")
    ax_bar.set_ylim(0, max(counts) * 1.18)
    ax_bar.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.45)
    ax_bar.tick_params(axis="x", labelbottom=False)

    pad = max(1, int(0.01 * max(counts)))
    for i, c in enumerate(counts):
        ax_bar.text(i, c + pad, str(c), ha="center", va="bottom", fontsize=FONT_SIZE)

    y_map = {ct: y for y, ct in enumerate(reversed(ct_cols))}
    for i, combo in enumerate(combos):
        ys = sorted(y_map[ct] for ct in combo)
        for y in y_map.values():
            ax_mat.plot(i, y, "o", color="#d9d9d9", markersize=2.2)
        for y in ys:
            ax_mat.plot(i, y, "o", color=COLOR_SIG if status_label == "sig" else COLOR_UNCH, markersize=2.8)
        if len(ys) >= 2:
            ax_mat.plot([i, i], [ys[0], ys[-1]], color=COLOR_SIG if status_label == "sig" else COLOR_UNCH, linewidth=1.3)

    ax_mat.set_yticks([y_map[ct] for ct in reversed(ct_cols)])
    ax_mat.set_yticklabels(list(reversed(ct_cols)), fontsize=FONT_SIZE)
    ax_mat.set_ylim(-0.7, len(ct_cols) - 0.3)
    ax_mat.set_xlabel("Combinations (sorted by size)")
    ax_mat.set_xticklabels([])
    ax_mat.grid(axis="x", linestyle=":", linewidth=0.4, alpha=0.35)


def build_figure_2() -> None:
    _configure_style()
    full_status_df = _read_full_status_sheet()

    plt.close("all")
    fig = plt.figure(figsize=(8.27, 11.69), constrained_layout=False)
    outer = GridSpec(
        4,
        2,
        figure=fig,
        width_ratios=[2.80, 0.50],
        height_ratios=[1.15, 1.15, 1.15, 1.15],
        wspace=0.38,
        hspace=0.34,
    )

    # Subgrid on the left lets us align the second plot under the center third
    # of the top 3-part panel (HH | HM | MM).
    left = outer[:, 0].subgridspec(
        4,
        3,
        height_ratios=[1.15, 1.15, 1.15, 1.15],
        wspace=0.14,
        hspace=0.34,
    )

    # Top row: recreate HH/HM/MM as 3 separate panels.
    ax_top_hh = fig.add_subplot(left[0, 0])
    ax_top_hm = fig.add_subplot(left[0, 1], sharey=ax_top_hh)
    ax_top_mm = fig.add_subplot(left[0, 2], sharey=ax_top_hh)
    _draw_top_three_panels([ax_top_hh, ax_top_hm, ax_top_mm])

    # Middle stacked percent panel aligned below the middle top segment (HM).
    ax_mid = fig.add_subplot(left[1, 1])
    _draw_middle_panel(ax_mid, full_status_df)

    # UpSet panels (high confidence only)
    _draw_upset_on_subspec(
        fig,
        left[2, :],
        full_status_df,
        status_label=STATUS_HIGH_CHANGE,
    )
    _draw_upset_on_subspec(
        fig,
        left[3, :],
        full_status_df,
        status_label=STATUS_HIGH_CONSERVED,
    )

    # Side heatmap: top 3 rows only so legend has room below
    ax_heat = fig.add_subplot(outer[:3, 1])
    n_keep = _draw_heatmap(ax_heat, full_status_df)

    fig.savefig(OUTPUT_PDF, format="pdf", dpi=1200, bbox_inches="tight")
    plt.close("all")

    print(f"Saved: {OUTPUT_PDF}")


def build_figure_s1() -> None:
    """Build supplemental figure S1 using all-cluster panels and A-only upsets."""
    _configure_style()
    full_status_df = _read_full_status_sheet()
    a_only_summary_df = _build_a_only_summary_df()

    plt.close("all")
    fig = plt.figure(figsize=(8.27, 11.69), constrained_layout=False)
    grid = GridSpec(
        4,
        3,
        figure=fig,
        height_ratios=[1.15, 1.15, 1.15, 1.15],
        wspace=0.14,
        hspace=0.34,
    )

    ax_top_hh = fig.add_subplot(grid[0, 0])
    ax_top_hm = fig.add_subplot(grid[0, 1], sharey=ax_top_hh)
    ax_top_mm = fig.add_subplot(grid[0, 2], sharey=ax_top_hh)
    _draw_top_three_panels_possible_as([ax_top_hh, ax_top_hm, ax_top_mm])

    ax_mid = fig.add_subplot(grid[1, 1])
    _draw_middle_panel_all_clusters(ax_mid, full_status_df)

    _draw_upset_from_summary_status(
        fig,
        grid[2, :],
        a_only_summary_df,
        status_label="sig",
    )
    _draw_upset_from_summary_status(
        fig,
        grid[3, :],
        a_only_summary_df,
        status_label="unchanged",
    )

    fig.savefig(OUTPUT_S1_PDF, format="pdf", dpi=1200, bbox_inches="tight")
    plt.close("all")
    print(f"Saved: {OUTPUT_S1_PDF}")


def main() -> None:
    build_figure_2()
    build_figure_s1()


if __name__ == "__main__":
    main()
