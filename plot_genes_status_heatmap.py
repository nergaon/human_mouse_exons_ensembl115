#!/usr/bin/env python3
"""Create a heatmap from genes_status table.

Rows = genes, columns = cell types.
Colors match the other plots in this project:
  not success        -> #f0f0f0 (light gray)
  not informative    -> #deebf7 (very light blue)
  splicing conserved -> #9ecae1 (medium blue)  [any confidence level]
  splicing change    -> #74c476 (green)         [any confidence level]
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

INPUT_FILE = Path(
    "/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026"
    "/human_immune_DS/genes_status"
)
OUTPUT_PNG = INPUT_FILE.parent / "genes_status_heatmap.png"
OUTPUT_SVG = INPUT_FILE.parent / "genes_status_heatmap.svg"

# ── colour scheme ────────────────────────────────────────────────────────────
STATUS_COLORS = {
    "not success":        "#f0f0f0",
    "not informative":    "#deebf7",
    "splicing conserved": "#9ecae1",
    "splicing change":    "#74c476",
}


def normalise_status(raw: str) -> str:
    """Collapse all confidence levels to the four canonical categories."""
    v = raw.strip().lower()
    if "splicing change" in v:
        return "splicing change"
    if "splicing conserved" in v:
        return "splicing conserved"
    if "not informative" in v:
        return "not informative"
    return "not success"


def main() -> None:
    df = pd.read_csv(INPUT_FILE, sep="\t", index_col=0)

    # Normalise all cell values
    norm = df.map(normalise_status)

    category_order = ["splicing change", "splicing conserved", "not informative", "not success"]

    # Sort genes: primary key = pattern across the 5 immune cell types of interest,
    # secondary key = full row pattern (remaining columns).
    SORT_COLS = ["CD4T", "CD8T", "NveB", "NK", "Mono"]
    cat_index = {c: i for i, c in enumerate(category_order)}
    sort_key_cols = [c for c in SORT_COLS if c in norm.columns]
    other_cols = [c for c in norm.columns if c not in sort_key_cols]
    sort_df = norm[sort_key_cols + other_cols].map(lambda x: cat_index[x])
    norm = norm.loc[sort_df.sort_values(sort_key_cols + other_cols).index]

    # Build numeric matrix for imshow
    num_matrix = norm.map(lambda x: cat_index[x]).values  # shape (genes, cells)

    cmap_colors = [STATUS_COLORS[c] for c in category_order]
    from matplotlib.colors import ListedColormap, BoundaryNorm
    cmap = ListedColormap(cmap_colors)
    bounds = list(range(len(category_order) + 1))
    bnorm = BoundaryNorm(bounds, cmap.N)

    n_genes, n_cells = num_matrix.shape
    cell_w = 1.0          # inches per column
    gene_h = 0.22         # inches per row
    fig_w = max(8, n_cells * cell_w + 3)
    fig_h = max(6, n_genes * gene_h + 2)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(
        num_matrix,
        cmap=cmap,
        norm=bnorm,
        aspect="auto",
        interpolation="none",
    )

    # Axes labels
    ax.set_xticks(range(n_cells))
    ax.set_xticklabels(norm.columns, rotation=35, ha="right", fontsize=9)
    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(norm.index, fontsize=7)

    ax.set_xlabel("Cell type", fontsize=11)
    ax.set_ylabel("Gene", fontsize=11)
    ax.set_title("Gene splicing status per cell type", fontsize=12, pad=10)

    # Grid lines between cells
    ax.set_xticks(np.arange(-0.5, n_cells, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_genes, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Legend
    legend_labels = {
        "splicing change":    "Splicing change",
        "splicing conserved": "Splicing conserved",
        "not informative":    "Not informative",
        "not success":        "Not success",
    }
    patches = [
        mpatches.Patch(facecolor=STATUS_COLORS[k], edgecolor="grey", linewidth=0.4, label=v)
        for k, v in legend_labels.items()
    ]
    ax.legend(
        handles=patches,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        borderaxespad=0,
        frameon=False,
        fontsize=9,
    )

    plt.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=200, bbox_inches="tight")
    fig.savefig(OUTPUT_SVG, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUTPUT_PNG}")
    print(f"Wrote {OUTPUT_SVG}")


if __name__ == "__main__":
    main()
