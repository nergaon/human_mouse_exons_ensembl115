#!/usr/bin/env python3
"""Plot selected clusters across all cells (Counts top, PSI bottom).

Requirements implemented:
- Plot only selected clusters.
- Top panel: Counts, bottom panel: PSI.
- Title: gene name only (no p-value, no deltapsi).
- Legend intron format: chr:start-end (not chr:start:end).
- Cell-group order on x-axis: CD4, CD8, NveB, NK, Mono, [small gap], Fibroblast (if present).
- Output format: PDF only.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['font.family'] = 'Liberation Sans'
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

BASE = Path("/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026")
VALUE_FILE_MAIN = BASE / "AS_clusters_value_HN6.txt"
VALUE_FILE_FIBRO = BASE / "AS_clusters_value_fibroblast_HN6.txt"
GENE_TABLE_MAIN = BASE / "leafcutter_GSE115736_GSE116177" / "clusters_sum_table_HN6.txt"
GENE_TABLE_FIBRO = BASE / "leafcutter_EMTAB5919H_EMTAB5919M" / "clusters_sum_table_HN6.txt"
OUTPUT_DIR = BASE / "genes_figs_selected_clusters"

TARGET_CLUSTERS = [
    #"chr21:clu_46835",
    #"chr9:clu_81561",
    #"chr9:clu_80597",
    #"chr8:clu_76270",
    #"clur3:clu_57873",
    #"chr6:clu_68620",
    #"chr3:clu_60710",
    #"chr3:clu_60717",
    #"chr14:clu_15672",
    #"chr3:clu_56081",
    #"chr1:clu_39954",
    #"chr20:clu_46642",
    "chr4:clu_62944"
]

CELL_GROUP_ORDER = ["CD4", "CD8", "NveB", "NK", "Mono", "Fibroblast"]
IMMUNE_CELL_GROUPS = ["CD4", "CD8", "NveB", "NK", "Mono"]
CELL_GROUP_DISPLAY = {
    "CD4": "T4",
    "CD8": "T8",
    "NveB": "B",
    "NK": "NK",
    "Mono": "Mo",
    "Fibroblast": "Fibroblasts",
}
DATASET_ORDER = ["GSE180020", "GSE116177", "GSE115736", "GSE60424", "MM", "HS"]
DATASET_LABELS = {
    "GSE180020": "M2",
    "GSE116177": "M1",
    "GSE115736": "H1",
    "GSE60424": "H2",
    "MM": "M",
    "HS": "H",
}


def _cluster_id_from_index(junction_id: str) -> str:
    if ":" not in junction_id:
        return ""
    junc_part, clu_suffix = junction_id.rsplit(":", 1)
    chrom = junc_part.split(":", 1)[0]
    return f"{chrom}:{clu_suffix}"


def _safe_cluster_candidates(cluster_id: str) -> List[str]:
    """Return fallback candidates for typo-prone cluster strings."""
    candidates = [cluster_id]
    if cluster_id.startswith("clur"):
        candidates.append("chr" + cluster_id[len("clur"):])
    return candidates


def _load_value_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df["cluster"] = [
        _cluster_id_from_index(str(idx))
        for idx in df.index.astype(str)
    ]
    return df


def _load_cluster_gene_map(paths: List[Path]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for p in paths:
        if not p.exists():
            continue
        try:
            df = pd.read_csv(p, sep="\t")
        except Exception:
            continue
        if "cluster" not in df.columns or "genes" not in df.columns:
            continue
        for _, row in df[["cluster", "genes"]].drop_duplicates(subset="cluster").iterrows():
            cl = str(row["cluster"])
            gn = str(row["genes"]) if pd.notna(row["genes"]) else ""
            if cl and cl not in out and gn:
                out[cl] = gn
    return out


def _first_gene_label(genes_text: str) -> str:
    txt = str(genes_text).strip()
    if not txt or txt.lower() == "nan":
        return ""
    return re.split(r"[;,|]", txt, maxsplit=1)[0].strip()


def _sample_cell_group(sample_name: str) -> str:
    s = sample_name
    upper = s.upper()

    if upper.startswith("HS") or upper.startswith("MM"):
        return "Fibroblast"
    if s.startswith(("Mono", "InfMono", "Monocytes", "Mo_6C+")):
        return "Mono"
    if s.startswith(("NveB", "BCell", "Bcells", "B_fo")):
        return "NveB"
    if s.startswith("NK"):
        return "NK"
    if s.startswith(("CD8T", "Cd8T", "MemCD8T", "NveCd8T", "CD8_", "T_8_")):
        return "CD8"
    if s.startswith(("CD4T", "EffCD4T", "MemCD4T", "NveCD4T", "CD4_", "T_4_")):
        return "CD4"
    return ""


def _natural_sort_key(text: str) -> List[object]:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", text)]


def extract_dataset_id(sample_name: str) -> str:
    upper_name = sample_name.upper()
    if upper_name.startswith("HS"):
        return "HS"
    if upper_name.startswith("MM"):
        return "MM"

    match = re.search(r"(GSE\d+)$", sample_name)
    return match.group(1) if match else ""


def _dataset_display_label(dataset_id: str) -> str:
    return DATASET_LABELS.get(dataset_id, dataset_id)


def _sort_samples(samples: List[str]) -> List[str]:
    order_map = {ct: i for i, ct in enumerate(CELL_GROUP_ORDER)}

    def key(s: str) -> Tuple[int, str]:
        grp = _sample_cell_group(s)
        ds = extract_dataset_id(s)
        ds_idx = DATASET_ORDER.index(ds) if ds in DATASET_ORDER else 999
        return (order_map.get(grp, 999), ds_idx, _natural_sort_key(s))

    return sorted(samples, key=key)


def _build_positions(
    samples_sorted: List[str],
    within_step: float = 0.30,
    dataset_gap: float = 0.5,
    cell_gap: float = 0.65,
    fibro_gap: float = 2.50,
) -> np.ndarray:
    if not samples_sorted:
        return np.array([])
    pos = np.zeros(len(samples_sorted), dtype=float)
    prev_cell = _sample_cell_group(samples_sorted[0])
    prev_ds = extract_dataset_id(samples_sorted[0])
    for i in range(1, len(samples_sorted)):
        cur_cell = _sample_cell_group(samples_sorted[i])
        cur_ds = extract_dataset_id(samples_sorted[i])
        pos[i] = pos[i - 1] + within_step
        if cur_ds != prev_ds:
            pos[i] += dataset_gap
        if cur_cell != prev_cell:
            # Make fibroblast visually separate from immune cell blocks.
            if "Fibroblast" in {cur_cell, prev_cell}:
                pos[i] += fibro_gap
            else:
                pos[i] += cell_gap
        prev_cell = cur_cell
        prev_ds = cur_ds
    return pos


def _build_xtick_labels(samples_sorted: List[str]) -> List[str]:
    labels = ["" for _ in samples_sorted]
    prev = ""
    for i, s in enumerate(samples_sorted):
        grp = _sample_cell_group(s)
        if grp and grp != prev:
            labels[i] = grp
            prev = grp
    return labels


def _group_segments(samples_sorted: List[str]) -> List[Tuple[str, int, int]]:
    """Return contiguous sample segments as (cell_group, start_idx, end_idx)."""
    if not samples_sorted:
        return []

    segments: List[Tuple[str, int, int]] = []
    start = 0
    prev = _sample_cell_group(samples_sorted[0])
    for i in range(1, len(samples_sorted)):
        cur = _sample_cell_group(samples_sorted[i])
        if cur != prev:
            segments.append((prev, start, i - 1))
            start = i
            prev = cur
    segments.append((prev, start, len(samples_sorted) - 1))
    return segments


def _contiguous_dataset_blocks(samples_sorted: List[str], positions: np.ndarray) -> List[Tuple[str, float, float, float]]:
    """Return contiguous dataset blocks as (dataset_id, start_x, end_x, center_x)."""
    blocks: List[Tuple[str, float, float, float]] = []
    if not samples_sorted:
        return blocks

    start = 0
    prev = extract_dataset_id(samples_sorted[0])
    for i in range(1, len(samples_sorted) + 1):
        end_block = i == len(samples_sorted) or extract_dataset_id(samples_sorted[i]) != prev
        if end_block:
            sx = float(positions[start])
            ex = float(positions[i - 1])
            cx = float((sx + ex) / 2.0)
            blocks.append((prev, sx, ex, cx))
            if i < len(samples_sorted):
                start = i
                prev = extract_dataset_id(samples_sorted[i])
    return blocks


def _has_fibro_data(df_cluster: pd.DataFrame) -> bool:
    """Return True only when fibroblast sample columns have at least one real value."""
    numeric_cols = [c for c in df_cluster.columns if c != "cluster" and pd.api.types.is_numeric_dtype(df_cluster[c])]
    fibro_cols = [c for c in numeric_cols if _sample_cell_group(c) == "Fibroblast"]
    if not fibro_cols:
        return False
    fibro_vals = df_cluster[fibro_cols]
    return bool(fibro_vals.notna().any().any())


def _legend_label_from_junction(junction_id: str) -> str:
    # expected: chr:start:end:clu_x or chr:start:end
    parts = str(junction_id).split(":")
    if len(parts) >= 3:
        return f"{parts[0]}:{parts[1]}-{parts[2]}"
    return str(junction_id)


def _plot_cluster(df_cluster: pd.DataFrame, sample_cols: List[str], gene_title: str, cluster_id: str, out_dir: Path) -> None:
    if df_cluster.empty or not sample_cols:
        return

    df_cluster = df_cluster.copy()
    psi = df_cluster[sample_cols].div(df_cluster[sample_cols].sum(axis=0), axis=1).fillna(0.0)
    counts = df_cluster[sample_cols].copy()

    intron_rows = df_cluster.index.astype(str).tolist()
    legend_labels = [_legend_label_from_junction(j) for j in intron_rows]

    cmap = plt.get_cmap("tab20")
    colors = [cmap(i % cmap.N) for i in range(len(intron_rows))]
    legend_handles = [Patch(facecolor=colors[i], edgecolor="none", label=legend_labels[i]) for i in range(len(intron_rows))]

    samples_sorted = _sort_samples(sample_cols)
    bar_w = 0.19
    # Keep samples from the same dataset contiguous (no extra intra-dataset gap).
    positions = _build_positions(samples_sorted, within_step=bar_w)
    segments = _group_segments(samples_sorted)
    ds_blocks = _contiguous_dataset_blocks(samples_sorted, positions)

    psi_vals = psi[samples_sorted].to_numpy(dtype=float)
    cnt_vals = counts[samples_sorted].to_numpy(dtype=float)

    bottom_cnt = np.zeros(len(samples_sorted), dtype=float)
    bottom_psi = np.zeros(len(samples_sorted), dtype=float)

    plt.close("all")
    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        # A4 width x half-A4 height in inches.
        figsize=(11.69, 4.13),
        gridspec_kw={"height_ratios": [1, 1]},
        sharex=False,
    )

    for i in range(len(intron_rows)):
        c = colors[i]
        ax_top.bar(positions, cnt_vals[i], width=bar_w, bottom=bottom_cnt, color=c, edgecolor="none")
        bottom_cnt += cnt_vals[i]

        ax_bottom.bar(positions, psi_vals[i], width=bar_w, bottom=bottom_psi, color=c, edgecolor="none")
        bottom_psi += psi_vals[i]

    # Keep identical x-limits across panels and mirror boundary spacing for
    # first/last cell groups so Mono/Fibroblast have edge space like inner groups.
    left_pad = bar_w * 0.55
    right_pad = bar_w * 0.55
    if len(segments) > 1:
        first_end = segments[0][2]
        second_start = segments[1][1]
        left_pad = max(left_pad, (positions[second_start] - positions[first_end]) / 2.0)

        prev_end = segments[-2][2]
        last_start = segments[-1][1]
        right_pad = max(right_pad, (positions[last_start] - positions[prev_end]) / 2.0)

    x_min = float(np.min(positions) - left_pad)
    x_max = float(np.max(positions) + right_pad)
    ax_top.set_xlim(x_min, x_max)
    ax_bottom.set_xlim(x_min, x_max)

    # Top: counts
    if segments:
        centers = [float(np.mean(positions[s:e + 1])) for _, s, e in segments]
        group_labels = [CELL_GROUP_DISPLAY.get(grp, grp) for grp, _, _ in segments]
        ax_top.set_xticks(centers)
        ax_top.set_xticklabels(group_labels, fontsize=7, rotation=0, ha="center")
    else:
        ax_top.set_xticks([])
    ax_top.tick_params(axis="x", length=0, labelbottom=True)
    ax_top.set_ylabel("JSR counts")

    # Bottom: PSI
    ax_bottom.set_xticks([center for _, _, _, center in ds_blocks])
    ax_bottom.set_xticklabels([_dataset_display_label(ds) for ds, _, _, _ in ds_blocks], fontsize=8, rotation=0, ha="center")
    ax_bottom.tick_params(axis="x", length=0)
    ax_bottom.set_ylabel("PSI")
    ax_bottom.set_ylim(0, 1)

    # Dashed vertical lines between different datasets (within cell groups).
    for i in range(1, len(samples_sorted)):
        prev_ds = extract_dataset_id(samples_sorted[i - 1])
        cur_ds = extract_dataset_id(samples_sorted[i])
        prev_cell = _sample_cell_group(samples_sorted[i - 1])
        cur_cell = _sample_cell_group(samples_sorted[i])
        if cur_ds != prev_ds and cur_cell == prev_cell:
            x_sep = (positions[i - 1] + positions[i]) / 2.0
            ax_top.axvline(x_sep, color="#aaaaaa", linewidth=0.7, linestyle="--", zorder=0)
            ax_bottom.axvline(x_sep, color="#aaaaaa", linewidth=0.7, linestyle="--", zorder=0)

    # Separate cell-type blocks with vertical guide lines.
    for i, (grp, s, e) in enumerate(segments[:-1]):
        next_grp = segments[i + 1][0]
        x_sep = (positions[e] + positions[e + 1]) / 2.0
        is_fibro_boundary = "Fibroblast" in {grp, next_grp}
        if is_fibro_boundary:
            band_half = bar_w * 1.2
            ax_top.axvspan(x_sep - band_half, x_sep + band_half, color="#f0f0f0", alpha=0.9, zorder=0)
            ax_bottom.axvspan(x_sep - band_half, x_sep + band_half, color="#f0f0f0", alpha=0.9, zorder=0)
        sep_color = "#4d4d4d" if is_fibro_boundary else "#999999"
        sep_lw = 1.8 if is_fibro_boundary else 0.8
        ax_top.axvline(x_sep, color=sep_color, linewidth=sep_lw)
        ax_bottom.axvline(x_sep, color=sep_color, linewidth=sep_lw)

    fig.suptitle(gene_title if gene_title else cluster_id, fontsize=10)
    fig.legend(legend_handles, legend_labels, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.03), ncol=4, frameon=False)
    plt.subplots_adjust(hspace=0.30)

    safe_name = re.sub(r"[^A-Za-z0-9._=+-]+", "_", gene_title if gene_title else cluster_id)
    safe_name = re.sub(r"_+", "_", safe_name)

    out_pdf = out_dir / f"{safe_name}_{cluster_id.replace(':', '_')}.pdf"
    fig.savefig(out_pdf, bbox_inches="tight", dpi=300)
    plt.close("all")


def run() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df_main = _load_value_table(VALUE_FILE_MAIN) if VALUE_FILE_MAIN.exists() else None
    df_fibro = _load_value_table(VALUE_FILE_FIBRO) if VALUE_FILE_FIBRO.exists() else None

    gene_map = _load_cluster_gene_map([GENE_TABLE_MAIN, GENE_TABLE_FIBRO])

    if df_main is None and df_fibro is None:
        raise FileNotFoundError("Missing both AS value files.")

    for requested in TARGET_CLUSTERS:
        candidates = _safe_cluster_candidates(requested)

        parts: List[pd.DataFrame] = []
        found_cluster = ""
        fibro_success = False

        for cid in candidates:
            got_any = False
            if df_main is not None:
                sub = df_main.loc[df_main["cluster"] == cid].copy()
                if not sub.empty:
                    parts.append(sub)
                    got_any = True
            if df_fibro is not None:
                sub_f = df_fibro.loc[df_fibro["cluster"] == cid].copy()
                if not sub_f.empty:
                    parts.append(sub_f)
                    got_any = True
                    if _has_fibro_data(sub_f):
                        fibro_success = True
            if got_any:
                found_cluster = cid
                break

        if not parts:
            print(f"Skipping {requested}: cluster not found in value tables")
            continue

        # Merge value tables by index/columns so fibroblast sample columns are retained.
        merged = parts[0].copy()
        for nxt in parts[1:]:
            merged = merged.combine_first(nxt)

        numeric_cols = [c for c in merged.columns if c != "cluster"]
        sample_cols = [c for c in numeric_cols if pd.api.types.is_numeric_dtype(merged[c])]

        # Keep only requested cell groups. Include fibroblast samples only if cluster
        # is present in the fibroblast success table.
        allowed_groups = [g for g in CELL_GROUP_ORDER if fibro_success or g != "Fibroblast"]
        sample_cols = [c for c in sample_cols if _sample_cell_group(c) in allowed_groups]
        sample_cols = _sort_samples(sample_cols)

        if not sample_cols:
            print(f"Skipping {requested}: no matching sample columns for requested cell groups")
            continue

        gene = _first_gene_label(gene_map.get(found_cluster, ""))
        _plot_cluster(merged, sample_cols, gene, found_cluster, OUTPUT_DIR)
        print(f"Plotted {found_cluster} -> {gene if gene else '(no gene)'}")


if __name__ == "__main__":
    run()
