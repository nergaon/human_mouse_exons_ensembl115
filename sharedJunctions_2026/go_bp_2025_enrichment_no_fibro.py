#!/usr/bin/env python3
"""
GO Biological Process 2025 enrichment for no-fibroblast cell types.

Builds 4 gene lists from cluster-level status across 3 comparisons:
  1) sig_in_any_hm_cell
  2) unchanged_in_any_hm_cell
  3) sig_hm_and_unchanged_hh_mm
  4) unchanged_in_all_three

Background genes:
  genes from clusters that are successful in at least one HM cell type
  (success = cluster has at least one numeric p.adjust in HM for that cell type).

Fibroblast data are intentionally excluded.
"""

from __future__ import annotations

from pathlib import Path
import json
import multiprocessing as mp
import re

import pandas as pd


BASE = Path("/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026")
FILE_HM = BASE / "leafcutter_GSE115736_GSE116177" / "clusters_sum_table_HN6.txt"
FILE_HH = BASE / "leafcutter_GSE115736_GSE60424" / "clusters_sum_table_HN6.txt"
FILE_MM = BASE / "leafcutter_GSE116177_GSE180020" / "clusters_sum_table_HN6.txt"
OUT_DIR = BASE / "go_bp_2025_no_fibro_results"

THRESHOLD = 0.05
SIG_DELTAPSI_THRESHOLD = 0.1
UNCHANGED_DELTAPSI_THRESHOLD = 0.05
ENRICHR_TIMEOUT_SEC = 180

# canonical CT name -> (prefix_in_HM, prefix_in_HH, prefix_in_MM)
CELL_TYPES: dict[str, tuple[str, str, str]] = {
    "CD4T": ("CD4T", "CD4T", "CD4T"),
    "CD8T": ("CD8T", "CD8T", "Cd8T"),
    "NveB": ("NveB", "NveB", "BCell"),
    "NK": ("NK", "NK", "NK"),
    "Mono": ("Mono", "Mono", "Mono"),
}


def padj_col(prefix: str) -> str:
    return f"{prefix}_p.adjust"


def abs_deltapsi_col(prefix: str) -> str:
    return f"{prefix}_abs_deltapsi"


def cluster_min_padj_lookup(df: pd.DataFrame, prefix: str) -> pd.Series:
    col = padj_col(prefix)
    if col not in df.columns:
        return pd.Series(dtype=float)
    return (
        pd.to_numeric(df[col], errors="coerce")
        .groupby(df["cluster"], dropna=False)
        .min()
        .astype(float)
    )


def cluster_max_abs_deltapsi_lookup(df: pd.DataFrame, prefix: str) -> pd.Series:
    col = abs_deltapsi_col(prefix)
    if col not in df.columns:
        return pd.Series(dtype=float)
    vals = pd.to_numeric(df[col], errors="coerce").abs()
    return vals.groupby(df["cluster"], dropna=False).max().astype(float)


def cluster_success_clusters(df: pd.DataFrame, prefix: str) -> set[str]:
    col = padj_col(prefix)
    if col not in df.columns:
        return set()
    vals = pd.to_numeric(df[col], errors="coerce")
    has_value = vals.groupby(df["cluster"], dropna=False).apply(lambda s: s.notna().any())
    return set(has_value[has_value].index.astype(str))


def classify_cluster_status(min_padj: float | None, max_abs_dpsi: float | None) -> str:
    if min_padj is None or max_abs_dpsi is None:
        return ""
    if pd.isna(min_padj) or pd.isna(max_abs_dpsi):
        return ""
    if (min_padj <= THRESHOLD) and (max_abs_dpsi >= SIG_DELTAPSI_THRESHOLD):
        return "sig"
    if (max_abs_dpsi < UNCHANGED_DELTAPSI_THRESHOLD) or (min_padj > THRESHOLD):
        return "unchanged"
    return ""


def split_genes(genes_text: str) -> list[str]:
    text = str(genes_text).strip()
    if not text or text.lower() == "nan":
        return []
    parts = re.split(r"[;,|]", text)
    genes = []
    for p in parts:
        g = p.strip()
        if g and g.lower() != "nan":
            genes.append(g)
    return genes


def clusters_to_gene_set(clusters: set[str], cluster_to_genes: pd.Series) -> set[str]:
    out: set[str] = set()
    for cl in clusters:
        for gene in split_genes(cluster_to_genes.get(cl, "")):
            out.add(gene)
    return out


def write_gene_list(path: Path, genes: set[str]) -> None:
    path.write_text("\n".join(sorted(genes)) + ("\n" if genes else ""))


def _enrichr_worker(
    genes: list[str],
    background: list[str],
    run_dir: str,
    out_tsv: str,
    queue: mp.Queue,
) -> None:
    try:
        import gseapy as gp

        enr = gp.enrichr(
            gene_list=genes,
            gene_sets="GO_Biological_Process_2025",
            organism="Human",
            background=background,
            outdir=run_dir,
            no_plot=True,
            cutoff=1.0,
        )

        results = getattr(enr, "results", pd.DataFrame())
        if isinstance(results, pd.DataFrame) and not results.empty:
            results.to_csv(out_tsv, sep="\t", index=False)
            queue.put({"status": "ok", "rows": int(len(results))})
        else:
            queue.put({"status": "empty"})
    except Exception as exc:
        queue.put({"status": "error", "message": str(exc)})


def run_enrichr_with_background(label: str, genes: set[str], background: set[str], out_dir: Path) -> None:
    run_dir = out_dir / label
    run_dir.mkdir(parents=True, exist_ok=True)

    if not genes:
        (run_dir / "EMPTY_GENE_LIST.txt").write_text("No genes in this list.\n")
        return

    try:
        import gseapy  # noqa: F401  # pre-check dependency before launching subprocess
    except Exception as exc:
        (run_dir / "ERROR_missing_gseapy.txt").write_text(
            "gseapy is required for enrichment. Install with: pip install gseapy\n"
            f"Import error: {exc}\n"
        )
        return

    queue: mp.Queue = mp.Queue()
    out_tsv = out_dir / f"{label}_GO_Biological_Process_2025.tsv"
    proc = mp.Process(
        target=_enrichr_worker,
        args=(sorted(genes), sorted(background), str(run_dir), str(out_tsv), queue),
        daemon=True,
    )
    proc.start()
    proc.join(ENRICHR_TIMEOUT_SEC)

    if proc.is_alive():
        proc.terminate()
        proc.join()
        (run_dir / "ERROR_enrichr_timeout.txt").write_text(
            "Enrichr request timed out and was terminated.\n"
            f"Timeout: {ENRICHR_TIMEOUT_SEC} seconds\n"
        )
        return

    if queue.empty():
        (run_dir / "ERROR_enrichr_request_failed.txt").write_text(
            "Enrichr process exited without returning a result.\n"
        )
        return

    result = queue.get()
    status = result.get("status", "error")
    if status == "ok":
        return
    if status == "empty":
        (run_dir / "NO_ENRICHMENT_RESULTS.txt").write_text("No terms returned by Enrichr.\n")
    else:
        (run_dir / "ERROR_enrichr_request_failed.txt").write_text(
            "Enrichr request failed (likely no network access to maayanlab.cloud).\n"
            f"Error: {result.get('message', 'Unknown error')}\n"
        )


def main() -> None:
    for p in (FILE_HM, FILE_HH, FILE_MM):
        if not p.exists():
            raise FileNotFoundError(f"Missing input file: {p}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df_hm = pd.read_csv(FILE_HM, sep="\t")
    df_hh = pd.read_csv(FILE_HH, sep="\t")
    df_mm = pd.read_csv(FILE_MM, sep="\t")

    cluster_to_genes = (
        df_hm.drop_duplicates(subset="cluster")[["cluster", "genes"]]
        .assign(cluster=lambda x: x["cluster"].astype(str))
        .set_index("cluster")["genes"]
    )

    sig_any_hm_clusters: set[str] = set()
    unchanged_any_hm_clusters: set[str] = set()
    sig_hm_unch_hh_mm_clusters: set[str] = set()
    unchanged_all_three_clusters: set[str] = set()
    background_success_clusters: set[str] = set()

    for ct, (pfx_hm, pfx_hh, pfx_mm) in CELL_TYPES.items():
        hm_min = cluster_min_padj_lookup(df_hm, pfx_hm)
        hm_max = cluster_max_abs_deltapsi_lookup(df_hm, pfx_hm)
        hh_min = cluster_min_padj_lookup(df_hh, pfx_hh)
        hh_max = cluster_max_abs_deltapsi_lookup(df_hh, pfx_hh)
        mm_min = cluster_min_padj_lookup(df_mm, pfx_mm)
        mm_max = cluster_max_abs_deltapsi_lookup(df_mm, pfx_mm)

        background_success_clusters |= cluster_success_clusters(df_hm, pfx_hm)

        clusters_with_all = (
            set(hm_min.index.astype(str))
            & set(hm_max.index.astype(str))
            & set(hh_min.index.astype(str))
            & set(hh_max.index.astype(str))
            & set(mm_min.index.astype(str))
            & set(mm_max.index.astype(str))
        )

        for cl in clusters_with_all:
            s_hm = classify_cluster_status(hm_min.get(cl), hm_max.get(cl))
            s_hh = classify_cluster_status(hh_min.get(cl), hh_max.get(cl))
            s_mm = classify_cluster_status(mm_min.get(cl), mm_max.get(cl))

            if s_hm == "sig":
                sig_any_hm_clusters.add(cl)
            if s_hm == "unchanged":
                unchanged_any_hm_clusters.add(cl)
            if (s_hm == "sig") and (s_hh == "unchanged") and (s_mm == "unchanged"):
                sig_hm_unch_hh_mm_clusters.add(cl)
            if (s_hm == "unchanged") and (s_hh == "unchanged") and (s_mm == "unchanged"):
                unchanged_all_three_clusters.add(cl)

    background_genes = clusters_to_gene_set(background_success_clusters, cluster_to_genes)

    gene_sets: dict[str, set[str]] = {
        "sig_in_any_hm_cell": clusters_to_gene_set(sig_any_hm_clusters, cluster_to_genes),
        "unchanged_in_any_hm_cell": clusters_to_gene_set(unchanged_any_hm_clusters, cluster_to_genes),
        "sig_hm_and_unchanged_hh_mm": clusters_to_gene_set(sig_hm_unch_hh_mm_clusters, cluster_to_genes),
        "unchanged_in_all_three": clusters_to_gene_set(unchanged_all_three_clusters, cluster_to_genes),
    }

    write_gene_list(OUT_DIR / "background_success_genes.txt", background_genes)
    for label, genes in gene_sets.items():
        write_gene_list(OUT_DIR / f"{label}.txt", genes)

    summary = {
        "background_success_clusters": len(background_success_clusters),
        "background_success_genes": len(background_genes),
    }
    for label, genes in gene_sets.items():
        summary[f"{label}_genes"] = len(genes)

    (OUT_DIR / "summary_counts.json").write_text(json.dumps(summary, indent=2) + "\n")

    total_sets = len(gene_sets)
    for idx, (label, genes) in enumerate(gene_sets.items(), start=1):
        print(f"[{idx}/{total_sets}] Running Enrichr for: {label}", flush=True)
        run_enrichr_with_background(label, genes, background_genes, OUT_DIR)
        print(f"[{idx}/{total_sets}] Finished: {label}", flush=True)

    print("Done. Results written to:")
    print(f"  {OUT_DIR}")


if __name__ == "__main__":
    main()
