#!/usr/bin/env python3
"""
GO Biological Process 2025 enrichment across all requested cell types.

Uses a single workbook:
    unique_sig_clusters_HN6.xlsx (sheet: all_cluster_status)

For each cell type, gene lists are built from status labels:
    - splicing unchanged
    - differentially spliced

Background genes are taken from all rows of column B (genes) in
all_cluster_status.
"""

from __future__ import annotations

from pathlib import Path
import json
import multiprocessing as mp
import re

import pandas as pd


BASE = Path("/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026")
FILE_XLSX = BASE / "leafcutter_GSE115736_GSE116177" / "unique_sig_clusters_HN6.xlsx"
SHEET_ALL_CLUSTER_STATUS = "all_cluster_status"
OUT_DIR = BASE / "go_bp_2025_with_fibro_neu_results"
ENRICHMENT_FILE_SUFFIX = "_GO_BP_enrichment.tsv"

THRESHOLD = 0.05
ENRICHR_TIMEOUT_SEC = 180
GO_LIBRARIES: tuple[str, ...] = ("GO_Biological_Process_2025", "GO_Biological_Process_2023")

REQUESTED_CELL_TYPES: tuple[str, ...] = (
    "CD4T",
    "CD8T",
    "NveB",
    "NK",
    "Mono",
    "Fibroblast",
    "Neu",
)

CELL_TYPE_CANDIDATES: dict[str, tuple[str, ...]] = {
    "CD4T": ("CD4T",),
    "CD8T": ("CD8T",),
    "NveB": ("NveB", "BCell"),
    "NK": ("NK",),
    "Mono": ("Mono",),
    "Fibroblast": ("Fibroblast", "Fibroblasts", "fibroblast", "fibroblasts"),
    "Neu": ("Neu", "Neut", "Neutrophils"),
}

CELL_TYPE_DISPLAY: dict[str, str] = {
    "CD4T": "CD4T",
    "CD8T": "CD8T",
    "NveB": "NaveB",
    "NK": "NK",
    "Mono": "Mono",
    "Fibroblast": "Fibroblast",
    "Neu": "Neu",
}


def resolve_cell_type_columns(df: pd.DataFrame) -> list[tuple[str, str]]:
    """Map requested cell type labels to actual column names present in the sheet."""
    resolved: list[tuple[str, str]] = []
    missing: list[str] = []

    for ct in REQUESTED_CELL_TYPES:
        candidates = CELL_TYPE_CANDIDATES.get(ct, (ct,))
        actual = next((col for col in candidates if col in df.columns), None)
        if actual is None:
            missing.append(f"{ct} (tried: {', '.join(candidates)})")
            continue
        resolved.append((ct, actual))

    if missing:
        raise ValueError(
            f"Missing expected cell type columns in sheet '{SHEET_ALL_CLUSTER_STATUS}': "
            f"{'; '.join(missing)}"
        )

    return resolved


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


def genes_series_to_set(genes_series: pd.Series) -> set[str]:
    out: set[str] = set()
    for genes_text in genes_series:
        for gene in split_genes(genes_text):
            out.add(gene)
    return out


def _status_mask(status_series: pd.Series, keyword: str) -> pd.Series:
    s = status_series.fillna("").astype(str).str.strip().str.lower()
    return s.str.contains(keyword, regex=False)


def write_gene_list(path: Path, genes: set[str]) -> None:
    path.write_text("\n".join(sorted(genes)) + ("\n" if genes else ""))


def _go_worker(
    genes: list[str],
    background: list[str],
    library: str,
    run_dir: str,
    out_tsv: str,
    queue: mp.Queue,
) -> None:
    try:
        import gseapy as gp

        enr = gp.enrichr(
            gene_list=genes,
            gene_sets=library,
            organism="human",
            outdir=None,
            background=background,
            cutoff=1.0,
            no_plot=True,
        )

        results = getattr(enr, "results", None)
        if not isinstance(results, pd.DataFrame):
            results = pd.DataFrame() if results is None else pd.DataFrame(results)

        if isinstance(results, pd.DataFrame) and not results.empty:
            req_cols = {"Overlap", "Odds Ratio"}
            if req_cols.issubset(set(results.columns)):
                results["Enrichment Ratio"] = pd.to_numeric(results["Odds Ratio"], errors="coerce")

            results.to_csv(out_tsv, sep="\t", index=False)
            queue.put({"status": "ok", "rows": int(len(results))})
        else:
            queue.put({"status": "empty"})
    except Exception as exc:
        queue.put({"status": "error", "message": str(exc)})


def run_go_with_background(label: str, genes: set[str], background: set[str], library: str, out_dir: Path) -> None:
    run_dir = out_dir / label
    run_dir.mkdir(parents=True, exist_ok=True)
    for stale_name in (
        "ERROR_missing_gseapy.txt",
        "ERROR_go_request_failed.txt",
        "ERROR_go_timeout.txt",
        "NO_ENRICHMENT_RESULTS.txt",
        "EMPTY_GENE_LIST.txt",
    ):
        stale_path = run_dir / stale_name
        if stale_path.exists():
            stale_path.unlink()

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
    out_tsv = out_dir / f"{label}{ENRICHMENT_FILE_SUFFIX}"
    proc = mp.Process(
        target=_go_worker,
        args=(sorted(genes), sorted(background), library, str(run_dir), str(out_tsv), queue),
        daemon=True,
    )
    proc.start()
    proc.join(ENRICHR_TIMEOUT_SEC)

    if proc.is_alive():
        proc.terminate()
        proc.join()
        (run_dir / "ERROR_go_timeout.txt").write_text(
            "GO request timed out and was terminated.\n"
            f"Timeout: {ENRICHR_TIMEOUT_SEC} seconds\n"
        )
        return

    if queue.empty():
        (run_dir / "ERROR_go_request_failed.txt").write_text(
            "GO process exited without returning a result.\n"
        )
        return

    result = queue.get()
    status = result.get("status", "error")
    if status == "ok":
        return
    if status == "empty":
        (run_dir / "NO_ENRICHMENT_RESULTS.txt").write_text("No terms returned by Enrichr.\n")
    else:
        (run_dir / "ERROR_go_request_failed.txt").write_text(
            "GO request failed (likely no network access to Enrichr).\n"
            f"Error: {result.get('message', 'Unknown error')}\n"
        )


def _first_existing_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def build_go_table(rows_meta: list[dict[str, str]], out_dir: Path, significant_only: bool) -> pd.DataFrame:
    out_rows: list[dict[str, object]] = []
    for meta in rows_meta:
        tsv_path = out_dir / f"{meta['label']}{ENRICHMENT_FILE_SUFFIX}"
        if not tsv_path.exists():
            continue

        try:
            df = pd.read_csv(tsv_path, sep="\t")
        except Exception:
            continue
        if df.empty:
            continue

        term_col = _first_existing_col(df, ["Term", "term", "name"])
        pval_col = _first_existing_col(df, ["P-value", "P value", "p_value", "pval", "PValue"])
        fdr_col = _first_existing_col(df, ["Adjusted P-value", "Adjusted P value", "adj_p", "FDR", "fdr", "p_value"])
        enr_col = _first_existing_col(df, ["Odds Ratio", "Odds_Ratio", "odds_ratio", "Enrichment Ratio", "enrichment_ratio"])

        if term_col is None:
            continue

        use_cols = [term_col]
        if pval_col:
            use_cols.append(pval_col)
        if fdr_col:
            use_cols.append(fdr_col)
        if enr_col:
            use_cols.append(enr_col)

        tmp = df[use_cols].copy()
        if pval_col:
            tmp[pval_col] = pd.to_numeric(tmp[pval_col], errors="coerce")
        else:
            tmp["_pvalue"] = pd.NA
            pval_col = "_pvalue"

        if fdr_col:
            tmp[fdr_col] = pd.to_numeric(tmp[fdr_col], errors="coerce")
        else:
            tmp["_fdr"] = pd.NA
            fdr_col = "_fdr"

        if significant_only:
            tmp = tmp[tmp[fdr_col] <= THRESHOLD]
        if tmp.empty:
            continue

        if enr_col is None:
            tmp["_enrichment_ratio"] = pd.NA
            enr_use = "_enrichment_ratio"
        else:
            tmp[enr_col] = pd.to_numeric(tmp[enr_col], errors="coerce")
            enr_use = enr_col

        for _, row in tmp.iterrows():
            out_rows.append(
                {
                    "library": meta["library"],
                    "cell_type": meta["cell_type"],
                    "type_of_list": meta["type_of_list"],
                    "go_term": row[term_col],
                    "p_value": row[pval_col],
                    "enrichment_ratio": row[enr_use],
                    "fdr": row[fdr_col],
                }
            )

    if not out_rows:
        return pd.DataFrame(columns=["library", "cell_type", "type_of_list", "go_term", "p_value", "enrichment_ratio", "fdr"])

    out_df = pd.DataFrame(out_rows)
    out_df = out_df.sort_values(["library", "cell_type", "type_of_list", "fdr", "p_value", "go_term"], kind="stable").reset_index(drop=True)
    return out_df


def _collect_run_errors(rows_meta: list[dict[str, str]], out_dir: Path) -> list[str]:
    errors: list[str] = []
    for meta in rows_meta:
        run_dir = out_dir / meta["label"]
        for err_name in ("ERROR_missing_gseapy.txt", "ERROR_go_request_failed.txt", "ERROR_go_timeout.txt"):
            err_path = run_dir / err_name
            if err_path.exists():
                try:
                    msg = err_path.read_text().strip().replace("\n", " | ")
                except Exception:
                    msg = "(could not read error file)"
                errors.append(f"{meta['label']}: {msg}")
                break
    return errors


def main() -> None:
    if not FILE_XLSX.exists():
        raise FileNotFoundError(f"Missing input file: {FILE_XLSX}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df_status = pd.read_excel(FILE_XLSX, sheet_name=SHEET_ALL_CLUSTER_STATUS)

    if df_status.shape[1] < 2:
        raise ValueError(
            f"Sheet '{SHEET_ALL_CLUSTER_STATUS}' must include at least 2 columns; "
            "column B should contain genes."
        )

    genes_col = df_status.columns[1]
    cell_type_columns = resolve_cell_type_columns(df_status)

    background_genes = genes_series_to_set(df_status[genes_col])
    all_diff_genes: set[str] = set()
    all_unchanged_genes: set[str] = set()

    rows_meta: list[dict[str, str]] = []
    summary: dict[str, object] = {}

    for ct, status_col in cell_type_columns:
        unchanged_mask = _status_mask(df_status[status_col], "splicing unchanged")
        diff_mask = _status_mask(df_status[status_col], "differentially spliced")

        unchanged_genes = genes_series_to_set(df_status.loc[unchanged_mask, genes_col])
        diff_genes = genes_series_to_set(df_status.loc[diff_mask, genes_col])
        all_diff_genes.update(diff_genes)
        all_unchanged_genes.update(unchanged_genes)

        ct_dir = OUT_DIR / ct
        ct_dir.mkdir(parents=True, exist_ok=True)
        write_gene_list(ct_dir / "background_genes_all_cluster_status_colB.txt", background_genes)
        write_gene_list(ct_dir / "differentially_spliced.txt", diff_genes)
        write_gene_list(ct_dir / "splicing_unchanged.txt", unchanged_genes)

        for library in GO_LIBRARIES:
            label_diff = f"{ct}_differentially_spliced_{library}"
            label_unchanged = f"{ct}_splicing_unchanged_{library}"

            print(f"[{ct}] Running {library} enrichment for differentially spliced", flush=True)
            run_go_with_background(label_diff, diff_genes, background_genes, library, OUT_DIR)
            print(f"[{ct}] Running {library} enrichment for splicing unchanged", flush=True)
            run_go_with_background(label_unchanged, unchanged_genes, background_genes, library, OUT_DIR)

            rows_meta.append(
                {
                    "label": label_diff,
                    "library": library,
                    "cell_type": CELL_TYPE_DISPLAY.get(ct, ct),
                    "type_of_list": "differentially spliced",
                }
            )
            rows_meta.append(
                {
                    "label": label_unchanged,
                    "library": library,
                    "cell_type": CELL_TYPE_DISPLAY.get(ct, ct),
                    "type_of_list": "splicing unchanged",
                }
            )

        summary[ct] = {
            "background_genes_all_cluster_status_colB": len(background_genes),
            "differentially_spliced_clusters": int(diff_mask.sum()),
            "differentially_spliced_genes": len(diff_genes),
            "splicing_unchanged_clusters": int(unchanged_mask.sum()),
            "splicing_unchanged_genes": len(unchanged_genes),
        }

    (OUT_DIR / "summary_counts.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_gene_list(OUT_DIR / "all_differentially_spliced_genes.txt", all_diff_genes)
    write_gene_list(OUT_DIR / "all_splicing_unchanged_genes.txt", all_unchanged_genes)

    final_all_table = build_go_table(rows_meta, OUT_DIR, significant_only=False)
    if final_all_table.empty:
        run_errors = _collect_run_errors(rows_meta, OUT_DIR)
        preview = "\n".join(run_errors[:8]) if run_errors else "No per-run error files found."
        raise RuntimeError(
            "No enrichment terms were collected for any cell type/library. "
            "Check dependency/network errors in run directories.\n"
            f"Error preview:\n{preview}"
        )

    final_all_tsv = OUT_DIR / "go_bp_2023_2025_all_terms_table.tsv"
    final_all_csv = OUT_DIR / "go_bp_2023_2025_all_terms_table.csv"
    final_all_table.to_csv(final_all_tsv, sep="\t", index=False)
    final_all_table.to_csv(final_all_csv, index=False)

    final_sig_table = build_go_table(rows_meta, OUT_DIR, significant_only=True)
    final_sig_tsv = OUT_DIR / "go_bp_2023_2025_significant_terms_table.tsv"
    final_sig_csv = OUT_DIR / "go_bp_2023_2025_significant_terms_table.csv"
    final_sig_table.to_csv(final_sig_tsv, sep="\t", index=False)
    final_sig_table.to_csv(final_sig_csv, index=False)

    if not final_sig_table.empty:
        sig_counts = (
            final_sig_table.groupby(["library", "cell_type", "type_of_list"], as_index=False)
            .size()
            .rename(columns={"size": "n_significant_terms_fdr_le_0_05"})
            .sort_values(["library", "cell_type", "type_of_list"], kind="stable")
        )
    else:
        sig_counts = pd.DataFrame(
            columns=["library", "cell_type", "type_of_list", "n_significant_terms_fdr_le_0_05"]
        )
    sig_counts.to_csv(OUT_DIR / "significant_counts_by_library.tsv", sep="\t", index=False)
    sig_counts.to_csv(OUT_DIR / "significant_counts_by_library.csv", index=False)

    all_terms_txt = OUT_DIR / "go_bp_2023_2025_all_terms_table.txt"
    significant_terms_txt = OUT_DIR / "go_bp_2023_2025_significant_terms_table.txt"
    significant_counts_txt = OUT_DIR / "significant_counts_by_library.txt"

    print("Done. Results written to:")
    print(f"  {OUT_DIR}")
    print("Final all-terms GO table:")
    print(f"  {final_all_tsv}")
    print("Final significant GO table:")
    print(f"  {final_sig_tsv}")
    all_terms_txt.write_text(final_all_table.to_string(index=False) + "\n")
    significant_terms_txt.write_text(final_sig_table.to_string(index=False) + "\n")
    significant_counts_txt.write_text(sig_counts.to_string(index=False) + "\n")
    print("Saved full table-style outputs:")
    print(f"  {all_terms_txt}")
    print(f"  {significant_terms_txt}")
    print(f"  {significant_counts_txt}")


if __name__ == "__main__":
    main()
