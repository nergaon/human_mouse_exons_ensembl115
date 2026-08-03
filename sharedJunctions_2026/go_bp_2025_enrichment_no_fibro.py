#!/usr/bin/env python3
"""
GO Biological Process 2025 enrichment for one global SU list and one global DS list.

Uses a single workbook:
    unique_sig_clusters_HN6.xlsx (sheet: all_cluster_status)

Gene lists are built across the requested immune cell types:
    - SU: genes seen only in splicing unchanged rows
    - DS: genes seen only in differentially spliced rows

Background genes are taken from all rows of column B (genes) in
all_cluster_status.
"""

from __future__ import annotations

from pathlib import Path
import multiprocessing as mp
import re
import math

import pandas as pd
from gseapy import Msigdb


BASE = Path("/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026")
FILE_XLSX = BASE / "leafcutter_GSE115736_GSE116177" / "unique_sig_clusters_HN6.xlsx"
SHEET_ALL_CLUSTER_STATUS = "all_cluster_status"
OUT_DIR = BASE / "go_bp_2025"
ENRICHMENT_FILE_SUFFIX = "_GO_BP_enrichment.tsv"

THRESHOLD = 0.05
ENRICHR_TIMEOUT_SEC = 180
GO_LIBRARIES: tuple[str, ...] = ("GO_Biological_Process_2025", "hallmark", "Reactome_2022")

MSIGDB_LIBRARY_CATEGORY: dict[str, str] = {
    "GO_Biological_Process_2025": "c5.go.bp",
    "hallmark": "h.all",
    "Reactome_2022": "c2.cp.reactome",
}
MSIGDB_DB_VERSION = "2024.1.Hs"

REQUESTED_CELL_TYPES: tuple[str, ...] = (
    "CD4T",
    "CD8T",
    "NveB",
    "NK",
    "Mono",
)

PRIMARY_IMMUNE_CELL_TYPES: tuple[str, ...] = ("CD4T", "CD8T", "NveB", "NK", "Mono")

CELL_TYPE_CANDIDATES: dict[str, tuple[str, ...]] = {
    "CD4T": ("CD4T",),
    "CD8T": ("CD8T",),
    "NveB": ("NveB", "BCell"),
    "NK": ("NK",),
    "Mono": ("Mono",),
}

CELL_TYPE_DISPLAY: dict[str, str] = {
    "CD4T": "CD4T",
    "CD8T": "CD8T",
    "NveB": "B",
    "NK": "NK",
    "Mono": "Mono",
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
    parts = re.split(r"[;|]", text)
    genes = []
    for p in parts:
        # Some entries are "GENE,LOC..." aliases; keep only the first gene symbol.
        g = p.split(",", 1)[0].strip()
        if g and g.lower() != "nan":
            genes.append(g)
    return genes


def genes_series_to_set(genes_series: pd.Series) -> set[str]:
    out: set[str] = set()
    for genes_text in genes_series:
        for gene in split_genes(genes_text):
            out.add(gene)
    return out


def split_result_genes(genes_text: object) -> list[str]:
    text = str(genes_text).strip()
    if not text or text.lower() == "nan":
        return []
    parts = re.split(r"[;,|]", text)
    genes: list[str] = []
    for p in parts:
        g = p.strip()
        if g and g.lower() != "nan":
            genes.append(g)
    return genes


def join_genes(genes: set[str]) -> str:
    return ";".join(sorted(genes))


def _status_mask(status_series: pd.Series, keyword: str) -> pd.Series:
    s = status_series.fillna("").astype(str).str.strip().str.lower()
    return s.str.contains(keyword, regex=False)


def _go_worker(
    genes: list[str],
    background: list[str],
    library: str,
    out_tsv: str,
    queue: mp.Queue,
) -> None:
    try:
        import gseapy as gp

        category = MSIGDB_LIBRARY_CATEGORY.get(library)
        if category is None:
            queue.put({"status": "error", "message": f"Unsupported library: {library}"})
            return

        try:
            msig = Msigdb()
            gmt = msig.get_gmt(category=category, dbver=MSIGDB_DB_VERSION)
            results = gp.enrich(
                gene_list=genes,
                gene_sets=gmt,
                background=background,
                outdir=None,
            )
        except Exception as exc:
            queue.put({"status": "error", "message": f"{library} loading failed: {str(exc)}"})
            return

        result_df = getattr(results, "results", None)
        if not isinstance(result_df, pd.DataFrame):
            result_df = pd.DataFrame() if result_df is None else pd.DataFrame(result_df)

        if isinstance(result_df, pd.DataFrame) and not result_df.empty:
            req_cols = {"Overlap", "Odds Ratio"}
            if req_cols.issubset(set(result_df.columns)):
                result_df["Enrichment Ratio"] = pd.to_numeric(result_df["Odds Ratio"], errors="coerce")

            # Remove terms with infinite ratio values; these come from zero denominators
            # in odds-ratio calculations and are not informative for ranking.
            drop_mask = pd.Series(False, index=result_df.index)
            for ratio_col in ("Odds Ratio", "Enrichment Ratio", "odds_ratio", "enrichment_ratio"):
                if ratio_col in result_df.columns:
                    vals = pd.to_numeric(result_df[ratio_col], errors="coerce")
                    drop_mask = drop_mask | (~vals.isna() & ~vals.map(math.isfinite))
            if drop_mask.any():
                result_df = result_df.loc[~drop_mask].copy()

            if result_df.empty:
                queue.put({"status": "empty"})
            else:
                result_df.to_csv(out_tsv, sep="\t", index=False)
                queue.put({"status": "ok", "rows": int(len(result_df))})
        else:
            queue.put({"status": "empty"})
    except Exception as exc:
        queue.put({"status": "error", "message": str(exc)})


def run_go_with_background(label: str, genes: set[str], background: set[str], library: str, out_dir: Path) -> None:
    # Keep all status/error artifacts as flat files in OUT_DIR to avoid creating
    # per-run folders that are often empty.
    stale_paths = [
        out_dir / f"{label}.ERROR_missing_gseapy.txt",
        out_dir / f"{label}.ERROR_go_request_failed.txt",
        out_dir / f"{label}.ERROR_go_timeout.txt",
        out_dir / f"{label}.NO_ENRICHMENT_RESULTS.txt",
        out_dir / f"{label}.EMPTY_GENE_LIST.txt",
    ]
    for stale_path in stale_paths:
        if stale_path.exists():
            stale_path.unlink()

    if not genes:
        (out_dir / f"{label}.EMPTY_GENE_LIST.txt").write_text("No genes in this list.\n")
        return

    try:
        import gseapy  # noqa: F401  # pre-check dependency before launching subprocess
    except Exception as exc:
        (out_dir / f"{label}.ERROR_missing_gseapy.txt").write_text(
            "gseapy is required for enrichment. Install with: pip install gseapy\n"
            f"Import error: {exc}\n"
        )
        return

    queue: mp.Queue = mp.Queue()
    out_tsv = out_dir / f"{label}{ENRICHMENT_FILE_SUFFIX}"
    # Avoid carrying over old enrichment tables when a new run fails or times out.
    if out_tsv.exists():
        out_tsv.unlink()
    proc = mp.Process(
        target=_go_worker,
        args=(sorted(genes), sorted(background), library, str(out_tsv), queue),
        daemon=True,
    )
    proc.start()
    proc.join(ENRICHR_TIMEOUT_SEC)

    if proc.is_alive():
        proc.terminate()
        proc.join()
        (out_dir / f"{label}.ERROR_go_timeout.txt").write_text(
            "GO request timed out and was terminated.\n"
            f"Timeout: {ENRICHR_TIMEOUT_SEC} seconds\n"
        )
        return

    if queue.empty():
        (out_dir / f"{label}.ERROR_go_request_failed.txt").write_text(
            "GO process exited without returning a result.\n"
        )
        return

    result = queue.get()
    status = result.get("status", "error")
    if status == "ok":
        return
    if status == "empty":
        (out_dir / f"{label}.NO_ENRICHMENT_RESULTS.txt").write_text("No enrichment terms were returned.\n")
    else:
        (out_dir / f"{label}.ERROR_go_request_failed.txt").write_text(
            "GO request failed while loading or evaluating the requested gene set library.\n"
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
        term_genes_col = _first_existing_col(df, ["Genes", "genes", "Lead_genes", "lead_genes"])

        if term_col is None:
            continue

        use_cols = [term_col]
        if pval_col:
            use_cols.append(pval_col)
        if fdr_col:
            use_cols.append(fdr_col)
        if enr_col:
            use_cols.append(enr_col)
        if term_genes_col and term_genes_col not in use_cols:
            use_cols.append(term_genes_col)

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

        if term_genes_col is None:
            tmp["_term_genes"] = ""
            term_genes_use = "_term_genes"
        else:
            tmp[term_genes_col] = tmp[term_genes_col].fillna("").astype(str)
            term_genes_use = term_genes_col

        for _, row in tmp.iterrows():
            enr_val = row[enr_use]
            if pd.notna(enr_val) and isinstance(enr_val, (int, float)) and not math.isfinite(float(enr_val)):
                continue
            term_genes = split_result_genes(row[term_genes_use])
            out_rows.append(
                {
                    "library": meta["library"],
                    "cell_type": meta["cell_type"],
                    "type_of_list": meta["type_of_list"],
                    "background_gene_count": meta["background_gene_count"],
                    "list_gene_count": meta["list_gene_count"],
                    "list_genes": meta["list_genes"],
                    "go_term": row[term_col],
                    "p_value": row[pval_col],
                    "enrichment_ratio": enr_val,
                    "source_fdr": row[fdr_col],
                    "term_genes": ";".join(term_genes),
                    "term_genes_count": len(set(term_genes)),
                }
            )

    if not out_rows:
        return pd.DataFrame(
            columns=[
                "library",
                "cell_type",
                "type_of_list",
                "background_gene_count",
                "list_gene_count",
                "list_genes",
                "go_term",
                "p_value",
                "enrichment_ratio",
                "source_fdr",
                "term_genes",
                "term_genes_count",
            ]
        )

    out_df = pd.DataFrame(out_rows)
    out_df = out_df.sort_values(
        ["library", "cell_type", "type_of_list", "source_fdr", "p_value", "go_term"],
        kind="stable",
    ).reset_index(drop=True)
    return out_df


def _bh_fdr(pvals: pd.Series) -> pd.Series:
    """Benjamini-Hochberg FDR correction for a single family of tests."""
    numeric = pd.to_numeric(pvals, errors="coerce")
    valid = [(idx, float(p)) for idx, p in numeric.items() if pd.notna(p)]
    out = pd.Series(pd.NA, index=pvals.index, dtype="Float64")
    if not valid:
        return out

    valid.sort(key=lambda x: x[1])
    m = len(valid)
    raw_adj: list[float] = []
    for i, (_, p) in enumerate(valid, start=1):
        raw_adj.append(min(1.0, (p * m) / i))

    monotone_adj = [1.0] * m
    min_so_far = 1.0
    for i in range(m - 1, -1, -1):
        min_so_far = min(min_so_far, raw_adj[i])
        monotone_adj[i] = min_so_far

    for (idx, _), adj in zip(valid, monotone_adj):
        out.loc[idx] = adj
    return out


def apply_library_specific_fdr(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["fdr"] = pd.Series(pd.NA, index=out.index, dtype="Float64")
    for library, grp in out.groupby("library", sort=False):
        _ = library
        out.loc[grp.index, "fdr"] = _bh_fdr(grp["p_value"])
    return out


def _collect_run_errors(rows_meta: list[dict[str, object]], out_dir: Path) -> list[str]:
    errors: list[str] = []
    for meta in rows_meta:
        label = str(meta["label"])
        for suffix in (
            "ERROR_missing_gseapy.txt",
            "ERROR_go_request_failed.txt",
            "ERROR_go_timeout.txt",
        ):
            err_path = out_dir / f"{label}.{suffix}"
            if err_path.exists():
                try:
                    msg = err_path.read_text().strip().replace("\n", " | ")
                except Exception:
                    msg = "(could not read error file)"
                errors.append(f"{label}: {msg}")
                break
    return errors


def write_gene_lists_to_excel(
    excel_path: Path,
    background_genes: set[str],
    su_genes: set[str],
    ds_genes: set[str],
) -> None:
    """Write one workbook with the 3 requested global lists."""
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        pd.DataFrame({"gene": sorted(background_genes)}).to_excel(writer, sheet_name="back_gound", index=False)
        pd.DataFrame({"gene": sorted(su_genes)}).to_excel(writer, sheet_name="SU", index=False)
        pd.DataFrame({"gene": sorted(ds_genes)}).to_excel(writer, sheet_name="DS", index=False)


def write_library_result_workbook(excel_path: Path, table: pd.DataFrame) -> None:
    """Write one workbook for a single library with DS and SU sheets."""
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        for list_type in ("DS", "SU"):
            list_df = table[table["type_of_list"] == list_type].copy()
            list_df = list_df.sort_values(["fdr", "p_value", "go_term"], kind="stable")
            list_df.to_excel(writer, sheet_name=list_type, index=False)


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

    rows_meta: list[dict[str, object]] = []
    summary: dict[str, object] = {}

    for ct, status_col in cell_type_columns:
        unchanged_mask = _status_mask(df_status[status_col], "splicing unchanged")
        diff_mask = _status_mask(df_status[status_col], "differentially spliced")

        unchanged_genes = genes_series_to_set(df_status.loc[unchanged_mask, genes_col])
        diff_genes = genes_series_to_set(df_status.loc[diff_mask, genes_col])
        all_diff_genes.update(diff_genes)
        all_unchanged_genes.update(unchanged_genes)

        summary[ct] = {
            "background_genes_all_cluster_status_colB": len(background_genes),
            "differentially_spliced_clusters": int(diff_mask.sum()),
            "differentially_spliced_genes": len(diff_genes),
            "splicing_unchanged_clusters": int(unchanged_mask.sum()),
            "splicing_unchanged_genes": len(unchanged_genes),
        }

    ds_genes = all_diff_genes - all_unchanged_genes
    su_genes = all_unchanged_genes - all_diff_genes

    for library in GO_LIBRARIES:
        ds_label = f"all_differentially_spliced_{library}"
        su_label = f"all_splicing_unchanged_{library}"

        print(f"[global] Running {library} enrichment for DS", flush=True)
        run_go_with_background(ds_label, ds_genes, background_genes, library, OUT_DIR)
        print(f"[global] Running {library} enrichment for SU", flush=True)
        run_go_with_background(su_label, su_genes, background_genes, library, OUT_DIR)

        rows_meta.append(
            {
                "label": ds_label,
                "library": library,
                "cell_type": "all",
                "type_of_list": "DS",
                "background_gene_count": len(background_genes),
                "list_gene_count": len(ds_genes),
                "list_genes": join_genes(ds_genes),
            }
        )
        rows_meta.append(
            {
                "label": su_label,
                "library": library,
                "cell_type": "all",
                "type_of_list": "SU",
                "background_gene_count": len(background_genes),
                "list_gene_count": len(su_genes),
                "list_genes": join_genes(su_genes),
            }
        )

    final_all_table = build_go_table(rows_meta, OUT_DIR, significant_only=False)
    if final_all_table.empty:
        run_errors = _collect_run_errors(rows_meta, OUT_DIR)
        preview = "\n".join(run_errors[:8]) if run_errors else "No per-run error files found."
        raise RuntimeError(
            "No enrichment terms were collected for any cell type/library. "
            "Check dependency/network errors in run directories.\n"
            f"Error preview:\n{preview}"
        )

    final_all_table = apply_library_specific_fdr(final_all_table)
    final_all_table = final_all_table.sort_values(
        ["library", "cell_type", "type_of_list", "fdr", "p_value", "go_term"],
        kind="stable",
    ).reset_index(drop=True)

    # Exactly 3 result files: one per dataset/library, each containing FDR.
    library_output_names = {
        "GO_Biological_Process_2025": "results_go_biological_process_2025.tsv",
        "hallmark": "results_hallmark.tsv",
        "Reactome_2022": "results_reactome_2022.tsv",
    }
    for library, out_name in library_output_names.items():
        lib_df = final_all_table[final_all_table["library"] == library].copy()
        lib_df = lib_df.sort_values(["fdr", "p_value", "cell_type", "type_of_list", "go_term"], kind="stable")
        lib_df.to_csv(OUT_DIR / out_name, sep="\t", index=False)

        lib_stem = out_name.replace("results_", "").replace(".tsv", "")
        for stale_path in (
            OUT_DIR / f"results_{lib_stem}_DS_only.tsv",
            OUT_DIR / f"results_{lib_stem}_SU_only.tsv",
        ):
            if stale_path.exists():
                stale_path.unlink()
        write_library_result_workbook(OUT_DIR / f"results_{lib_stem}.xlsx", lib_df)

    # Gene lists workbook with requested sheet order.
    gene_lists_excel = OUT_DIR / "gene_lists.xlsx"
    write_gene_lists_to_excel(
        gene_lists_excel,
        background_genes,
        su_genes,
        ds_genes,
    )

    print("Done. Results written to:")
    print(f"  {OUT_DIR}")
    print("Gene list workbook:")
    print(f"  {gene_lists_excel}")
    print("Result files:")
    print(f"  {OUT_DIR / 'results_go_biological_process_2025.tsv'}")
    print(f"  {OUT_DIR / 'results_hallmark.tsv'}")
    print(f"  {OUT_DIR / 'results_reactome_2022.tsv'}")
    print("Per-library workbooks:")
    print(f"  {OUT_DIR / 'results_go_biological_process_2025.xlsx'}")
    print(f"  {OUT_DIR / 'results_hallmark.xlsx'}")
    print(f"  {OUT_DIR / 'results_reactome_2022.xlsx'}")


if __name__ == "__main__":
    main()
