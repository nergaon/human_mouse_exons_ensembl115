#!/usr/bin/env python3
"""
Merge E-MTAB-5919 fibroblast junction counts with AS_clusters_value_HN6.txt.

Left join: all rows from AS_clusters_value_HN6.txt are kept.
Join key : chr:start:end  (first 3 colon-separated parts of the junction_id,
           stripping the trailing :clu_XXXX or :. field present in each file).
Missing  : filled with 0.

Outputs
-------
  AS_clusters_value_fibroblast_HN6.txt  – merged table
  leafcutter_EMTAB5919H_EMTAB5919M/groups_file.txt  – HS / MM groups

Notes
-----
This script does not classify clusters as sig/unchanged; it only prepares
the merged fibroblast count matrix and groups file. Sig/unchanged criteria
are applied downstream in merge_leafcutter_results.py.
"""

import csv
import sys
from pathlib import Path

BASE = Path("/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026")

AS_CLUSTERS = BASE / "AS_clusters_value_HN6.txt"
H_FILE      = BASE / "E-MTAB-5919_H_orthologs_junctions.tsv"
M_FILE      = BASE / "E-MTAB-5919_M_orthologs_junctions.tsv"
OUT_FILE    = BASE / "AS_clusters_value_fibroblast_HN6.txt"
LEAF_DIR    = BASE / "leafcutter_EMTAB5919H_EMTAB5919M"


def junction_key(junction_id: str) -> str:
    """chr:start:end  (strips the 4th colon-delimited field)."""
    return ":".join(junction_id.split(":")[:3])


def read_junction_file(
    path: Path,
    exclude_cols: set | None = None,
) -> tuple[list[str], dict[str, list[str]]]:
    """
    Read a junction TSV file.

    Parameters
    ----------
    path         : path to the TSV
    exclude_cols : column names to drop (e.g. ``{"original_junction_id"}``)

    Returns
    -------
    sample_cols  : ordered list of retained sample column names
    data         : mapping  junction_key -> [value, …]  for every data row
    """
    exclude_cols = exclude_cols or set()
    data: dict[str, list[str]] = {}
    sample_cols: list[str] = []

    with path.open("r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        keep_indices = [
            i for i, col in enumerate(header[1:], 1) if col not in exclude_cols
        ]
        sample_cols = [header[i] for i in keep_indices]

        for row in reader:
            if not row:
                continue
            key = junction_key(row[0])
            data[key] = [row[i] if i < len(row) else "0" for i in keep_indices]

    return sample_cols, data


def main() -> None:
    for f in (AS_CLUSTERS, H_FILE, M_FILE):
        if not f.exists():
            print(f"ERROR: required file not found: {f}", file=sys.stderr)
            sys.exit(1)

    # ── read junction files ────────────────────────────────────────────────
    h_cols, h_data = read_junction_file(H_FILE)
    m_cols, m_data = read_junction_file(M_FILE, exclude_cols={"original_junction_id"})

    h_zeros = ["0"] * len(h_cols)
    m_zeros = ["0"] * len(m_cols)

    # ── write merged output ────────────────────────────────────────────────
    matched_h = matched_m = total = 0

    with (
        AS_CLUSTERS.open("r", newline="") as fh_in,
        OUT_FILE.open("w", newline="") as fh_out,
    ):
        reader = csv.reader(fh_in, delimiter="\t")
        writer = csv.writer(fh_out, delimiter="\t", lineterminator="\n")

        as_header = next(reader)
        writer.writerow([as_header[0]] + h_cols + m_cols)

        for row in reader:
            if not row:
                continue
            total += 1
            key = junction_key(row[0])
            h_vals = h_data.get(key)
            m_vals = m_data.get(key)
            if h_vals is not None:
                matched_h += 1
            if m_vals is not None:
                matched_m += 1
            writer.writerow([row[0]] + (h_vals or h_zeros) + (m_vals or m_zeros))

    print(f"Merged  : {OUT_FILE}")
    print(f"Rows    : {total} total | {matched_h} matched H | {matched_m} matched M")

    # ── write leafcutter groups file ───────────────────────────────────────
    LEAF_DIR.mkdir(parents=True, exist_ok=True)
    groups_file = LEAF_DIR / "groups_file.txt"

    with groups_file.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for sample in h_cols:
            writer.writerow([sample, "HS"])
        for sample in m_cols:
            writer.writerow([sample, "MM"])

    print(f"Groups  : {groups_file}")
    print(f"  HS ({len(h_cols)}): {', '.join(h_cols)}")
    print(f"  MM ({len(m_cols)}): {', '.join(m_cols)}")


if __name__ == "__main__":
    main()
