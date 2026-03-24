#!/usr/bin/env python3
"""Filter and optionally remap junction count tables between mouse and human positions.

The script can either:
1. filter-only: keep junctions whose two breakpoints are present in a selected
    column from unique_points_HN6.txt, preserving the original junction_id
2. remap: keep junctions whose two breakpoints are present in a selected source
    column and replace them with mapped coordinates from a selected target column

Junction format is preserved as chr:start:end:strand.
"""

import argparse
import csv
from pathlib import Path


def load_position_map(points_path: Path, source_col: str, target_col: str) -> dict[str, str]:
    """Load point mapping from source_col to target_col from unique_points_HN6.txt."""
    mapping: dict[str, str] = {}
    with points_path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = set(reader.fieldnames or [])
        if source_col not in fieldnames:
            raise ValueError(f"Input points file is missing required column: {source_col}")
        if target_col not in fieldnames:
            raise ValueError(f"Input points file is missing required column: {target_col}")

        for row in reader:
            source_pos = (row.get(source_col) or "").strip()
            target_pos = (row.get(target_col) or "").strip()
            if source_pos and target_pos:
                mapping[source_pos] = target_pos
    return mapping


def remap_junction_id(junction_id: str, position_map: dict[str, str]) -> str | None:
    """Return remapped junction_id or None if junction cannot be mapped."""
    parts = junction_id.strip().split(":")
    if len(parts) < 4:
        return None

    chrom_m, start_m, end_m, strand = parts[0], parts[1], parts[2], parts[3]
    key_start_m = f"{chrom_m}:{start_m}"
    key_end_m = f"{chrom_m}:{end_m}"

    pos_start_h = position_map.get(key_start_m)
    pos_end_h = position_map.get(key_end_m)
    if not pos_start_h or not pos_end_h:
        return None

    h_start_parts = pos_start_h.split(":")
    h_end_parts = pos_end_h.split(":")
    if len(h_start_parts) < 2 or len(h_end_parts) < 2:
        return None

    chrom_h_start, start_h = h_start_parts[0], h_start_parts[1]
    chrom_h_end, end_h = h_end_parts[0], h_end_parts[1]

    # Require both mapped breakpoints to be on the same human chromosome.
    if chrom_h_start != chrom_h_end:
        return None

    return f"{chrom_h_start}:{start_h}:{end_h}:{strand}"


def filter_junction_id(junction_id: str, valid_points: set[str]) -> bool:
    """Return True if both breakpoints are present in valid_points."""
    parts = junction_id.strip().split(":")
    if len(parts) < 4:
        return False

    chrom, start, end = parts[0], parts[1], parts[2]
    return f"{chrom}:{start}" in valid_points and f"{chrom}:{end}" in valid_points


def filter_and_remap_junctions(
    input_tsv: Path,
    output_tsv: Path,
    position_map: dict[str, str],
) -> tuple[int, int]:
    kept = 0
    total = 0

    with input_tsv.open("r", newline="") as fin, output_tsv.open("w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")

        header = next(reader, None)
        if header is None:
            raise ValueError(f"Empty input file: {input_tsv}")
        if not header or header[0] != "junction_id":
            raise ValueError(f"First column must be junction_id in file: {input_tsv}")

        writer.writerow(header)

        for row in reader:
            total += 1
            if not row:
                continue

            remapped_id = remap_junction_id(row[0], position_map)
            if remapped_id is None:
                continue

            row[0] = remapped_id
            writer.writerow(row)
            kept += 1

    return total, kept


def filter_junctions(
    input_tsv: Path,
    output_tsv: Path,
    valid_points: set[str],
) -> tuple[int, int]:
    kept = 0
    total = 0

    with input_tsv.open("r", newline="") as fin, output_tsv.open("w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")

        header = next(reader, None)
        if header is None:
            raise ValueError(f"Empty input file: {input_tsv}")
        if not header or header[0] != "junction_id":
            raise ValueError(f"First column must be junction_id in file: {input_tsv}")

        writer.writerow(header)

        for row in reader:
            total += 1
            if not row:
                continue
            if filter_junction_id(row[0], valid_points):
                writer.writerow(row)
                kept += 1

    return total, kept


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Filter junction TSVs by a source coordinate column and replace with "
            "mapped coordinates from a target column."
        )
    )
    parser.add_argument("--points", required=True, help="Path to unique_points_HN6.txt")
    parser.add_argument(
        "--mode",
        choices=["remap", "filter-only"],
        default="remap",
        help="Whether to remap junction coordinates or only filter matching junctions",
    )
    parser.add_argument(
        "--source-column",
        choices=["position_m", "position_h"],
        default="position_m",
        help="Source coordinates expected in input junctions (default: position_m)",
    )
    parser.add_argument(
        "--target-column",
        choices=["position_h", "position_m"],
        default="position_h",
        help="Target coordinates to write to output junctions (default: position_h)",
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="One or more input junction count TSV files",
    )
    parser.add_argument(
        "--outputs",
        nargs="+",
        required=True,
        help="Output TSV files (same number and order as --inputs)",
    )
    args = parser.parse_args()

    if len(args.inputs) != len(args.outputs):
        raise ValueError("--inputs and --outputs must have the same number of files")
    if args.mode == "remap" and args.source_column == args.target_column:
        raise ValueError("--source-column and --target-column must be different")

    points_path = Path(args.points)
    if args.mode == "remap":
        position_map = load_position_map(points_path, args.source_column, args.target_column)
    else:
        valid_points = set(load_position_map(points_path, args.source_column, args.source_column).keys())

    for in_path_str, out_path_str in zip(args.inputs, args.outputs):
        in_path = Path(in_path_str)
        out_path = Path(out_path_str)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        if args.mode == "remap":
            total, kept = filter_and_remap_junctions(in_path, out_path, position_map)
        else:
            total, kept = filter_junctions(in_path, out_path, valid_points)
        print(f"{in_path} -> {out_path}: kept {kept}/{total} junctions")


if __name__ == "__main__":
    main()
