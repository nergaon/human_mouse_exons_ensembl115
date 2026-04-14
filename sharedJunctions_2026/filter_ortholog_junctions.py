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


def _parse_count(value: str) -> float:
    """Parse a junction count value, supporting plain numbers and 'count/total'."""
    text = value.strip()
    if not text:
        return 0.0
    if "/" in text:
        text = text.split("/", 1)[0].strip()
    try:
        return float(text)
    except ValueError:
        return 0.0


def passes_expression_filter(row: list[str], min_reads: float, min_samples: int) -> bool:
    """Return True when at least min_samples sample columns have >= min_reads."""
    if min_samples <= 0 or min_reads <= 0:
        return True
    sample_values = row[1:]
    passing = sum(1 for value in sample_values if _parse_count(value) >= min_reads)
    return passing >= min_samples


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
    express_output_tsv: Path | None,
    min_reads: float,
    min_samples: int,
) -> tuple[int, int, int]:
    kept = 0
    express_kept = 0
    total = 0

    with input_tsv.open("r", newline="") as fin, output_tsv.open("w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
        express_fout = express_output_tsv.open("w", newline="") if express_output_tsv else None
        express_writer = (
            csv.writer(express_fout, delimiter="\t", lineterminator="\n") if express_fout else None
        )

        header = next(reader, None)
        if header is None:
            raise ValueError(f"Empty input file: {input_tsv}")
        if not header or header[0] != "junction_id":
            raise ValueError(f"First column must be junction_id in file: {input_tsv}")

        output_header = [*header, "original_junction_id"]
        writer.writerow(output_header)
        if express_writer:
            express_writer.writerow(output_header)

        for row in reader:
            total += 1
            if not row:
                continue

            original_junction_id = row[0]
            remapped_id = remap_junction_id(original_junction_id, position_map)
            if remapped_id is None:
                continue

            out_row = row.copy()
            out_row[0] = remapped_id
            out_row.append(original_junction_id)

            writer.writerow(out_row)
            kept += 1

            if express_writer and passes_expression_filter(row, min_reads, min_samples):
                express_writer.writerow(out_row)
                express_kept += 1

        if express_fout:
            express_fout.close()

    return total, kept, express_kept


def filter_junctions(
    input_tsv: Path,
    output_tsv: Path,
    valid_points: set[str],
    express_output_tsv: Path | None,
    min_reads: float,
    min_samples: int,
) -> tuple[int, int, int]:
    kept = 0
    express_kept = 0
    total = 0

    with input_tsv.open("r", newline="") as fin, output_tsv.open("w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
        express_fout = express_output_tsv.open("w", newline="") if express_output_tsv else None
        express_writer = (
            csv.writer(express_fout, delimiter="\t", lineterminator="\n") if express_fout else None
        )

        header = next(reader, None)
        if header is None:
            raise ValueError(f"Empty input file: {input_tsv}")
        if not header or header[0] != "junction_id":
            raise ValueError(f"First column must be junction_id in file: {input_tsv}")

        writer.writerow(header)
        if express_writer:
            express_writer.writerow(header)

        for row in reader:
            total += 1
            if not row:
                continue
            if filter_junction_id(row[0], valid_points):
                writer.writerow(row)
                kept += 1
                if express_writer and passes_expression_filter(row, min_reads, min_samples):
                    express_writer.writerow(row)
                    express_kept += 1

        if express_fout:
            express_fout.close()

    return total, kept, express_kept


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
    parser.add_argument(
        "--express-outputs",
        nargs="+",
        help="Optional expression-filtered outputs (same number and order as --inputs)",
    )
    parser.add_argument(
        "--min-reads",
        type=float,
        default=10.0,
        help="Minimum reads in a sample for expression filtering (default: 10)",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=2,
        help="Minimum number of samples meeting --min-reads (default: 2)",
    )
    args = parser.parse_args()

    if len(args.inputs) != len(args.outputs):
        raise ValueError("--inputs and --outputs must have the same number of files")
    if args.express_outputs and len(args.inputs) != len(args.express_outputs):
        raise ValueError("--inputs and --express-outputs must have the same number of files")
    if args.mode == "remap" and args.source_column == args.target_column:
        raise ValueError("--source-column and --target-column must be different")

    points_path = Path(args.points)
    if args.mode == "remap":
        position_map = load_position_map(points_path, args.source_column, args.target_column)
    else:
        valid_points = set(load_position_map(points_path, args.source_column, args.source_column).keys())

    for idx, (in_path_str, out_path_str) in enumerate(zip(args.inputs, args.outputs)):
        in_path = Path(in_path_str)
        out_path = Path(out_path_str)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        express_out_path = None
        if args.express_outputs:
            express_out_path = Path(args.express_outputs[idx])
            express_out_path.parent.mkdir(parents=True, exist_ok=True)

        if args.mode == "remap":
            total, kept, express_kept = filter_and_remap_junctions(
                in_path,
                out_path,
                position_map,
                express_out_path,
                args.min_reads,
                args.min_samples,
            )
        else:
            total, kept, express_kept = filter_junctions(
                in_path,
                out_path,
                valid_points,
                express_out_path,
                args.min_reads,
                args.min_samples,
            )
        print(f"{in_path} -> {out_path}: kept {kept}/{total} junctions")
        if express_out_path:
            print(
                f"{in_path} -> {express_out_path}: kept {express_kept}/{kept} "
                f"junctions with >= {args.min_reads} reads in >= {args.min_samples} samples"
            )


if __name__ == "__main__":
    main()
