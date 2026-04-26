#!/usr/bin/env python3
"""Find junctions present in all files with >=10 reads in >=2 samples in >=2 datasets."""

import csv
import sys
from pathlib import Path


def normalize_junction_id(junction_id):
    """Normalize junction by sorting start/end so reversed breakpoints match."""
    parts = str(junction_id).split(":")
    if len(parts) < 3:
        return str(junction_id)

    chrom, start, end = parts[0], parts[1], parts[2]
    tail = parts[3:]

    try:
        start_val = int(start)
        end_val = int(end)
        if start_val <= end_val:
            ordered_start, ordered_end = start, end
        else:
            ordered_start, ordered_end = end, start
    except ValueError:
        if start <= end:
            ordered_start, ordered_end = start, end
        else:
            ordered_start, ordered_end = end, start

    return ":".join([chrom, ordered_start, ordered_end] + tail)


def _parse_count(value):
    """Parse numeric count values, supporting plain numbers and 'count/total'."""
    text = str(value).strip()
    if not text:
        return 0
    if "/" in text:
        text = text.split("/", 1)[0].strip()
    return int(float(text))


def _is_count_like(value):
    """Return True if value can be parsed as a count (or is empty)."""
    text = str(value).strip()
    if not text:
        return True
    try:
        _parse_count(text)
        return True
    except ValueError:
        return False


def load_junctions(file_path):
    """Load junctions and return normalized-key counts plus original labels."""
    junctions = {}
    original_labels = {}
    with open(file_path, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise ValueError(f"Empty file: {file_path}")
        if not header or header[0] != "junction_id":
            raise ValueError(f"First column must be junction_id in file: {file_path}")

        rows = list(reader)

        # Keep only columns that are count-like in every row; this drops metadata
        # columns such as original_junction_id appended at the end.
        numeric_col_indices = []
        for idx in range(1, len(header)):
            is_numeric_col = True
            for row in rows:
                value = row[idx] if idx < len(row) else ""
                if not _is_count_like(value):
                    is_numeric_col = False
                    break
            if is_numeric_col:
                numeric_col_indices.append(idx)

        sample_names = [header[idx] for idx in numeric_col_indices]

        for row in rows:
            if not row:
                continue
            junction_id = row[0]
            normalized_id = normalize_junction_id(junction_id)
            read_counts = []
            for idx in numeric_col_indices:
                value = row[idx] if idx < len(row) else ""
                try:
                    count = _parse_count(value)
                except ValueError:
                    count = 0
                read_counts.append(count)
            if normalized_id not in junctions:
                junctions[normalized_id] = read_counts
                original_labels[normalized_id] = junction_id
    
    return junctions, sample_names, original_labels


def has_min_coverage(read_counts, min_reads=10, min_samples=2):
    """Check if read_counts has at least min_reads in at least min_samples."""
    samples_with_min = sum(1 for count in read_counts if count >= min_reads)
    return samples_with_min >= min_samples


def find_shared_junctions(input_dir, min_reads=10, min_samples=2, output_file=None):
    """Find junctions present in all files with sufficient coverage.
    
    Args:
        input_dir: Path to directory with junction TSV files
        min_reads: Minimum read count (default: 10)
        min_samples: Minimum number of samples with min_reads (default: 2)
        output_file: Optional output file path
    
    Returns:
        Tuple of (shared_junctions dict, all_filenames)
    """
    input_path = Path(input_dir)
    tsv_files = sorted(
        p for p in input_path.glob("GSE*_orthologs_junctions.tsv")
        if "express" not in p.name
    )
    
    if not tsv_files:
        print(f"Error: No TSV files found in {input_dir}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Processing {len(tsv_files)} files:", file=sys.stderr)
    for f in tsv_files:
        print(f"  {f.name}", file=sys.stderr)
    
    # Load all junctions from all files
    all_junctions = {}  # file_name -> {normalized_junction_id -> [read_counts]}
    all_original_labels = {}  # file_name -> {normalized_junction_id -> original junction_id}
    all_sample_names = {}  # file_name -> [sample_names]
    
    for tsv_file in tsv_files:
        print(f"Loading {tsv_file.name}...", file=sys.stderr)
        junctions, samples, original_labels = load_junctions(tsv_file)
        all_junctions[tsv_file.name] = junctions
        all_original_labels[tsv_file.name] = original_labels
        all_sample_names[tsv_file.name] = samples
        print(f"  {len(junctions)} junctions, {len(samples)} samples", file=sys.stderr)
    
    # Find normalized junctions present in all files
    file_names = list(all_junctions.keys())
    reference_file = next((name for name in file_names if name.startswith("GSE115736")), file_names[0])
    first_file_junctions = set(all_junctions[reference_file].keys())
    
    shared_junction_keys = first_file_junctions.copy()
    for fname in file_names:
        if fname == reference_file:
            continue
        shared_junction_keys &= set(all_junctions[fname].keys())
    
    print(f"\nJunctions in all {len(file_names)} files: {len(shared_junction_keys)}", file=sys.stderr)
    
    # Filter by coverage criterion: >= min_reads in >= min_samples in at least two files
    filtered_shared = {}
    for junction_key in shared_junction_keys:
        # Count how many datasets pass the per-dataset coverage threshold.
        passing_dataset_count = 0
        for fname in file_names:
            read_counts = all_junctions[fname][junction_key]
            if has_min_coverage(read_counts, min_reads, min_samples):
                passing_dataset_count += 1
        
        if passing_dataset_count >= 2:
            reference_junction_id = all_original_labels[reference_file].get(junction_key, junction_key)
            # Store the junction with all its data from all files
            filtered_shared[reference_junction_id] = {
                fname: all_junctions[fname][junction_key]
                for fname in file_names
            }
    
    print(
        f"Shared junctions with >={min_reads} reads in >={min_samples} samples in >=2 datasets: "
        f"{len(filtered_shared)}",
          file=sys.stderr)
    
    # Write output
    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, "w", newline="") as fh:
            # Write header
            all_samples = []
            for fname in file_names:
                all_samples.extend(all_sample_names[fname])
            
            writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
            header = ["junction_id"] + [
                                        f"{fname.replace('.tsv', '').replace('_orthologs_junctions', '')}_{s}"
                                        for fname in file_names
                                        for s in all_sample_names[fname]]
            writer.writerow(header)
            
            # Write junctions
            for junction_id in sorted(filtered_shared.keys()):
                row = [junction_id]
                for fname in file_names:
                    row.extend(filtered_shared[junction_id][fname])
                writer.writerow(row)
        
        print(f"Output written to {output_path}", file=sys.stderr)
    
    return filtered_shared, file_names


def main():
    input_dir = "/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026"
    output_file = "/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/shared_junctions_filtered.tsv"
    
    min_reads = 10  # Minimum read count
    min_samples = 2  # Minimum number of samples
    
    shared_junctions, files = find_shared_junctions(
        input_dir,
        min_reads=min_reads,
        min_samples=min_samples,
        output_file=output_file
    )
    return

if __name__ == "__main__":
    main()