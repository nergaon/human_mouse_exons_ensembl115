#!/usr/bin/env python3
"""Find junctions present in all files with >=10 reads in at least 2 samples each."""

import csv
import sys
from pathlib import Path
from collections import defaultdict


def load_junctions(file_path):
    """Load junctions from a TSV file and return dict: junction_id -> [read_counts]."""
    junctions = {}
    with open(file_path, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty file: {file_path}")
        
        sample_names = reader.fieldnames[1:]  # Skip junction_id column
        
        for row in reader:
            junction_id = row["junction_id"]
            read_counts = []
            for sample in sample_names:
                try:
                    count = int(row[sample])
                except (ValueError, KeyError):
                    count = 0
                read_counts.append(count)
            junctions[junction_id] = read_counts
    
    return junctions, sample_names


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
    tsv_files = sorted(input_path.glob("GSE*.tsv"))
    
    if not tsv_files:
        print(f"Error: No TSV files found in {input_dir}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Processing {len(tsv_files)} files:", file=sys.stderr)
    for f in tsv_files:
        print(f"  {f.name}", file=sys.stderr)
    
    # Load all junctions from all files
    all_junctions = {}  # file_name -> {junction_id -> [read_counts]}
    all_sample_names = {}  # file_name -> [sample_names]
    
    for tsv_file in tsv_files:
        print(f"Loading {tsv_file.name}...", file=sys.stderr)
        junctions, samples = load_junctions(tsv_file)
        all_junctions[tsv_file.name] = junctions
        all_sample_names[tsv_file.name] = samples
        print(f"  {len(junctions)} junctions, {len(samples)} samples", file=sys.stderr)
    
    # Find junctions present in all files
    file_names = list(all_junctions.keys())
    first_file_junctions = set(all_junctions[file_names[0]].keys())
    
    shared_junction_ids = first_file_junctions.copy()
    for fname in file_names[1:]:
        shared_junction_ids &= set(all_junctions[fname].keys())
    
    print(f"\nJunctions in all {len(file_names)} files: {len(shared_junction_ids)}", file=sys.stderr)
    
    # Filter by coverage criterion: >= min_reads in >= min_samples in at least one file
    filtered_shared = {}
    for junction_id in shared_junction_ids:
        # Check if this junction passes the coverage threshold in any file
        passes = False
        for fname in file_names:
            read_counts = all_junctions[fname][junction_id]
            if has_min_coverage(read_counts, min_reads, min_samples):
                passes = True
                break
        
        if passes:
            # Store the junction with all its data from all files
            filtered_shared[junction_id] = {
                fname: all_junctions[fname][junction_id]
                for fname in file_names
            }
    
    print(f"Shared junctions with >={min_reads} reads in >={min_samples} samples: {len(filtered_shared)}", 
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
            header = ["junction_id"] + [f"{fname.replace('.tsv', '')}_{s}" 
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