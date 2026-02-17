#!/usr/bin/env python3
"""
Analyze shared junctions across immune cell types from leafcutter results.
Creates bar plots showing the number of clusters (defined by junctions) 
that are shared across 1, 2, 3, 4, or 5 immune cell types.

Two analyses:
1. All junctions
2. Junctions with |ΔΨ| >= 0.2 AND p-value <= 0.05
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from collections import defaultdict

def extract_junction(junction_str):
    """
    Extract junction coordinates from a junction string.
    Input: 'chr1:944800:945057:clu_1'
    Output: 'chr1:944800:945057'
    """
    parts = junction_str.split(':')
    if len(parts) >= 3:
        # Return chr:start:end (without the cluster name)
        return f"{parts[0]}:{parts[1]}:{parts[2]}"
    return junction_str

def read_leafcutter_clusters(file_path):
    """
    Read leafcutter AS_clusters file and extract junction coordinates.
    Returns a set of unique junctions.
    """
    junctions = set()
    try:
        # Read the file, skipping the header row
        df = pd.read_csv(file_path, sep='\t', skiprows=1, header=None)
        # First column contains the junction:cluster info
        for junction_str in df[0]:
            junction = extract_junction(str(junction_str))
            junctions.add(junction)
        print(f"  Found {len(junctions)} unique junctions")
        return junctions
    except Exception as e:
        print(f"  Error reading {file_path}: {e}")
        return set()

def read_leafcutter_clusters_filtered(effect_sizes_path, significance_path, deltapsi_threshold=0.2, pval_threshold=0.05):
    """
    Read leafcutter effect_sizes and cluster significance files, and extract junction coordinates 
    filtered by |deltapsi| >= deltapsi_threshold and p-value <= pval_threshold.
    Returns a set of unique junctions.
    """
    junctions = set()
    try:
        # Read the effect_sizes file
        df_effect = pd.read_csv(effect_sizes_path, sep='\t')
        
        # Read the cluster significance file
        df_sig = pd.read_csv(significance_path, sep='\t')
        
        # Extract cluster ID from intron string: chr:start:end:clu_X -> chr:clu_X
        def extract_cluster_id(intron):
            parts = str(intron).split(':')
            if len(parts) >= 4:
                return f"{parts[0]}:{parts[3]}"
            return None
        
        df_effect['cluster'] = df_effect['intron'].apply(extract_cluster_id)
        
        # Merge on cluster
        df_merged = df_effect.merge(df_sig[['cluster', 'p']], on='cluster', how='left')
        
        # Filter by deltapsi and p-value thresholds
        df_filtered = df_merged[(abs(df_merged['deltapsi']) >= deltapsi_threshold) & 
                                (df_merged['p'] <= pval_threshold)]
        
        print(f"  Total junctions: {len(df_effect)}")
        print(f"  Filtered junctions (|deltapsi| >= {deltapsi_threshold} AND p-value <= {pval_threshold}): {len(df_filtered)}")
        
        # Extract junctions
        for junction_str in df_filtered['intron']:
            junction = extract_junction(str(junction_str))
            junctions.add(junction)
        
        print(f"  Unique junctions: {len(junctions)}")
        return junctions
    except Exception as e:
        print(f"  Error reading files: {e}")
        return set()

def analyze_shared_clusters(base_dir, cell_types, use_deltapsi=False, deltapsi_threshold=0.2, pval_threshold=0.05):
    """
    Analyze shared junctions across immune cell types.
    
    Parameters:
    -----------
    base_dir : str
        Base directory containing cell type folders
    cell_types : list
        List of cell type folder names
    use_deltapsi : bool
        If True, filter by deltapsi threshold
    deltapsi_threshold : float
        Minimum absolute deltapsi value to include (if use_deltapsi=True)
    pval_threshold : float
        Maximum p-value to include for significant clusters
    
    Returns:
    --------
    dict : Counts of junctions shared by N cell types (N = 1 to 5)
    """
    cell_junctions = {}
    
    if use_deltapsi:
        print(f"Reading leafcutter results (|deltapsi| >= {deltapsi_threshold} AND p-value <= {pval_threshold})...\n")
    else:
        print("Reading leafcutter results for each cell type...\n")
    
    for cell_type in cell_types:
        print(f"Processing {cell_type}...")
        cell_dir = Path(base_dir) / cell_type
        
        if use_deltapsi:
            # Use effect_sizes and significance files for deltapsi and p-value filtering
            effect_sizes_file = cell_dir / "leafcutter_ds_effect_sizes.txt"
            significance_file = cell_dir / "leafcutter_ds_cluster_significance.txt"
            
            if not effect_sizes_file.exists():
                file_list = sorted(cell_dir.glob("leafcutter_ds_effect_sizes*.txt"), reverse=True)
                if file_list:
                    effect_sizes_file = file_list[0]
            
            if not significance_file.exists():
                file_list = sorted(cell_dir.glob("leafcutter_ds_cluster_significanc*.txt"), reverse=True)
                if file_list:
                    significance_file = file_list[0]
            
            if not effect_sizes_file.exists() or not significance_file.exists():
                print(f"  Warning: Missing files in {cell_dir}")
                continue
            
            junctions = read_leafcutter_clusters_filtered(effect_sizes_file, significance_file, 
                                                         deltapsi_threshold, pval_threshold)
        else:
            # Use AS_clusters file for all junctions
            file_path = cell_dir / "AS_clusters_psi_HN6.txt"
            if not file_path.exists():
                file_list = sorted(cell_dir.glob("AS_clusters_psi_HN*.txt"), reverse=True)
                if file_list:
                    file_path = file_list[0]
                else:
                    print(f"  Warning: No AS_clusters file found in {cell_dir}")
                    continue
            junctions = read_leafcutter_clusters(file_path)
        
        cell_junctions[cell_type] = junctions
    
    # Count sharing across cell types
    print("\n" + "="*60)
    print("Analyzing junction sharing...")
    print("="*60)
    
    # Get all unique junctions
    all_junctions = set()
    for junctions in cell_junctions.values():
        all_junctions.update(junctions)
    
    print(f"Total unique junctions across all cell types: {len(all_junctions)}")
    
    # Count how many cell types each junction appears in
    sharing_counts = defaultdict(int)
    for junction in all_junctions:
        count = sum(1 for junctions in cell_junctions.values() if junction in junctions)
        sharing_counts[count] += 1
    
    # Print summary
    print("\nSharing summary:")
    for num_cell_types in sorted(sharing_counts.keys()):
        print(f"  Junctions in {num_cell_types} cell type(s): {sharing_counts[num_cell_types]}")
    
    return dict(sharing_counts)

def plot_sharing_results(sharing_counts, cell_types, output_file, title_suffix=""):
    """
    Create a bar plot showing the number of clusters shared across cell types.
    
    Parameters:
    -----------
    sharing_counts : dict
        Dictionary with keys = number of cell types, values = number of junctions
    cell_types : list
        List of cell type names (for annotation)
    output_file : str
        Output file path for the plot
    title_suffix : str
        Additional text to append to the title
    """
    n_cell_types = len(cell_types)
    
    # Ensure all categories are present (fill with 0 if missing)
    counts = [sharing_counts.get(i, 0) for i in range(1, n_cell_types + 1)]
    x_labels = [str(i) for i in range(1, n_cell_types + 1)]
    x_pos = np.arange(len(x_labels))
    
    # Create the bar plot
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(x_pos, counts, color='steelblue', edgecolor='navy', alpha=0.7)
    
    # Add value labels on top of bars
    for i, (bar, count) in enumerate(zip(bars, counts)):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(count)}',
                ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    # Customize the plot
    ax.set_xlabel('Number of Immune Cell Types', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Shared Clusters', fontsize=12, fontweight='bold')
    ax.set_title(f'Shared Junction Clusters{title_suffix}', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # Add cell type names in legend or subtitle
    cell_types_str = ', '.join(cell_types)
    fig.text(0.5, 0.02, f'Cell types: {cell_types_str}', 
             ha='center', fontsize=10, style='italic')
    
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")
    plt.close()

def main():
    """Main function."""
    base_dir = "/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/leafcutter_GSE115736_GSE116177"
    cell_types = ["CD4T", "CD8T", "Mono", "Neut", "NveB"]
    
    print("="*70)
    print("Analyzing Shared Clusters Across Immune Cell Types")
    print("="*70)
    print(f"Base directory: {base_dir}")
    print(f"Cell types: {', '.join(cell_types)}")
    print("="*70)
    
    # Analysis 1: All junctions
    print("\n" + "="*70)
    print("ANALYSIS 1: All Junctions")
    print("="*70 + "\n")
    sharing_counts_all = analyze_shared_clusters(base_dir, cell_types, use_deltapsi=False)
    output_file_all = os.path.join(base_dir, "shared_clusters_across_immune_cells.png")
    plot_sharing_results(sharing_counts_all, cell_types, output_file_all)
    
    # Analysis 2: Filtered by deltapsi and p-value
    print("\n" + "="*70)
    print("ANALYSIS 2: Junctions with |ΔΨ| >= 0.2 AND p-value <= 0.05")
    print("="*70 + "\n")
    deltapsi_threshold = 0.2
    pval_threshold = 0.05
    sharing_counts_filtered = analyze_shared_clusters(base_dir, cell_types, 
                                                      use_deltapsi=True, 
                                                      deltapsi_threshold=deltapsi_threshold,
                                                      pval_threshold=pval_threshold)
    output_file_filtered = os.path.join(base_dir, "shared_clusters_across_immune_cells_deltapsi_0.2_pval_0.05.png")
    plot_sharing_results(sharing_counts_filtered, cell_types, output_file_filtered, 
                        f" (|ΔΨ| ≥ {deltapsi_threshold}, p ≤ {pval_threshold})")
    
    print("\n" + "="*70)
    print("All analyses complete!")
    print("="*70)

if __name__ == "__main__":
    main()
