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

# Set matplotlib backend to non-interactive
plt.switch_backend('Agg')

def extract_cluster(junction_str):
    """
    Extract the cluster identifier from a junction string.
    Input: 'chr1:944800:945057:clu_1'
    Output: 'chr1:944800:945057:clu_1' (full cluster string)
    This keeps coordinates and cluster name together so that each cluster
    is treated as a unique entity. We avoid stripping off the cluster
    since we want to count shared *clusters* rather than individual
    junction coordinates.
    """
    # In many files the cluster string is already in the format we want
    # so simply return the input unchanged (ensures consistent behaviour)
    return str(junction_str)

def read_leafcutter_clusters(file_path):
    """
    Read leafcutter AS_clusters file and extract cluster identifiers.
    Returns a set of unique clusters.
    """
    clusters = set()
    try:
        # Read the file, skipping the header row
        df = pd.read_csv(file_path, sep='\t', skiprows=1, header=None)
        # First column contains the junction:cluster info
        for junction_str in df[0]:
            cluster = extract_cluster(str(junction_str))
            clusters.add(cluster)
        print(f"  Found {len(clusters)} unique clusters")
        return clusters
    except Exception as e:
        print(f"  Error reading {file_path}: {e}")
        return set()

def read_leafcutter_clusters_filtered(effect_sizes_path, significance_path, deltapsi_threshold=0.2, pval_threshold=0.05):
    """
    Read leafcutter effect_sizes and cluster significance files, and extract cluster identifiers
    filtered by |deltapsi| >= deltapsi_threshold and p-value <= pval_threshold.
    Returns a set of unique clusters.
    """
    clusters = set()
    try:
        # Read the effect_sizes file
        df_effect = pd.read_csv(effect_sizes_path, sep='\t')
        
        # Read the cluster significance file
        df_sig = pd.read_csv(significance_path, sep='\t')
        
        # Extract cluster ID from intron string: chr:start:end:clu_X -> chr:clu_X
        def extract_cluster_id(intron):
            parts = str(intron).split(':')
            if len(parts) >= 4:
                return f"{parts[0]}:{parts[-1]}"  # chr:clu_X
            return str(intron)
        
        df_effect['cluster'] = df_effect['intron'].apply(extract_cluster_id)
        
        # Merge on cluster
        df_merged = df_effect.merge(df_sig[['cluster', 'p.adjust']], on='cluster', how='left')
        
        # Filter by deltapsi and p-value thresholds
        df_filtered = df_merged[(abs(df_merged['deltapsi']) >= deltapsi_threshold) & 
                                (df_merged['p.adjust'] <= pval_threshold)]
        
        print(f"  Total rows: {len(df_effect)}")
        print(f"  Filtered rows (|deltapsi| >= {deltapsi_threshold} AND p-value <= {pval_threshold}): {len(df_filtered)}")
        
        # Extract clusters
        for cluster_str in df_filtered['cluster']:
            clusters.add(cluster_str)
        
        print(f"  Unique clusters: {len(clusters)}")
        return clusters
    except Exception as e:
        print(f"  Error reading files: {e}")
        return set()

def analyze_shared_clusters(base_dir, cell_types, use_deltapsi=False, deltapsi_threshold=0.2, pval_threshold=0.05):
    """
    Analyze shared clusters across immune cell types.
    
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
    dict : Counts of clusters shared by N cell types (N = 1 to 5)
    """
    cell_clusters = {}
    
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
            
            clusters = read_leafcutter_clusters_filtered(effect_sizes_file, significance_file, 
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
            clusters = read_leafcutter_clusters(file_path)
        
        cell_clusters[cell_type] = clusters
    
    # Count sharing across cell types
    print("\n" + "="*60)
    print("Analyzing cluster sharing...")
    print("="*60)
    
    # Get all unique clusters
    all_clusters = set()
    for clusters in cell_clusters.values():
        all_clusters.update(clusters)
    
    print(f"Total unique clusters across all cell types: {len(all_clusters)}")
    
    # Count how many cell types each cluster appears in
    sharing_counts = defaultdict(int)
    cluster_names_by_sharing = defaultdict(list)
    for cluster in all_clusters:
        count = sum(1 for clusters in cell_clusters.values() if cluster in clusters)
        sharing_counts[count] += 1
        cluster_names_by_sharing[count].append(cluster)
    
    # Print summary
    print("\nSharing summary:")
    for num_cell_types in sorted(sharing_counts.keys()):
        print(f"  Clusters in {num_cell_types} cell type(s): {sharing_counts[num_cell_types]}")
    
    return dict(sharing_counts), dict(cluster_names_by_sharing)

def plot_sharing_results_combined(sharing_counts_all, sharing_counts_filtered, cell_types, output_file, deltapsi_threshold, pval_threshold):
    """
    Create a grouped bar plot showing both analyses side by side.
    
    Parameters:
    -----------
    sharing_counts_all : dict
        Dictionary with sharing counts for all clusters
    sharing_counts_filtered : dict
        Dictionary with sharing counts for filtered clusters
    cell_types : list
        List of cell type names (for annotation)
    output_file : str
        Output file path for the plot
    deltapsi_threshold : float
        Delta PSI threshold used for filtering
    pval_threshold : float
        P-value threshold used for filtering
    """
    n_cell_types = len(cell_types)
    
    # Ensure all categories are present (fill with 0 if missing)
    counts_all = [sharing_counts_all.get(i, 0) for i in range(1, n_cell_types + 1)]
    counts_filtered = [sharing_counts_filtered.get(i, 0) for i in range(1, n_cell_types + 1)]
    x_labels = [str(i) for i in range(1, n_cell_types + 1)]
    x_pos = np.arange(len(x_labels))
    
    # Width of each bar
    bar_width = 0.35
    
    # Create the grouped bar plot
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Bars for all clusters
    bars_all = ax.bar(x_pos - bar_width/2, counts_all, bar_width, 
                      color='steelblue', edgecolor='navy', alpha=0.7, label='All Clusters')
    
    # Bars for filtered clusters
    bars_filtered = ax.bar(x_pos + bar_width/2, counts_filtered, bar_width,
                          color='darkorange', edgecolor='darkred', alpha=0.7, 
                          label=f'|ΔΨ| ≥ {deltapsi_threshold}, p ≤ {pval_threshold}')
    
    # Add value labels on top of bars
    for bars, counts in [(bars_all, counts_all), (bars_filtered, counts_filtered)]:
        for i, (bar, count) in enumerate(zip(bars, counts)):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(count)}',
                    ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Customize the plot
    ax.set_xlabel('Number of Immune Cell Types', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Shared Clusters', fontsize=12, fontweight='bold')
    ax.set_title('Shared Clusters Across Immune Cell Types', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # Add legend
    ax.legend(loc='upper right', fontsize=11)
    
    # Add cell type names in subtitle
    cell_types_str = ', '.join(cell_types)
    fig.text(0.5, 0.02, f'Cell types: {cell_types_str}', 
             ha='center', fontsize=10, style='italic')
    
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nCombined plot saved to: {output_file}")
    plt.close()

def save_cluster_names_to_file(cluster_names_all, cluster_names_filtered, base_dir, cell_types, deltapsi_threshold, pval_threshold):
    """
    Save cluster names grouped by sharing level to a text file.
    Each column represents clusters found in N cell types, with each cluster on a separate row.
    
    Parameters:
    -----------
    cluster_names_all : dict
        Dictionary with cluster names for all clusters by sharing level
    cluster_names_filtered : dict
        Dictionary with cluster names for filtered clusters by sharing level
    base_dir : str
        Base directory to save the file
    cell_types : list
        List of cell type names
    deltapsi_threshold : float
        Delta PSI threshold used for filtering
    pval_threshold : float
        P-value threshold used for filtering
    """
    output_file = os.path.join(base_dir, "cluster_names_by_sharing_level.txt")
    
    # Extract cluster IDs (clu_XXXX) from full cluster names
    def extract_cluster_id(cluster_name):
        # Handle both formats: full junction string or just cluster ID
        if ':' in cluster_name:
            return cluster_name.split(':')[-1]  # Extract clu_XXXX from end
        return cluster_name
    
    # Process cluster names - extract IDs and organize by sharing level
    processed_all = {}
    processed_filtered = {}
    
    for sharing_level in cluster_names_all:
        processed_all[sharing_level] = sorted([extract_cluster_id(name) for name in cluster_names_all[sharing_level]])
    
    for sharing_level in cluster_names_filtered:
        processed_filtered[sharing_level] = sorted([extract_cluster_id(name) for name in cluster_names_filtered[sharing_level]])
    
    # Find the maximum number of cell types
    max_sharing = max(max(processed_all.keys(), default=0), max(processed_filtered.keys(), default=0))
    
    with open(output_file, 'w') as f:
        f.write("Cluster names by sharing level (each cluster on separate row)\n")
        f.write("="*80 + "\n\n")
        
        # Write header - one column for each sharing level
        header_parts = []
        for level in range(1, max_sharing + 1):
            header_parts.append(f"{level}_cells_all")
        for level in range(1, max_sharing + 1):
            header_parts.append(f"{level}_cells_filtered")
        header = "\t".join(header_parts)
        f.write(header + "\n")
        f.write("-" * len(header) + "\n")
        
        # Find the maximum number of rows we need (longest column)
        max_rows = max(
            max(len(processed_all.get(level, [])) for level in range(1, max_sharing + 1)),
            max(len(processed_filtered.get(level, [])) for level in range(1, max_sharing + 1))
        )
        
        # Write data row by row
        for row_idx in range(max_rows):
            row_parts = []
            
            # All clusters columns (1 cell, 2 cells, etc.)
            for level in range(1, max_sharing + 1):
                clusters_at_level = processed_all.get(level, [])
                if row_idx < len(clusters_at_level):
                    row_parts.append(clusters_at_level[row_idx])
                else:
                    row_parts.append("")
            
            # Filtered clusters columns (1 cell, 2 cells, etc.)
            for level in range(1, max_sharing + 1):
                clusters_at_level = processed_filtered.get(level, [])
                if row_idx < len(clusters_at_level):
                    row_parts.append(clusters_at_level[row_idx])
                else:
                    row_parts.append("")
            
            f.write("\t".join(row_parts) + "\n")
    
    print(f"\nCluster names saved to: {output_file}")

def main():
    base_dir = "/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/leafcutter_GSE115736_GSE116177"
    cell_types = ["CD4T", "CD8T", "NveB", "NK", "Mono", "Neut"]
    
    print("="*70)
    print("Analyzing Shared Clusters Across Immune Cell Types")
    print("="*70)
    print(f"Base directory: {base_dir}")
    print(f"Cell types: {', '.join(cell_types)}")
    print("="*70)
    
    # Analysis 1: All junctions
    print("\n" + "="*70)
    print("ANALYSIS 1: All Clusters")
    print("="*70 + "\n")
    sharing_counts_all, cluster_names_all = analyze_shared_clusters(base_dir, cell_types, use_deltapsi=False)
    
    # Analysis 2: Filtered by deltapsi and p-value
    print("\n" + "="*70)
    print("ANALYSIS 2: Clusters with |ΔΨ| >= 0.2 AND p-value <= 0.05")
    print("="*70 + "\n")
    deltapsi_threshold = 0.2
    pval_threshold = 0.05
    sharing_counts_filtered, cluster_names_filtered = analyze_shared_clusters(base_dir, cell_types, 
                                                      use_deltapsi=True, 
                                                      deltapsi_threshold=deltapsi_threshold,
                                                      pval_threshold=pval_threshold)
    
    # Save cluster names to text file
    save_cluster_names_to_file(cluster_names_all, cluster_names_filtered, base_dir, cell_types, deltapsi_threshold, pval_threshold)
    
    # Create combined plot
    output_file_combined = os.path.join(base_dir, "shared_clusters_across_immune_cells_combined.png")
    plot_sharing_results_combined(sharing_counts_all, sharing_counts_filtered, cell_types, 
                                 output_file_combined, deltapsi_threshold, pval_threshold)
    
    print("\n" + "="*70)
    print("All analyses complete!")
    print("="*70)

if __name__ == "__main__":
    main()
