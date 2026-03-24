#!/usr/bin/env python3
"""Assign clusters to shared junctions based on overlapping breakpoints."""

import networkx as nx
import pandas as pd

def assign_clusters_to_groups(df, junctionCol):    
    """Assign junctions to clusters based on shared breakpoints."""
    print("Assigning clusters...")
    df = df.reset_index(drop=False)
    
    # Extract chromosome part
    df.loc[:, 'chrom'] = df[junctionCol].astype(str).str.split(':').str[0]
    
    results = []
    cluster_counter = 1  # Global cluster counter    
    
    # Process each chromosome separately
    # Only regular chromosomes (length <= 5) to avoid scaffolds and unplaced contigs
    for chrom in sorted(df[df['chrom'].str.len() <= 5]['chrom'].unique()):
        print(f"Processing chromosome {chrom}...")
        df_chrom = df[df['chrom'] == chrom].copy()      
        
        # Build graph: connect junctions that share breakpoints
        G = nx.Graph()
        for idx, row in df_chrom.iterrows():
            parts = row[junctionCol].split(':')
            start, end = parts[1], parts[2]
            # Create nodes for breakpoints and connect to junctions
            G.add_edge(idx, f'start_{start}')
            G.add_edge(idx, f'end_{end}')
        
        # Find connected components (groups of junctions with shared breakpoints)
        components = list(nx.connected_components(G))
        print(f"  {len(components)} initial components")
       
        for comp in components:
            # Extract only junction nodes (filter out breakpoint nodes)
            junction_nodes = sorted([node for node in comp if node in df_chrom.index])
            if not junction_nodes:
                continue
            
            df_cluster = df_chrom.loc[junction_nodes].copy()
            
            # Assign cluster ID
            df_cluster['cluster'] = f'clu_{cluster_counter}'
            cluster_counter += 1
            results.append(df_cluster)
    
    # Combine all results
    if not results:
        print("No clusters found!")
        return pd.DataFrame(), pd.DataFrame()
    
    df_result = pd.concat(results).sort_index()
    df_result.drop(columns=['chrom'], inplace=True)
    
    print(f"Total junctions: {len(df_result)}")
    print(f"Total clusters: {df_result['cluster'].nunique()}")
    
    # Keep only clusters with more than 1 junction (alternative splicing events)
    df_result_AS = df_result[df_result['cluster'].duplicated(keep=False)].copy()
    print(f"AS clusters (>1 junction): {df_result_AS['cluster'].nunique()}")

    df_result.set_index(junctionCol, inplace=True, drop=True)
    df_result_AS.set_index(junctionCol, inplace=True, drop=True)
    
    return df_result, df_result_AS

def main():
    """Main entry point."""
    ensembl_dir = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/'
    
    # Load pre-filtered shared junctions
    input_file = ensembl_dir + "sharedJunctions_2026/shared_junctions_filtered.tsv"
    print(f"Loading {input_file}...")
    merged_df = pd.read_csv(input_file, sep="\t", index_col='junction_id')
    
    print(f"Loaded {len(merged_df)} junctions with {len(merged_df.columns)} samples\n")
    
    # Run clustering
    junction_sum_final, junction_AS = assign_clusters_to_groups(merged_df, 'junction_id')
    
    if junction_AS.empty:
        print("No alternative splicing clusters found!")
        return
    
    # Sort by cluster, then by junction_id
    junction_AS = junction_AS.reset_index().sort_values(['cluster', 'junction_id']).set_index('junction_id')
    
    # Write outputs
    output_file = ensembl_dir + "sharedJunctions_2026/shared_junctions_cluster_all.tsv"
    junction_sum_final.to_csv(output_file, sep="\t")
    print(f"\nWrote all clusters to: {output_file}")
    
    output_file_AS = ensembl_dir + "sharedJunctions_2026/shared_junctions_cluster_AS.tsv"
    junction_AS.to_csv(output_file_AS, sep="\t")
    print(f"Wrote AS junctions to: {output_file_AS}")
    
    num_clusters = junction_AS['cluster'].nunique()
    num_as_junctions = len(junction_AS)
    print(f"\nSummary: {num_clusters} unique clusters, {num_as_junctions} junctions in AS clusters")

if __name__ == "__main__":
    main()
