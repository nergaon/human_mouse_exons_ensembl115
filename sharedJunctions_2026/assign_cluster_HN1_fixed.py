import networkx as nx
import pandas as pd

def assign_clusters_to_groups(df, junctionCol, startNumCol, min_psi=0.0):    
    print("Assigning clusters...")
    df = df.reset_index(drop=False)
    # Extract chromosome part
    df.loc[:, 'chrom'] = df[junctionCol].astype(str).str.split(':').str[0]
    # Sample columns: all columns from startNumCol onwards, excluding 'chrom'
    sample_cols = [col for col in df.columns[startNumCol:] if col != 'chrom' and col != junctionCol]
    print(f"Sample columns: {len(sample_cols)} columns")
    
    # Initialize an empty list to store results
    results = []
    cluster_counter = 1  # Global cluster counter    
    # Process each chromosome separately. only regular chromosomes (length <= 5) are processed, to avoid issues with scaffolds and unplaced contigs
    for chrom in df[df['chrom'].str.len() <= 5]['chrom'].unique():
        print(f"Processing chromosome {chrom}...")
        # Filter rows for the current chromosome
        df_chrom = df[df['chrom'] == chrom].copy()      
        # Step 1: Build initial graph and clusters
        G = nx.Graph()
        # Add edges between index and their start and end positions
        for idx, row in df_chrom.iterrows():
            parts = row[junctionCol].split(':')
            start, end = parts[1], parts[2]
            G.add_edge(idx, f'start_{start}')
            G.add_edge(idx, f'end_{end}')
        # Find connected components
        components = list(nx.connected_components(G))
        print(f"  {len(components)} initial components")
       
        for comp in components:
            # Get only junction nodes (i.e., DataFrame indices)
            junction_nodes = [node for node in comp if node in df_chrom.index]
            df_cluster = df_chrom.loc[junction_nodes].copy()
    
            # Step 2: Filter out low-expression junctions (<min_psi of total in cluster in all samples)
            # Calculate total per sample (column) within the cluster
            sample_totals = df_cluster[sample_cols].sum()
            # Create a boolean DataFrame: True if a junction is ≥ min_psi of the sample total
            high_expr_mask = (df_cluster[sample_cols] >= (min_psi * sample_totals))
            # Keep rows (junctions) with at least one True in any sample
            keep_mask = high_expr_mask.any(axis=1)
            df_cluster_filtered = df_cluster[keep_mask].copy()
            if df_cluster_filtered.empty:
                continue  # Skip if nothing left
    
            # Step 3: Rebuild graph for filtered cluster to check for disconnection
            G_filtered = nx.Graph()
            for idx, row in df_cluster_filtered.iterrows():
                parts = row[junctionCol].split(':')
                start, end = parts[1], parts[2]
                G_filtered.add_edge(idx, f'start_{start}')
                G_filtered.add_edge(idx, f'end_{end}')
    
            subcomponents = list(nx.connected_components(G_filtered))
    
            for subcomp in subcomponents:
                sub_junctions = [node for node in subcomp if node in df_cluster_filtered.index]
                if not sub_junctions:
                    continue
                df_subcluster = df_cluster_filtered.loc[sub_junctions].copy()
                df_subcluster['cluster'] = f'clu_{cluster_counter}'
                cluster_counter += 1
                results.append(df_subcluster)
    
    # Combine all results
    if not results:
        print("No clusters found!")
        return pd.DataFrame(), pd.DataFrame()
    
    df_result = pd.concat(results).sort_index()
    df_result.drop(columns=['chrom'], inplace=True)
    print(f"Total clusters created: {df_result['cluster'].nunique()}")
    
    # Keep only clusters with more than 1 junction
    df_result_AS = df_result[df_result['cluster'].duplicated(keep=False)]
    print(f"AS clusters (>1 junction): {df_result_AS['cluster'].nunique()}")

    df_result.set_index(junctionCol, inplace=True, drop=True)
    df_result_AS.set_index(junctionCol, inplace=True, drop=True)
    return df_result, df_result_AS

def main():
    """Assign clusters to shared junctions from the filtered junction file."""
    ensembl_dir = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/'
    
    # Load the pre-filtered shared junctions
    input_file = ensembl_dir + "sharedJunctions_2026/shared_junctions_filtered.tsv"
    print(f"Loading {input_file}...")
    merged_df = pd.read_csv(input_file, sep="\t", index_col='junction_id')
    
    print(f"Loaded {len(merged_df)} junctions")
    print(f"Junctions shape: {merged_df.shape}")
    print(f"Sample columns: {merged_df.columns.tolist()[:5]}... (showing first 5)")
    
    # All columns are sample data (no metadata columns to skip)
    startNumCol = 0
    
    # Run clustering on junction_id
    junction_sum_final, junction_AS = assign_clusters_to_groups(
        merged_df, 
        'junction_id',  # Will be reset_index'd, so this becomes a column
        startNumCol,
        min_psi=0.0
    )
    
    # Sort by cluster, then by junction_id
    junction_AS = junction_AS.reset_index().sort_values(['cluster', 'junction_id']).set_index('junction_id')
    
    # Write outputs
    output_file = ensembl_dir + "sharedJunctions_2026/shared_junctions_cluster_all.tsv"
    junction_sum_final.to_csv(output_file, sep="\t")
    print(f"Wrote all clusters to {output_file}")
    
    output_file_AS = ensembl_dir + "sharedJunctions_2026/shared_junctions_cluster_AS.tsv"
    junction_AS.to_csv(output_file_AS, sep="\t")
    print(f"Wrote AS junctions (clusters with >1 junction) to {output_file_AS}")
    
    num_clusters = junction_AS['cluster'].nunique()
    num_as_junctions = len(junction_AS)
    print(f"\n{num_clusters} unique clusters, {num_as_junctions} AS (alternative splicing) junctions")
    return
        
if __name__ == "__main__":
    main()
