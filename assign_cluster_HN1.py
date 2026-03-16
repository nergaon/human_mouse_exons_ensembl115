import networkx as nx
import pandas as pd

def assign_clusters_to_groups(df, junctionCol, startNumCol, min_psi=0.0):    
    print("Assigning cluster")
    df = df.reset_index(drop=False)
    # Extract chromosome part
    df.loc[:, 'chrom'] = df[junctionCol].astype(str).str.split(':').str[0]
    # Sample columns (exclude metadata columns)
    sample_cols = df.columns[startNumCol:-1]
    print(sample_cols)
    #sample_cols = [col for col in df.columns if col not in [junctionCol, 'chrom']]
    # Initialize an empty list to store results
    results = []
    cluster_counter = 1  # Global cluster counter    
    # Process each chromosome separately
    for chrom in df['chrom'].unique():
        print(chrom)
        # Filter rows for the current chromosome
        df_chrom = df[df['chrom'] == chrom].copy()      
        # Step 1: Build initial graph and clusters
        G = nx.Graph()
        # Add edges between index and their start and end positions
        for idx, row in df_chrom.iterrows():
            _, start, end = row[junctionCol].split(':')
            G.add_edge(idx, f'start_{start}')
            G.add_edge(idx, f'end_{end}')
        # Find connected components
        components = list(nx.connected_components(G))
       
        for comp in components:
            # Get only junction nodes (i.e., DataFrame indices)
            junction_nodes = [node for node in comp if node in df_chrom.index]
            df_cluster = df_chrom.loc[junction_nodes].copy()
    
            # Step 2: Filter out low-expression junctions (<5% of total in cluster in all sample)
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
                _, start, end = row[junctionCol].split(':')
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
    df_result = pd.concat(results).sort_index()
    df_result.drop(columns=['chrom'], inplace=True)
    #keep only clusters with more than 1 junction
    df_result_AS = df_result[df_result['cluster'].duplicated(keep=False)]

    # remove clusters composed entirely of '-' junctions in either h_junction or m_junction
    # if a column does not exist it is ignored
    good_clusters = []
    for clu, grp in df_result.groupby('cluster'):
        bad = False
        for col in ('h_junction', 'm_junction'):
            if col in grp.columns and grp[col].astype(str).str.contains('-').all():
                bad = True
                break
        if not bad:
            good_clusters.append(clu)
    if good_clusters:
        df_result = df_result[df_result['cluster'].isin(good_clusters)]
        df_result_AS = df_result_AS[df_result_AS['cluster'].isin(good_clusters)]
    else:
        # no good clusters, return empty frames
        df_result = df_result.iloc[0:0]
        df_result_AS = df_result_AS.iloc[0:0]

    df_result.set_index(junctionCol, inplace=True, drop=True)
    df_result_AS.set_index(junctionCol, inplace=True, drop=True)
    return df_result, df_result_AS

def main():
    #if len(sys.argv) != 2:
    #    print("Usage: python get_input.py <value>")
    #    sys.exit(1)
    #version = sys.argv[1] #version of the results
    version = "HN6"
    #min_psi = 0.0 #remove junctions that are less than 5% from all the junctions in a clusters, in all samples
    min_jsr = 10 #remove junctions with less than 10 reads in at least minSamples samples from the same cell type group   
    minSamples = 2 #remove junctions with less than minSamples samples with minReads reads
    main_dir = '/gpfs0/tals/projects/Analysis/human_mouse_exons/'
    #human_mouse
    groups = ['GSE115736', 'GSE116177']
    species = ['h', 'm']
    #human_human
    #groups = ['GSE115736', 'GSE60424']
    #species = ['h', 'h']
    #mouse_mouse
    #groups = ['GSE116177', 'GSE180020']
    #species = ['m', 'm']
    ensembl_dir = main_dir + 'ensembl115/'
    input_file = ensembl_dir + "junctions_merge_" + groups[0] + "_" + groups[1] + "_" + version + ".txt"
    merged_df = pd.read_csv(input_file, sep="\t", index_col=0)
    #remove junctions with less than min_jsr reads in at least minSamples samples from the same cell type group
    cell_type_group_1 = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono', 'Neut']
    cell_type_group_2 = ['CD4T', 'Cd8T', 'BCell', 'NK', 'Mono', 'Neut']
    
    # Filter junctions: keep only those with at least min_jsr reads in at least minSamples samples from at least one cell type
    filtered_junctions = set()
    for i in range(len(cell_type_group_1)):
        df_filtered = merged_df[[col for col in merged_df.columns if cell_type_group_1[i] in col or cell_type_group_2[i] in col]]
        # Check each junction in df_filtered
        for idx, row in df_filtered.iterrows():
            count = (row >= min_jsr).sum()
            if count >= minSamples:
                filtered_junctions.add(idx)
    
    merged_df = merged_df.loc[list(filtered_junctions)]
    print(merged_df.columns)
    print("After filtering for min_jsr and minSamples:", len(merged_df), "junctions remain.")
    startNumCol = 7
    if species[0] == species[1]:
        junction_sum_final, junction_AS = assign_clusters_to_groups(merged_df, 'junction_id', startNumCol)
    else:
        junction_sum_final, junction_AS = assign_clusters_to_groups(merged_df, 'h_junction', startNumCol)
    
    # Sort by cluster, then by h_junction
    junction_AS = junction_AS.reset_index().sort_values(['cluster', 'h_junction']).set_index('h_junction')
    
    output_file = ensembl_dir + "junctions_cluster_" + groups[0] + "_" + groups[1] + "_" + version + ".txt"
    junction_sum_final.to_csv(output_file, sep="\t")
    output_file_AS = ensembl_dir + "junctions_cluster_AS_" + groups[0] + "_" + groups[1] + "_" + version + ".txt"
    junction_AS.to_csv(output_file_AS, sep="\t")
    print(junction_AS['cluster'].nunique(), 'unique clusters,', len(junction_AS), 'AS junctions', groups[0], groups[1])
    return
        
if __name__ == "__main__":
    main() 