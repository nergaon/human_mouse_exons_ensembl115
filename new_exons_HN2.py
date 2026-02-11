#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 15:04:37 2025

@author: nergaon
"""

import pandas as pd
import networkx as nx
import sys

# Function to find triplet indices within a chromosome group
def find_triplet_indices(group_df, junction_col):
    # Split into chr, start, end
    group_df[['chr', 'start', 'end']] = group_df[junction_col].str.split(':', expand=True)
    group_df['start'] = pd.to_numeric(group_df['start'], errors='coerce')
    group_df['end'] = pd.to_numeric(group_df['end'], errors='coerce')
    triplets = []
    for i, row in group_df.iterrows():
        # skip invalid coordinates
        if row["start"] < 0 or row["end"] < 0:
            continue
        start_matches = group_df[
            (group_df["chr"] == row["chr"]) &
            (group_df["start"] == row["start"]) &
            (group_df[junction_col] != row[junction_col]) &
            (group_df["start"] >= 0) &
            (group_df["end"] >= 0)
        ]
        end_matches = group_df[
            (group_df["chr"] == row["chr"]) &
            (group_df["end"] == row["end"]) &
            (group_df[junction_col] != row[junction_col]) &
            (group_df["start"] >= 0) &
            (group_df["end"] >= 0)
        ]
        for _, s_row in start_matches.iterrows():
            for _, e_row in end_matches.iterrows():
                # longest junction condition
                if (row["start"] <= s_row["start"] and row["start"] <= e_row["start"] and
                    row["end"] >= s_row["end"] and row["end"] >= e_row["end"]):
                    triplets.append([i, s_row.name, e_row.name])
    return triplets

def mark_triplet_types(df: pd.DataFrame, triplets: list, col_name):
    for triplet in triplets:
        if len(triplet) != 3:
            continue  # skip malformed triplets
        se_idx = triplet[0]
        ie_indices = triplet[1:]
        # mark SE (the first index in each triplet)
        df.loc[se_idx, col_name] = "SE"
        # mark IE (the remaining indices)
        df.loc[ie_indices, col_name] = "IE"
    return df

def assign_clusters_to_groups(df):    
    # Adding a new column 'Sequential_Index' with sequential numbers
    #df.reset_index(inplace=True, drop=False)               
    # Extract chromosome part
    df.loc[:, 'chrom'] = df['h_junction'].astype(str).str.split(':').str[0]
    # Sample columns (exclude metadata columns)
    #sample_cols = [col for col in df.columns if col not in ['h_junction', 'chrom']]
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
            _, start, end = row['h_junction'].split(':')
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
            # sample_totals = df_cluster[sample_cols].sum()
            # # Create a boolean DataFrame: True if a junction is ≥ min_psi of the sample total
            # high_expr_mask = (df_cluster[sample_cols] >= (min_psi * sample_totals))
            # # Keep rows (junctions) with at least one True in any sample
            # keep_mask = high_expr_mask.any(axis=1)
            # df_cluster_filtered = df_cluster[keep_mask].copy()
            # if df_cluster_filtered.empty:
            #     continue  # Skip if nothing left
    
            # Step 3: Rebuild graph for filtered cluster to check for disconnection
            G_filtered = nx.Graph()
            for idx, row in df_cluster.iterrows():
                _, start, end = row['h_junction'].split(':')
                G_filtered.add_edge(idx, f'start_{start}')
                G_filtered.add_edge(idx, f'end_{end}')
    
            subcomponents = list(nx.connected_components(G_filtered))
    
            for subcomp in subcomponents:
                sub_junctions = [node for node in subcomp if node in df_cluster.index]
                if not sub_junctions:
                    continue
                df_subcluster = df_cluster.loc[sub_junctions].copy()
                df_subcluster['cluster'] = f'clu_{cluster_counter}'
                cluster_counter += 1
                results.append(df_subcluster)
    
    # Combine all results
    df_result = pd.concat(results).sort_index()
    df_result.drop(columns=['chrom'], inplace=True)
    #df_result.set_index('h_junction', inplace=True, drop=True)
    return df_result

def main():
    if len(sys.argv) != 2:
        print("Usage: python get_input.py <value>")
        sys.exit(1)
    version = sys.argv[1] #version of the results
    main_dir = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/'
    group_1 = 'GSE115736'
    group_2 = 'GSE116177'
    clusters_input = main_dir + 'junctions_sum_' + group_1 + "_" + group_2 + "_" + version + ".txt"
    df = pd.read_csv(clusters_input, sep='\t', index_col=0)   
    #selecting rows where at least two numeric columns have values greater than or equal to 10
    numeric_df = df.select_dtypes(include='number')    
    # Create a boolean mask: True where value >= 10
    mask = numeric_df >= 10    
    # Count how many columns per row have value >= 10
    count = mask.sum(axis=1)    
    # Filter rows where count >= 2
    express_df = df[count >= 2]
    df_cluster = assign_clusters_to_groups(express_df)    
   
    #clusters with at least 2 rows
    filtered_df = df_cluster.groupby("cluster").filter(lambda x: len(x) >= 2)
    filtered_df = filtered_df.sort_values(by="cluster")
    #find skip exons
    filtered_df["h_type"] = ""
    filtered_df["m_type"] = ""
    #look for SE in each clsuter
    SE_df = filtered_df.groupby("cluster").filter(lambda x: len(x) >= 3)
    for cluster_name, group in SE_df.groupby('cluster'):
        print(cluster_name)
        triplet_indices_h = find_triplet_indices(group, 'h_junction')
        filtered_df = mark_triplet_types(filtered_df, triplet_indices_h, "h_type")
       
        triplet_indices_m = find_triplet_indices(group,'m_junction')
        filtered_df = mark_triplet_types(filtered_df, triplet_indices_m, "m_type")
        
    output = main_dir + 'junctions_with_type_' + version + '.txt'
    filtered_df.to_csv(output, sep="\t", index=False)
    print("all junctions", len(df))
    print("express junctions", len(express_df))
    return
          
if __name__ == "__main__":
      main() 