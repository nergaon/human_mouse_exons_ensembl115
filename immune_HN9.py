#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import networkx as nx
import sys
from general_def import assign_clusters_to_groups

def get_table (df, cell_type, species, main_dir):
    #from the leafcutter expression table, get the col of the cell i want
    
    # Select columns that start with cell_type
    #cell_df = df[[col for col in df.columns if col.startswith(cell_type)]]   
    # Selecting columns that include cell type
    cell_df = df.filter(regex=re.compile(cell_type, re.IGNORECASE))
    # Filter rows where all values are zeros
    cell_df = cell_df.loc[(cell_df != 0).any(axis=1)]
    # Remove '.sort.bam' from column names
    species = "_" + species
    cell_df.columns = cell_df.columns.str.replace('.sort.bam', species)
    # Apply the function to modify the index values
    #cell_df.index = cell_df.index.map(modify_index)
    return(cell_df)

def plot_histogram(df, cell, main_dir):
     # Flatten the DataFrame to get a single series of values
     numeric_values = df.values.flatten()   
     # Filter the values to include only those <= 100
     filtered_values = numeric_values[(numeric_values > 0) & (numeric_values <= 100)]
     output_fig = main_dir + 'JSR_distribution/' + cell + '_JSR_distribution.pdf'
     with PdfPages(output_fig) as pdf:
         plt.figure(figsize=(20, 8))
     # Create a histogram
         plt.hist(filtered_values, bins=100, edgecolor='black')
         plt.ylim(0, 35000)  # Set the x-axis limit
         plt.xticks(range(0, 101, 10))  # Set x-axis ticks at intervals of 10
         title_name = 'JSR distribution in ' + " " + cell + " cells"
         plt.title(title_name)
         plt.xlabel('JSR')
         plt.ylabel('Frequency')
         pdf.savefig()
         plt.show()
     return

def create_group_file(cell_dir, cell_df, group_1, group_2):
    """
    Creates a file with column names and their respective group identifiers.
    Parameters:
    df (pd.DataFrame): The input DataFrame.
    cell_dir (str): The path to the output file.
    The function writes a tab-separated file with the column names in the first column
    column and group_1 or group_2 in the second column based on the 
    presence of these identifiers in the column names.
    """
    # Determine the group for each column
    groups = [(col, group_1 if group_1 in col else group_2) for col in cell_df.columns]
    
    # Create a DataFrame for groups
    groups_df = pd.DataFrame(groups, columns=['Column_Name', 'Group'])
    
    # Save the DataFrame to a text file
    output_file = cell_dir + '/groups_file.txt'
    groups_df.to_csv(output_file, sep='\t', index=False, header=False)
    return
   
def main():   #clu_7543, ptprc, need to be split to 2 in CD8T
    #if len(sys.argv) != 2:
    #    print("Usage: python get_input.py <value>")
    #    sys.exit(1)
    #version = sys.argv[1] #version of the results
    version = "HN6"
    min_jsr = 0 #remove junctions with less than 10 reads in at least 2 samples   
    minSamples = 0 #remove junctions with less than minSamples samples with minReads reads
    #min_psi = 0.0 #remove junctions that are less than min_psi (for example, 5%) from all the junctions in a clusters, in all samples
    main_dir = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/'
    # Initialize an empty Index
    #junction_index = pd.Index([])
    #human_mouse
    #groups = ['GSE115736', 'GSE116177']
    #cell_type_group_1 = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono', 'Neut']
    #cell_type_group_2 = ['CD4T', 'Cd8T', 'BCell', 'NK', 'Mono', 'Neut']
    #human_human
    #groups = ['GSE115736', 'GSE60424']
    #cell_type_group_1 = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono', 'Neut']
    #cell_type_group_2 = ['CD4', 'CD8', 'Bcells', 'NK', 'Monocytes', 'Neutrophils']
    #mouse mouse
    groups = ['GSE116177', 'GSE180020']
    cell_type_group_1 = ['CD4T', 'Cd8T', 'BCell', 'NK', 'Mono']
    cell_type_group_2 = ['T_4', 'T_8', 'B_fo', 'NK', 'Mo']
    #human_mouse fibroblasts
    #group_1 = 'GSE121052'
    #group_2 = 'GSE161648'
    #human_mouse fibroblasts
    #group_1 = 'E-MTAB-5919_human'
    #group_2 = 'E-MTAB-5919_mouse'
    #cell_type_group_1 = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono', 'Neut']
    #cell_type_group_2 = ['CD4T', 'Cd8T', 'BCell', 'NK', 'Mono', 'Neut']
    #cell_type_group_1 = ['Fibroblast']
    #cell_type_group_2 = ['Fibroblast']
    clusters_input = main_dir + "junctions_cluster_AS_" + groups[0] + "_" + groups[1] + "_" + version + ".txt"
    clusters_df = pd.read_csv(clusters_input, sep='\t', index_col=0)
    columns = ["All"] + cell_type_group_1
    index = ["Orthologues Clusters", "Orthologues Junctions", "Express Clusters", "Express Junctions","AS Clusters", "AS Junctions"]  
    sum_df = pd.DataFrame(index=index, columns=columns)
    unique_cluster_count = clusters_df['cluster'].nunique()
    print("all junctions", len(clusters_df))
    sum_df.loc["Orthologues Clusters", 'All'] = unique_cluster_count
    sum_df.loc["Orthologues Junctions", 'All'] = len(clusters_df)    
    #clusters_index_df = pd.DataFrame(clusters_df.index, columns=['h_junction'])
    #each cell type seperatly
    for i in range(len(cell_type_group_1)):
        print(cell_type_group_1[i])
        df_filtered = clusters_df[[col for col in clusters_df.columns if cell_type_group_1[i] in col or cell_type_group_2[i] in col]]
        plot_histogram(df_filtered, cell_type_group_1[i], main_dir)
        # Count the number of columns with values >= 10 for each row
        junctions_filter = df_filtered.iloc[:,:].apply(lambda x: (x >= min_jsr).sum(), axis=1)
        # Filter to keep only rows that have at least 2 such columns
        junctions_df = df_filtered[junctions_filter >= minSamples]
        junctions_df['cluster'] = clusters_df.loc[junctions_df.index, 'cluster']
        print("express junctions:", len(junctions_df))
        unique_cluster_count = junctions_df['cluster'].nunique()
        sum_df.loc["Express Clusters", cell_type_group_1[i]] = unique_cluster_count
        sum_df.loc["Express Junctions", cell_type_group_1[i]] = len(junctions_df)
        #remove clusters with only one junction
        AS_events_value = junctions_df.groupby('cluster').filter(lambda x: len(x) > 1)
        #calculate percentages
        AS_events_psi = pd.DataFrame()
        AS_events_value_final = pd.DataFrame()
        cluster_list = list(AS_events_value.cluster.unique())  
        first = 1
        for cluster in cluster_list:
            #cluster = "clu_7543" #for testing
            one_cluster = AS_events_value.loc[AS_events_value['cluster'] == cluster]
            # Remove junctions that don't share any start or end points with other junctions
            if len(one_cluster) > 1:
                junction_indices = one_cluster.index.tolist()
                # Extract start and end points from junction indices (format: chr:start:end)
                junction_points = {}
                for junction in junction_indices:
                    parts = junction.split(':')
                    start = parts[1]
                    end = parts[2]
                    junction_points[junction] = (start, end)
                # Keep only junctions that share at least one point with another junction
                junctions_to_keep = []
                for junction in junction_indices:
                    start, end = junction_points[junction]
                    has_shared_point = False
                    for other_junction in junction_indices:
                        if junction != other_junction:
                            other_start, other_end = junction_points[other_junction]
                            if start == other_start or start == other_end or end == other_start or end == other_end:
                                has_shared_point = True
                                break
                    if has_shared_point:
                        junctions_to_keep.append(junction)               
                one_cluster = one_cluster.loc[junctions_to_keep]
                if first:
                    AS_events_value_final = one_cluster.copy()
                else:
                    AS_events_value_final = pd.concat([AS_events_value_final, one_cluster])
            else:
                print("Cluster", cluster, "has only one junction and will be removed.")
            # Calculate percentages for each column separately
            for col in one_cluster.columns[:-1]: 
                #convert to perecent only if there are more than 0 reads in a col
                col_sum = one_cluster[col].sum()
                if col_sum > 0:
                    one_cluster[col] = one_cluster[col].astype(float)
                    one_cluster.loc[:, col] = one_cluster[col] / one_cluster[col].sum()
                else:
                    one_cluster.loc[:, col] = 0
            if first:
                AS_events_psi = one_cluster.copy()
                first = 0
            else:
                frames = [AS_events_psi, one_cluster]
                AS_events_psi = pd.concat(frames)
        unique_cluster_count = AS_events_value_final['cluster'].nunique()
        print("Number of human mouse leafcutter AS clusters", unique_cluster_count)
        print("Number of human mouse leafcutter AS junctions", len(AS_events_value_final))
        sum_df.loc["AS Clusters", cell_type_group_1[i]] = unique_cluster_count
        sum_df.loc["AS Junctions", cell_type_group_1[i]] = len(AS_events_value_final)
        #junction_index = junction_index.union(AS_events_value_final.index)

        AS_events_value_final.index = AS_events_value_final.index + ":" + AS_events_value_final['cluster'].astype(str)
        AS_events_value_final = AS_events_value_final.drop(columns=['cluster'])
        AS_events_psi.index = AS_events_psi.index + ":" + AS_events_psi['cluster'].astype(str)
        AS_events_psi = AS_events_psi.drop(columns=['cluster'])
        
        #convert to leafcutter table
        cell_df = AS_events_value_final.copy()
        # Remove the index name
        cell_df.index.name = None
        # Replace '.' with '_' and '#' with ''
        cell_df.columns = cell_df.columns.str.replace('.', '_').str.replace('#', '')
        cell_dir = main_dir + 'leafcutter_' + groups[0] + '_' + groups[1] + '/' + cell_type_group_1[i]
        df_output =  cell_dir  + '/AS_clusters_psi_' + version + '.txt' 
        AS_events_psi.columns = AS_events_psi.columns.str.replace('.', '_').str.replace('#', '')
        AS_events_psi.to_csv(df_output, sep="\t")
        df_output =  cell_dir + '/AS_clusters_value_' + version + '.txt' 
        cell_df.to_csv(df_output, sep=" ")
        #duplicates = cell_df[cell_df.duplicated(keep=False)]
        create_group_file(cell_dir, cell_df, groups[0], groups[1])
    
    output_file = main_dir + 'leafcutter_' + groups[0] + '_' + groups[1] + '/sum_table_' + version + '.txt'
    sum_df.to_csv(output_file, sep="\t")
    return
          
if __name__ == "__main__":
      main() 