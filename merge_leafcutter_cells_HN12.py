#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
matplotlib.use('Agg')  # Use a non-interactive backend
import sys

def plot_p(df, one_cell, main_dir, group_1, group_2):
    #output_fig = main_dir + 'leafcutter_0.2.9/' + one_cell + '_p_values.pdf'
    output_fig = main_dir + 'leafcutter_' + group_1 + '_' + group_2 + '/' + one_cell + '_p_values.pdf'
    # Drop rows with NaN values in any column
    df = df.dropna()
    df_sorted = df.sort_values(by='p.adjust', ascending=False)
    # for each cell type, define corresponding colors
    # Create the plot
    with PdfPages(output_fig) as pdf:
        plt.figure(figsize=(14, 8)) 
        plt.scatter(df_sorted['extracted_cluster'], -np.log10(df_sorted['p.adjust']), label=one_cell)
        # Add a horizontal line at -log10(0.05)
        plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='-log10(0.05)')
        # Adding labels and title
        plt.xlabel('cluster', fontsize=10)
        plt.ylabel('-log10(pAdjust)', fontsize=12)
        plt.title(one_cell,fontsize=14)
        plt.legend(fontsize=12)
        #plt.grid(True)
        
        # Show the plot
        pdf.savefig()
        plt.show()
    return

def plot_p_vs_deltapsi(ax, df, one_cell):
    ax.scatter(df['deltapsi'], -np.log10(df['p.adjust']))
    # Add labels
    ax.set_xlabel('abs(max DeltaPsi)')
    ax.set_ylabel('-log10(p.adjust)')
    #add a horizontal line at -log10(0.05)
    ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='-log10(0.05)')
    ax.set_title(one_cell,fontsize=14)
    # Add the names of the genes with the 0.5% lowest p.adjust
    threshold = np.percentile(df['p.adjust'], 0.5)
    for i in range(len(df)):
        if df['p.adjust'][i] <= threshold:
            try:
                ax.text(df['deltapsi'][i], -np.log10(df['p.adjust'][i]), df['genes'][i], fontsize=9)
            except:
                ax.text(df['deltapsi'][i], -np.log10(df['p.adjust'][i]), df['cluster'][i], fontsize=9)
    # Add the names of the genes with the 0.5% higest deltapsi
    threshold = np.percentile(df['deltapsi'], 99.5)
    for i in range(len(df)):
        if df['deltapsi'][i] >= threshold:
            try:
                ax.text(df['deltapsi'][i], -np.log10(df['p.adjust'][i]), df['genes'][i], fontsize=9)
            except:
                ax.text(df['deltapsi'][i], -np.log10(df['p.adjust'][i]), df['cluster'][i], fontsize=9)
    return
   
def main():
    #if len(sys.argv) != 2:
    #    print("Usage: python get_input.py <value>")
    #    sys.exit(1)
    #version = sys.argv[1] #version of the results
    version = "HN6"
    deltapsi = 0.2
    #human_mouse
    #cell_type = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono', 'Neut']
    #group_1 = 'GSE115736'
    #group_2 = 'GSE116177'
    #col_name = 'human_junction'
    #cell_type = ['Fibroblast']
    #group_1 = 'GSE121052'
    #group_2 = 'GSE161648'
    #merge_col = 'human_junction'
    #group_1 = 'E-MTAB-5919_human'
    #group_2 = 'E-MTAB-5919_mouse'
    #human_human
    cell_type = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono', 'Neut']
    group_1 = 'GSE115736'
    group_2 = 'GSE60424'
    col_name = 'junction_id' 
    #mouse_mouse. in this comparison there is no gene col. this is why there is try and except
    # cell_type = ['CD4T', 'Cd8T', 'BCell', 'NK', 'Mono']
    # group_1 = 'GSE116177'
    # group_2 = 'GSE180020'
    # merge_col = 'mouse_junction'
    #col_name = 'junction_id' 
    main_dir = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/'
    clusters_input = main_dir + 'junctions_cluster_' + group_1 + "_" + group_2 + "_" + version + '.txt'
    clusters_df = pd.read_csv(clusters_input, sep='\t')
    # Remove rows with clusters that have only 1 junction
    cluster_counts = clusters_df['cluster'].value_counts()
    clusters_to_keep = cluster_counts[cluster_counts > 1].index
    clusters_df = clusters_df[clusters_df['cluster'].isin(clusters_to_keep)]
    
    # Count the number of junctions and clusters
    num_junctions = len(clusters_df)
    num_clusters = clusters_df['cluster'].nunique()
    print(f"Number of junctions after filtering: {num_junctions}")
    print(f"Number of clusters after filtering: {num_clusters}")

    # Select only columns with dtype 'object' or 'string'
    clusters_df_data  = clusters_df.select_dtypes(include=['object', 'string']).copy()
    #to the sum data table, add leafcutter results
    sum_df_input = main_dir + 'leafcutter_' + group_1 + '_' + group_2 + '/sum_table_' + version + '.txt'
    sum_df = pd.read_csv(sum_df_input, sep='\t', index_col=0)
    output_fig = main_dir + '/leafcutter_' + group_1 + '_' + group_2 + '/' + group_1 + "_" + group_2 + '_p_deltapis.pdf'
    output_fig_jpeg = main_dir + '/leafcutter_' + group_1 + '_' + group_2 + '/' + group_1 + "_" + group_2 + '_p_deltapis.jpeg'
    
    # Add cluster information once before the loop (since clusters are the same for all cells)
    first_cell = cell_type[0]
    intron_input = main_dir + 'leafcutter_' + group_1 + '_' + group_2 + '/' + first_cell + '/leafcutter_ds_effect_sizes.txt'
    intron_sig = pd.read_csv(intron_input, sep='\t')
    intron_sig['junction'] = intron_sig['intron'].apply(lambda x: ':'.join(x.split(':')[:3]))
    intron_sig['cluster'] = intron_sig['intron'].apply(lambda x: x.split(':')[-1])
    #intron_sig_part = intron_sig[['junction','cluster']].drop_duplicates()
    intron_sig_part = intron_sig[['junction']].drop_duplicates()
    clusters_df_data = pd.merge(clusters_df_data, intron_sig_part, left_on=col_name, right_on='junction', how='left')
    clusters_df_data.drop(columns=['junction'], inplace=True)
    
    with PdfPages(output_fig) as pdf:
        # Create a single figure with 6 subplots
        fig, axs = plt.subplots(3, 2, figsize=(8, 12))  # 3 rows x 2 columns
        i = 0
        for one_cell in cell_type:
            print(one_cell)
            cluster_sig_input = main_dir + 'leafcutter_' + group_1 + '_' + group_2 + '/' + one_cell + '/leafcutter_ds_cluster_significance.txt'
            cluster_sig = pd.read_csv(cluster_sig_input, sep='\t')
            # Extract the part after the colon in the 'cluster' column
            cluster_sig['extracted_cluster'] = cluster_sig['cluster'].str.split(':').str[1]
            genes_df = cluster_sig[['extracted_cluster', 'p.adjust','genes']]
            cluster_sig_part = cluster_sig[['extracted_cluster', 'p.adjust']].drop_duplicates(subset=['extracted_cluster'], keep='first')
            # Count the number of rows with non-null values in the 'p.adjust' column
            non_null_count = cluster_sig_part['p.adjust'].notna().sum()
            # Count how many of those values are <= 0.05
            below_threshold_count = cluster_sig_part['p.adjust'].le(0.05).sum()
            print(f"Number success clusters: {non_null_count}")
            print(f"Number of clusters with 'p.adjust' <= 0.05: {below_threshold_count}")
            sum_df.loc["Leafcutter success clusters", one_cell] = non_null_count
            sum_df.loc["Leafcutter sig clusters (p<0.05)", one_cell] = below_threshold_count
            plot_p(cluster_sig_part, one_cell, main_dir, group_1, group_2)
            intron_input = main_dir + 'leafcutter_' + group_1 + '_' + group_2 + '/' + one_cell + '/leafcutter_ds_effect_sizes.txt'
            intron_sig = pd.read_csv(intron_input, sep='\t')
            #get the junction without the clu
            intron_sig['junction'] = intron_sig['intron'].apply(lambda x: ':'.join(x.split(':')[:3]))
            intron_sig['cluster'] = intron_sig['intron'].apply(lambda x: x.split(':')[-1])
            # Add a column for absolute deltapsi
            intron_sig['abs_deltapsi'] = intron_sig['deltapsi'].abs()
            # Get abs_deltapsi for each junction
            intron_sig_part = intron_sig[['junction','abs_deltapsi']].drop_duplicates(subset=['junction'], keep='first')
            clusters_df_data = pd.merge(clusters_df_data, intron_sig_part, left_on=col_name, right_on='junction', how='left')
            clusters_df_data.drop(columns=['junction'], inplace=True)
            # Merge p.adjust values using the cluster column - rename extracted_cluster to cluster for merging
            cluster_sig_part_renamed = cluster_sig_part.rename(columns={'extracted_cluster': 'cluster'})
            clusters_df_data = pd.merge(clusters_df_data, cluster_sig_part_renamed, left_on='cluster', right_on='cluster', how='left', suffixes=('', f'_{one_cell}'))
            clusters_df_data.drop_duplicates(inplace=True)
            #new_name = one_cell + "_deltapsi"
            # Rename columns (rename abs_deltapsi with cell type, keep cluster for next iteration)
            clusters_df_data.rename(columns={'abs_deltapsi': f'{one_cell}_abs_deltapsi', 'p.adjust': f'{one_cell}_p.adjust'}, inplace=True)
            #plot
            intron_sig_part = intron_sig[['cluster', 'deltapsi']]
            intron_sig_part.loc[:, 'deltapsi'] = intron_sig_part['deltapsi'].abs()
            # Group by 'cluster' and keep the row with the highest 'deltapsi' for each cluster
            intron_sig_part = intron_sig_part.loc[intron_sig_part.groupby('cluster')['deltapsi'].idxmax()]
            intron_sig_part = pd.merge(intron_sig_part, genes_df, left_on = 'cluster', right_on = 'extracted_cluster')
            ax = axs[i // 2, i % 2]  # Determine subplot position
            #lot_p_vs_deltapsi(ax,intron_sig_part, one_cell)
            i = i + 1
            # Count the number of rows with 'p.adjust <= 0.05' and 'deltapsi >= 0.1'
            count = intron_sig_part[(intron_sig_part['p.adjust'] <= 0.05) & (intron_sig_part['deltapsi'] >= deltapsi)].shape[0]
            print ('sig deltapsi', deltapsi)
            print(f"Number of clusters with p.adjust <= 0.05 and deltapsi >= sig deltapsi: {count}")
            sum_df.loc["Leafcutter sig and deltapsi>=0.2", one_cell] = count            
            count = intron_sig_part[(intron_sig_part['p.adjust'] <= 0.05) & (intron_sig_part['deltapsi'] >= 0.1)].shape[0]
            sum_df.loc["Leafcutter sig and deltapsi>=0.1", one_cell] = count
            count = intron_sig_part[(intron_sig_part['p.adjust'] <= 0.05) & (intron_sig_part['deltapsi'] >= 0.05)].shape[0]
            sum_df.loc["Leafcutter sig and deltapsi>=0.05", one_cell] = count            
            #add expression data
            input_value = main_dir + 'leafcutter_' + group_1 + '_' + group_2 + '/' + one_cell + '/AS_clusters_value_' + version + '.txt'
            AS_value = pd.read_csv(input_value, sep=' ')
            # Remove the 'clu_*' part from the first column
            AS_value['junction'] = AS_value['Unnamed: 0'].apply(lambda x: ':'.join(x.split(':')[:3]))
            # Automatically identify human and mouse columns based on their names
            group_1_cols = [col for col in AS_value.columns if group_1 in col]
            group_2_cols = [col for col in AS_value.columns if group_2 in col]
            # Calculate the average for human and mouse columns
            col_name_1 = group_1 + '_avg_' + one_cell
            AS_value[col_name_1] = AS_value[group_1_cols].mean(axis=1)
            col_name_2 = group_2 + '_avg_' + one_cell
            AS_value[col_name_2] = AS_value[group_2_cols].mean(axis=1)
            AS_value_part = AS_value[['junction', col_name_1, col_name_2]]
            clusters_df_data = pd.merge(clusters_df_data, AS_value_part, left_on = col_name, right_on = 'junction', how='left')
            clusters_df_data.drop_duplicates(inplace=True)
            clusters_df_data.drop(columns=['junction'], inplace=True)
        # Adjust layout
        fig.tight_layout()
        # Save to PDF
        pdf.savefig(fig)
        fig.savefig(output_fig_jpeg, format="jpeg", dpi=300)
        plt.close(fig)    
    # Apply the function to each row and create a new column 'clu'
    #clusters_df['clu'] = clusters_df.apply(lambda row: extract_unique_clusters(row, clusters_df), axis=1)
    # Remove columns with 'cluster' in their names
    #clusters_df = clusters_df.drop(columns=clusters_df.filter(like='cluster').columns)
    #clusters_df['cluster'] = clusters_df.apply(lambda row: get_unique_clusters(row, clusters_df), axis=1)
    output = main_dir + 'leafcutter_' + group_1 + '_' + group_2 + '/clusters_sum_table_' + version + '.txt'
    clusters_df_data.to_csv(output, sep="\t")
    output_file = main_dir + 'leafcutter_' + group_1 + '_' + group_2 + '/sum_table_' + version + '.txt'
    sum_df.to_csv(output_file, sep="\t")
    return
          
if __name__ == "__main__":
      main() 