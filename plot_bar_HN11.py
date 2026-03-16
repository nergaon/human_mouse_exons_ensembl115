# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 11:11:13 2020

@author: hadas
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
'''
plot count and PSI bar from all the clusters in the JSR counts table
'''

def plot_bar_2(persent_cluster, value_cluster, fig_title, main_dir, main_col, legend_color_map):
    '''
    input: datafaram with the clusters values in percent
    output: bar plot of the introns
    '''
    output_1 = main_dir + "genes_figs/" + fig_title + '.svg' #png, pdf, jpeg, , tiff, raw, bmp, eps or svg
    output_2 = main_dir + "genes_figs/" + fig_title + '.png'
    #make sure the junctions in the value and % tables are in the same order
    persent_cluster = persent_cluster.sort_index()
    value_cluster = value_cluster.sort_index()
    #get introns position as the legends
    intron_legend = persent_cluster.index
    intron_legend  = intron_legend.tolist() 
    legend_color_map = create_color_map(intron_legend, legend_color_map)

    #change col names to be only the species
    persent_cluster.columns = [col.split('_')[-1] for col in persent_cluster.columns[0:]]
    value_cluster.columns = [col.split('_')[-1] for col in value_cluster.columns[0:]]
    #value_cluster.drop(columns=[main_col], inplace = True)
    #persent_cluster.drop(columns=[main_col], inplace = True)
    index = []
    for i in range(0, len(persent_cluster.columns)-1):      
        if persent_cluster.columns[i] != persent_cluster.columns[i+1]:
            index.append(i+1)
    #add col with 0 to seperate different samples
    zeroCol = []
    runNum = 0
    for i in index:
        col_name = str(i)
        i = i + runNum
        persent_cluster.insert(i, col_name, 0)
        value_cluster.insert(i, col_name, 0)
        runNum += 1   
        zeroCol.append(i)
    #the names in x axis
    labels = persent_cluster.columns
    labels = labels.tolist()
    unique_names = []
    for i in range(len(labels)):
        if i in zeroCol:
            labels[i] = ""
        #elif i%5 != 1: #second sample will have the name
        elif labels[i] in unique_names: #first sample will have the name
            labels[i] = ' '
        else:
            unique_names.append(labels[i])
    N = int(len(labels)) 
    r = np.zeros(shape=(N))  # the x locations for the groups. 
    for i in range(0, N-1):
        r[i+1] = r[i] + 0.25 
    barWidth = 0.2  # the width of the bars 
    plt.close('all')
    #my_title = gene
    introns_persent = persent_cluster.values
    introns_counts = value_cluster.values
    bars_persent = np.zeros(shape=(N))
    bars_counts = np.zeros(shape=(N))
    #fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4,4))
    # Create subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 4), gridspec_kw={'height_ratios': [1, 1]})

    for i in range(0, len(introns_persent)):
        color = legend_color_map[intron_legend[i]]
        ax1.bar(r, introns_persent[i], width = barWidth, bottom=bars_persent, edgecolor = "none", color=color)
        bars_persent = np.add(bars_persent, introns_persent[i])
        ax2.bar(r, introns_counts[i], width = barWidth, bottom=bars_counts, edgecolor = "none", color=color)
        bars_counts = np.add(bars_counts, introns_counts[i])
    
    ax1.set_xticks(r) 
    ax1.set_xticklabels(labels, fontsize=0)  #no cell names in the upper fig
    ax1.tick_params(axis="x", which='both', length=0) 
    ax1.set_ylim(0, 1)
    ax1.set_yticks([0, 1]) 
    ax1.tick_params(axis="y", labelsize=10) 
    ax1.set_ylabel('PSI')
    
    ax2.set_xticks(r) 
    ax2.tick_params(axis="x", which='both', length=0) 
    #ax2.set_xticklabels(labels, fontsize=10, rotation=30) 
    ax2.set_xticklabels(labels, fontsize=10, rotation=0)
    ax2.tick_params(axis="y", labelsize=10)
    #ax2.legend(intron_legend, fontsize=8, markerscale=0.15)
    ax2.set_ylabel('Counts')
    # Add the legend below both plots
    fig.legend(intron_legend, fontsize=8, markerscale=0.15, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=1)
    fig.suptitle(fig_title, fontsize=12)
    # Adjust layout to make space for the legend
    plt.subplots_adjust(hspace=0.4)
    #plt.show()
    plt.savefig(output_1, bbox_inches='tight')
    plt.savefig(output_2, bbox_inches='tight')
    plt.close('all')
    return (legend_color_map)

def create_color_map(intron_legend, legend_color_map):
    """
    Create or update a color map for intron legends.
    - Keeps existing colors.
    - Assigns new colors only to new introns.
    """
    #legend_color_map = {}
    cmap = plt.get_cmap("tab20")  # or any other colormap
    color_index = 0
    used_color = [] #in the same cluster, don't use the same colors
    for intron in intron_legend:
        if intron in legend_color_map:
           used_color.append(legend_color_map[intron]) 
    for intron in intron_legend:
        if intron not in legend_color_map:
            while cmap(color_index % cmap.N) in used_color:
               color_index += 1 
            # Find the next unused color
            legend_color_map[intron] = cmap(color_index % cmap.N)
            color_index += 1
    return legend_color_map

def main():  
    legend_color_map = {} #color dict
    if len(sys.argv) != 2:
        print("Usage: python get_input.py <value>")
        sys.exit(1)
    version = sys.argv[1] #version of the results
    #version = "HN6"
    #fibroblasts
    #cell_type = ['Fibroblast']
    #group_1 = 'GSE121052'
    #group_2 = 'GSE161648'
    #main_col = 'junction'
    #fibroblasts
    #group_1 = 'E-MTAB-5919_human'
    #group_2 = 'E-MTAB-5919_mouse'
    #human_mouse
    cell_type = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono', 'Neut']
    main_col = 'h_junction'
    group_1 = 'GSE115736'
    group_2 = 'GSE116177'
    #human_human
    # cell_type = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono', 'Neut']
    # main_col = 'junction'
    # group_1 = 'GSE115736'
    # group_2 = 'GSE60424'
    #mouse_mouse
    # cell_type = ['CD4T', 'Cd8T', 'BCell', 'NK', 'Mono']
    # group_1 = 'GSE116177'
    # group_2 = 'GSE180020'
    # main_col = 'mouse_junction'
    
    main_dir = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/leafcutter_' + group_1 + '_' + group_2 + '/'
    sum_table_input = main_dir + 'clusters_sum_table_' + version + '.txt'
    sum_df = pd.read_csv(sum_table_input, sep='\t',index_col=0)
    
    for one_cell in cell_type:
        print(one_cell)
        cell_table = sum_df[[main_col] + sum_df.filter(like=one_cell).columns.tolist()]
        cell_col_p = one_cell + '_p.adjust'
        cell_col_deltaPSI = one_cell + '_abs_deltapsi'
        #cell_col_cluster = one_cell + '_cluster'
        sucess_clusters = cell_table.loc[cell_table[cell_col_p].notnull(), 'cluster'].unique()
        percent_table_input = main_dir + one_cell + "/AS_clusters_psi_" + version + ".txt"
        percent_table = pd.read_csv(percent_table_input, sep='\t')
        # Split the 'Unnamed: 0' column into two columns
        percent_table[['h_junction', 'cluster']] = percent_table['Unnamed: 0'].str.extract(r'(.+):(clu_.*)')
        # Set the 'index' column as the index of the DataFrame
        percent_table.set_index('h_junction', inplace=True)
        # Drop the original 'Unnamed: 0' column
        percent_table.drop(columns=['Unnamed: 0'], inplace=True)
        value_table_input = main_dir + one_cell + "/AS_clusters_value_" + version + ".txt"
        value_table = pd.read_csv(value_table_input, sep=' ')
        # Split the 'Unnamed: 0' column into two parts
        value_table[[main_col, 'cluster']] = value_table['Unnamed: 0'].str.rsplit(':', n=1, expand=True)
        # Drop the original 'Unnamed: 0' column
        value_table.drop(columns=['Unnamed: 0'], inplace=True)
        value_table.set_index('h_junction', inplace=True)
        count_fig = 0
        for cluster in sucess_clusters: #select clusters to plot
            print(cluster)
            #print(cluster)
            p_value = cell_table[cell_table['cluster'] == cluster][cell_col_p].unique()
            formatted_p_value = "{:.3f}".format(p_value[0])
            deltapsi = cell_table[cell_table['cluster'] == cluster][cell_col_deltaPSI].abs().max()
            formatted_deltapsi = "{:.3f}".format(deltapsi)
            comment = 'p=' + str(formatted_p_value) + "_deltapsi=" + str(formatted_deltapsi)
            count_fig = count_fig + 1
            #select the cluter and its introns in the same order in the count and persent db
            persent_cluster = percent_table.loc[percent_table['cluster'] == cluster]
            persent_cluster = persent_cluster.sort_values(by=[main_col]) 
            gene_name = sum_df.loc[sum_df['cluster'] == cluster, 'symbol_h'].values[0]
            #persent_cluster.drop(columns=['ensembl','symbol','mouse_junction','cluster'], inplace = True)
            persent_cluster.drop(columns=['cluster'], inplace = True)
            value_cluster = value_table[value_table.index.isin(persent_cluster.index)]
            value_cluster = value_cluster.drop(columns=['cluster'])
            #plot fig
            fig_title = one_cell + "_" + gene_name + "_" + cluster + "_" + comment
            legend_color_map = plot_bar_2(persent_cluster, value_cluster, fig_title, main_dir, main_col, legend_color_map) 
            #fig_title = one_cell + "_" + gene_name
            #plot_bar_1(persent_cluster, fig_title, main_dir, main_col)
        print(count_fig, 'success clusters in', one_cell)
          
    return
    
if __name__ == "__main__":
    main()   
    
