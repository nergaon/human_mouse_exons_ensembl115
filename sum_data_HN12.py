#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 11:06:52 2024

@author: nergaon
"""
import pandas as pd
import sys

def mark_transcripts(df, tra_num_input):
    #mark for each gene if it was done as should be or not
    df['Enough_exons'] = "yes"
    with open(tra_num_input, 'r') as f:
       tra_lines = f.readlines()
    for line in tra_lines:
        ensembl_m = line.strip().split(" ")[3]
        df.loc[df['ensembl_gene_id'] == ensembl_m, 'Enough_exons'] = 'no'
    print(len(tra_lines), 'genes with less than 2 orthologs exons')
    return(df)

# Group by 'position_h or _m' and select the best row per group
def select_best(df):
    # 1️⃣ Compute ratio-based score
    df['Ratio_based_score'] = df.apply(lambda r: r['Average %id'] * min(r['%id_human'], r['%id_mouse']) / max(r['%id_human'], r['%id_mouse']),axis=1)
    # 2️⃣ Count number of '*' in rank_h + rank_m
    df['star_count'] = df['rank_h'].astype(str).add(df['rank_m'].astype(str)).str.count(r'\*')
    # 3️⃣ Sort by position_h, score (descending), and star_count (ascending)
    df_sorted = df.sort_values(['position_h', 'Ratio_based_score', 'star_count'],ascending=[True, False, True])
    # 4️⃣ Pick best row per position_h
    best_rows = df_sorted.groupby('position_h').first().reset_index()
    best_rows = best_rows.sort_values(['position_m', 'Ratio_based_score', 'star_count'],ascending=[True, False, True])
    best_rows = best_rows.groupby('position_m').first().reset_index()
    return best_rows

def add_gene(input_file, gene, i, data_df, num_points):
    #read new gene data
    new_gene = pd.read_csv(input_file, sep="\t", index_col=0)
    new_gene['orthologous_points'] = num_points
    data_df = pd.concat([data_df, new_gene])
    return(data_df)

def main():
    if len(sys.argv) != 2:
        print("Usage: python get_input.py <value>")
        sys.exit(1)
    version = sys.argv[1] #version of the results
    #version = 'HN4'
    genes_folder = '/genes_' + version + '/'
    main_dir = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/'
    orthologs_input = main_dir + 'all/good_ortologs_df.txt'
    mouse_ortologs = pd.read_csv(orthologs_input, delimiter='\t', index_col=0)
    tra_num_input = main_dir + 'genesWithOneOrthologExon_' + version + '.txt'
    mouse_ortologs = mark_transcripts(mouse_ortologs, tra_num_input)    
    
    data_df_h = pd.DataFrame()
    data_df_m = pd.DataFrame()
    unique_points = pd.DataFrame()
    count_genes = 0 #genes with cds transcripts in human and mouse
    mouse_ortologs['Remarks'] = 'no folder'
    for i in mouse_ortologs.index:
    #for i in range(0, len(mouse_ortologs)):
        if i%500 == 0:
            print(i)
        gene = mouse_ortologs['hsapiens_homolog_ensembl_gene'][i]
        #print(gene)
        if isinstance(gene, float):
            gene = "novel_gene"
        input_h = main_dir + genes_folder + gene + "_" + str(i) + '/human_statistics.txt'
        input_m = main_dir + genes_folder + gene + "_" + str(i) + '/mouse_statistics.txt'
        #input_ortologs = main_dir + genes_folder + gene + "_" + str(i) + '/orthologs_unique_points.txt'
        input_ortologs = main_dir + genes_folder + gene + "_" + str(i) + '/orthologs_points.txt'
        try:
            #create one statistic df
           orthologs_df = pd.read_csv(input_ortologs, sep="\t", index_col=0)
           orthologs_df = select_best(orthologs_df)
           output_unique_junctions = main_dir + genes_folder + gene + "_" + str(i) + '/orthologs_unique_points_B.txt'
           orthologs_df.to_csv(output_unique_junctions, sep="\t")
           unique_points = pd.concat([unique_points, orthologs_df]) #data from all the genes
           mouse_ortologs.loc[mouse_ortologs['hsapiens_homolog_ensembl_gene'] == gene, 'Remarks'] = 'OK'  
           data_df_h = add_gene(input_h, gene, i, data_df_h, len(orthologs_df)) 
           data_df_m = add_gene(input_m, gene, i, data_df_m, len(orthologs_df)) 
           count_genes +=1         
        except:
            print('no files for gene', gene, i)
    #sum all the col
    print(count_genes, 'genes with more than 1 ortholog exons')
    data_df_h.loc['Total'] = data_df_h.iloc[:,1:].sum()
    data_df_m.loc['Total'] = data_df_m.iloc[:,1:].sum()
    output_h = main_dir + 'data_df_h_' + version + '.txt'
    data_df_h.to_csv(output_h, sep="\t")
    output_m = main_dir + 'data_df_m_' + version + '.txt'
    data_df_m.to_csv(output_m, sep="\t")
    output_exons = main_dir + 'unique_points_' + version + '.txt'
    unique_points = unique_points.reset_index(drop=True)
    unique_points.to_csv(output_exons, sep="\t")
    output_orthologs = main_dir + 'mouse_ortologs_remarks_' + version + '.txt'
    mouse_ortologs.to_csv(output_orthologs, sep="\t")
    return
          
if __name__ == "__main__":
      main() 