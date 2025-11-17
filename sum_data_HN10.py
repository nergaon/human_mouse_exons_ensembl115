#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 11:06:52 2024

@author: nergaon
"""
import pandas as pd
#import re
#from collections import Counter

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

def add_gene(input_file, gene, i, data_df, num_points):
    #read new gene data
    new_gene = pd.read_csv(input_file, sep="\t", index_col=0)
    new_gene['orthologous_points'] = num_points
    data_df = pd.concat([data_df, new_gene])
    return(data_df)

def get_unique_points(orthologs_df, identity_df, output_unique_junctions):
    #unique junctions are the junctions with the highest averatge id
    orthologs_df['average_id'] = 0.0
    for index, row in orthologs_df.iterrows():
        human_exon = orthologs_df['exon_h'][index]
        mouse_exon = orthologs_df['exon_m'][index]
        average_id = (identity_df[human_exon][mouse_exon] + identity_df[mouse_exon][human_exon]) / 2   
        #orthologs_df['average_id'][index] = average_id
        orthologs_df.loc[index, 'average_id'] = average_id
    orthologs_df['star_count'] = orthologs_df.apply(count_stars, axis=1)
    # Sort the DataFrame by 'average_id' (descending), 'human_exons' (ascending by number of '*'), 'mouse_exons' (ascending by number of '*')
    unique_junctions = orthologs_df.sort_values(by=['position_h','star_count','average_id'], ascending=[True, True, False])
    # Drop duplicates based on 'human_junction' and keep the first occurrence
    unique_junctions = unique_junctions.drop_duplicates(subset=['position_h'], keep='first')
    
    #do the same for mouse
    unique_junctions = unique_junctions.sort_values(by=['position_m','star_count','average_id'], ascending=[True, True, False])
    unique_junctions = unique_junctions.drop_duplicates(subset=['position_m'], keep='first')
    unique_junctions = unique_junctions.reset_index(drop=True)
    unique_junctions.to_csv(output_unique_junctions, sep="\t")
    unique_junctions.drop(columns=['average_id','star_count'], inplace=True)
    return(unique_junctions)

#count '*' in rank_h and rank_m
def count_stars(row):
    return (str(row['rank_h']).count('*') +
            str(row['rank_m']).count('*'))

def main():
    version = 'HN4'
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
    for i in range(0, len(mouse_ortologs)):
        if i%500 == 0:
            print(i)
        gene = mouse_ortologs['hsapiens_homolog_ensembl_gene'][i]
        if isinstance(gene, float):
            gene = "novel_gene"
        input_h = main_dir + genes_folder + gene + "_" + str(i) + '/human_statistics.txt'
        input_m = main_dir + genes_folder + gene + "_" + str(i) + '/mouse_statistics.txt'
        input_ortologs = main_dir + genes_folder + gene + "_" + str(i) + '/orthologs_unique_points.txt'
        input_identity = main_dir + genes_folder + gene + "_" + str(i) + '/alignment_identity.txt'
        try:
            #create one statistic df
           orthologs_df = pd.read_csv(input_ortologs, sep="\t", index_col=0)
           #one point in human or mouse is ortholog to more than 1 point
           if orthologs_df['position_h'].duplicated().any() or orthologs_df['position_m'].duplicated().any():
               identity_df = pd.read_csv(input_identity, sep="\t", index_col=0)
               output_unique_junctions = main_dir + genes_folder + gene + "_" + str(i) + '/orthologs_unique_points_B.txt'
               orthologs_df = get_unique_points(orthologs_df, identity_df, output_unique_junctions)
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