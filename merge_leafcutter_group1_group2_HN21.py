#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import json
import sys

def plot_histogram(df, species, cell, ensembl_dir):
     # Flatten the DataFrame to get a single series of values
     numeric_values = df.values.flatten()   
     # Filter the values to include only those <= 100
     filtered_values = numeric_values[(numeric_values > 0) & (numeric_values <= 100)]
     output_fig = ensembl_dir + 'JSR_distribution/' + cell + "_" + species + '_JSR_distribution.pdf'
     with PdfPages(output_fig) as pdf:
         plt.figure(figsize=(20, 8))
     # Create a histogram
         plt.hist(filtered_values, bins=100, edgecolor='black')
         plt.ylim(0, 35000)  # Set the x-axis limit
         plt.xticks(range(0, 101, 10))  # Set x-axis ticks at intervals of 10
         title_name = 'JSR distribution in ' + species + " " + cell + " cells"
         plt.title(title_name)
         plt.xlabel('JSR')
         plt.ylabel('Frequency')
         pdf.savefig()
         #plt.show()
     return

def get_rows(points_df, df, col_name):
    print('get_rows')
    species_dict = {"h": "m", "m": "h"}
    species = col_name[-1]
    second_col = 'position_' + species_dict[species]
    junction_col = species_dict[species] + "_junction"
    # Preprocess points_df into a dictionary for fast lookup
    points_dict = dict(zip(points_df[col_name], points_df[second_col]))
    # Prepare storage
    rows = []
    junction_dict = {}
    #run_index=0
    for i in df.index:
        # if run_index%1000 == 0:
        #     print(run_index)
        # run_index+=1
        chr_, start, end = i.split(":")
        point_a = f"{chr_}:{start}"
        point_b = f"{chr_}:{end}"
        value_a = points_dict.get(point_a, "")
        value_b = points_dict.get(point_b, "")
        # Skip rows where both values are missing
        if not value_a and not value_b:
            continue
        chrName = ""
        if value_a:
            chrName_a, genomic_position_a = value_a.split(":")
            junction_dict[point_a] = genomic_position_a
            chrName = chrName_a
        else:
            junction_dict[point_a] = -int(start)
        if value_b:
            chrName_b, genomic_position_b = value_b.split(":")
            junction_dict[point_b] = genomic_position_b
            chrName = chrName_b
        else:
            junction_dict[point_b] = -int(end)
        # Skip if both values exist but are from different chromosomes
        if value_a and value_b and chrName_a != chrName_b:
            continue
        # Create row and add metadata
        row = df.loc[[i]].copy()
        if int(junction_dict[point_a])>int(junction_dict[point_b]) and int(junction_dict[point_a])>0 and int(junction_dict[point_b])>0: #cd47 - different strand in human and mouse 
            row[junction_col] = f"{chrName}:{junction_dict[point_b]}:{junction_dict[point_a]}"
            junction_dict[point_a] = genomic_position_b
            junction_dict[point_b] = genomic_position_a
        else:
            row[junction_col] = f"{chrName}:{junction_dict[point_a]}:{junction_dict[point_b]}"
        rows.append(row)
    # Combine all rows into a single DataFrame
    orthologs_group = pd.concat(rows)
    # Modify index format
    orthologs_group.index = orthologs_group.index.map(modify_index)
    return orthologs_group, junction_dict

# Function to apply assign_clusters to each chromosome group with a unique global cluster count
def assign_clusters_to_groups(df):
    #remove :- from rows that have it in 'human_jucntion' col
    #df['human_junction'] = df['human_junction'].str.replace(':-', '', regex=False)
    #df['human_junction'] = df['human_junction'].str.replace(':+', '', regex=False)
    # Apply the function to create new columns with chromosome, start, and end positions
    df[['chr', 'start_pos', 'end_pos']] = df['h_junction'].apply(lambda x : pd.Series(extract_chr_start_end(x))) 
    # Initialize cluster index
    cluster_index = 0
    # Add 'cluster' column with initial value 0
    df['cluster'] = 0
    for name, group in df.groupby('chr'):
        print(name)     
        # Iterate over each row in the DataFrame
        while len(group) > 1:
            index = group.index[0]
            row = group.loc[index]
            # Get the start and end values of the current row
            start_list = [row['start_pos']]
            end_list = [row['end_pos']]
            # Increment cluster index and assign it to the current row
            cluster_index += 1
            df.at[index, 'cluster'] = cluster_index
            group.at[index, 'cluster'] = cluster_index
            # Initialize a flag to check if new rows are added to the cluster
            new_rows_added = True                
            while new_rows_added:
                new_rows_added = False                   
                # Iterate over each row again to find rows with shared start or end
                for idx, r in group.iterrows():
                    if r['cluster'] == 0 and (r['start_pos'] in start_list or r['end_pos'] in end_list):
                        # Assign the same cluster index to the row
                        df.at[idx, 'cluster'] = cluster_index
                        group.at[idx, 'cluster'] = cluster_index
                        # Add the start and end values to the lists
                        start_list.append(r['start_pos'])
                        end_list.append(r['end_pos'])                   
                        # Set flag to True indicating new rows are added
                        new_rows_added = True
            #remove the rows that were assign to clusters
            group = group[group['cluster'] == 0]
        if len(group) == 1:
            cluster_index += 1
            index = group.index[0]
            df.at[index, 'cluster'] = cluster_index
            group.at[index, 'cluster'] = cluster_index
    return df
    
# Function to extract chromosome, start, and end positions from junction
def extract_chr_start_end(junction):
    chr, start, end = junction.split(':')
    return chr, int(start), int(end)

# Function to modify the index
def modify_index(index):
    return ':'.join(index.split(':')[:3])

def get_exon_data(junction_sum, points_df):
    print('get_exon_data')
    # Preprocess points_df into dictionaries for fast lookup
    rank_h_dict = dict(zip(points_df['position_h'], points_df['rank_h']))
    rank_m_dict = dict(zip(points_df['position_m'], points_df['rank_m']))
    symbol_dict = dict(zip(points_df['position_h'], points_df['symbol']))
    ensembl_dict = dict(zip(points_df['position_h'], points_df['ensembl']))
    # Prepare lists to collect data
    rank_h_list = []
    rank_m_list = []
    symbol_h_list = []
    ensembl_h_list = []
    for index, row_data in junction_sum.iterrows():
        # if index % 1000 == 0:
        #     print(index, "get exon data")
        h_chr, h_start, h_end = row_data['h_junction'].split(':')
        m_chr, m_start, m_end = row_data['m_junction'].split(':')
        h_start_key = f"{h_chr}:{h_start}"
        h_end_key = f"{h_chr}:{h_end}"
        m_start_key = f"{m_chr}:{m_start}"
        m_end_key = f"{m_chr}:{m_end}"
        # Fast dictionary lookups
        rank_h_start = rank_h_dict.get(h_start_key, '*')
        rank_h_end = rank_h_dict.get(h_end_key, '*')
        rank_m_start = rank_m_dict.get(m_start_key, '*')
        rank_m_end = rank_m_dict.get(m_end_key, '*')
        rank_h_list.append(f"{rank_h_start}_{rank_h_end}")
        rank_m_list.append(f"{rank_m_start}_{rank_m_end}")
        symbol_h = symbol_dict.get(h_start_key) or symbol_dict.get(h_end_key) or ''
        ensembl_h = ensembl_dict.get(h_start_key) or ensembl_dict.get(h_end_key) or ''
        symbol_h_list.append(symbol_h)
        ensembl_h_list.append(ensembl_h)
    # Assign new columns
    junction_sum['rank_h'] = rank_h_list
    junction_sum['rank_m'] = rank_m_list
    junction_sum['symbol_h'] = symbol_h_list
    junction_sum['ensembl_h'] = ensembl_h_list
    return junction_sum

def duplicated_junctions(df, junction_dict, speciesName):
    species_dict = {"h": "m","m": "h"}
    mainJunction = speciesName + "_junction"
    mainRank = "rank_" + speciesName
    secondRank = "rank_" + species_dict[speciesName]
    # Dictionary to track negative replacements
    int_values = {k: int(v) for k, v in junction_dict.items()}
    negative_counter = min(int_values.values()) - 1
    # Group by h_junction to find duplicates
    duplicates = df[df.duplicated(mainJunction, keep=False)]
    problemIndex = 0
    while not duplicates.empty and problemIndex <= 10:
        problemIndex += 1
        # Process each group of duplicated h_junction
        for junction_value, group in duplicates.groupby(mainJunction):
            # Count stars in rank_h_combined and rank_m_combined
            star_counts = group.apply(lambda row: row[mainRank].count('*') + row[secondRank].count('*'), axis=1)
            # Identify the row with the most stars
            max_star_index = star_counts.idxmax()
            row = df.loc[max_star_index]
            # Determine which part of h_junction to replace based on rank_m_combined        
            first_junction_value = df.at[max_star_index, mainJunction]
            first_parts = str(first_junction_value).split(':')  # Ensure it's a string
            rank_second_value = row[secondRank]
            rank_first_value = df.at[max_star_index, mainRank]
            rank_first_split = str(rank_first_value).split('_')
            # Assign a unique negative value
            first_junction_value = df.at[max_star_index, mainJunction]
            rank_second_value = df.at[max_star_index, secondRank]
            #replacement_key = (first_junction_value, rank_second_value)
            #if replacement_key not in junction_dict:
            #    junction_dict[replacement_key] = negative_counter
            #    negative_counter -= 1
            #neg_value = junction_dict[replacement_key]
            # Replace start or end based on position of '*'
            if rank_second_value.startswith('*'):
                replacement_key = f"{first_parts[0]}:{first_parts[1]}"
                if replacement_key not in junction_dict:
                    junction_dict[replacement_key] = negative_counter
                    negative_counter -= 1
                new_junction = f"{first_parts[0]}:{junction_dict[replacement_key]}:{first_parts[2]}"
                new_rank = f"*_{rank_first_split[1]}"
            if rank_second_value.endswith('*'):
                replacement_key = f"{first_parts[0]}:{first_parts[2]}"
                if replacement_key not in junction_dict:
                    junction_dict[replacement_key] = negative_counter
                    negative_counter -= 1
                new_junction = f"{first_parts[0]}:{first_parts[1]}:{junction_dict[replacement_key]}"
                new_rank = f"{rank_first_split[0]}_*"
    
            # Update the DataFrame
            df.at[max_star_index, mainJunction] = new_junction
            df.at[max_star_index, mainRank] = new_rank
        duplicates = df[df.duplicated(mainJunction, keep=False)]
    if problemIndex >= 10:
        print("problem in duplication", problemIndex)
    return(df)

def main():
    minReads = 10 #remove junctions with less than 5 reads in all the samples
    main_dir = '/gpfs0/tals/projects/Analysis/human_mouse_exons/'
    if len(sys.argv) != 2:
        print("Usage: python get_input.py <value>")
        sys.exit(1)
    version = sys.argv[1] #version of the results
    #version = "HN4"
    ensembl_dir = main_dir + 'ensembl115/'
    points_input = ensembl_dir + 'unique_points_' + version + '.txt' #all the ortholgs points
    points_df = pd.read_csv(points_input, sep='\t',index_col=0)
    #human_mouse
    group_1 = 'GSE115736'
    group_2 = 'GSE116177'
    species_1 = "h"
    species_2 = "m"
    #human_mouse fibroblast 
    #group_1 = 'GSE121052'
    #group_2 = 'GSE161648'
    #human_mouse fibroblast E-MTAB-5919
    #group_1 = 'E-MTAB-5919_human'
    #group_2 = 'E-MTAB-5919_mouse'
    #human_human
    # group_1 = 'GSE115736'
    # group_2 = 'GSE60424'
    # species_1 = "h"
    # species_2 = "h"
    #mouse mouse
    # group_1 = 'GSE116177'
    # group_2 = 'GSE180020'
    #fibroblast_dir = '/gpfs0/tals/projects/Analysis/human_mouse_exons/fibroblast/'
    #group_1_input = fibroblast_dir + group_1 + "/leafcutter_0.2.9/" + group_1 + "_JSR_junction_counts.tsv"
    #group_1_input = fibroblast_dir + group_1[0:11] + "/leafcutter_0.2.9/" + group_1 + "_JSR_junction_counts.tsv"
    group_1_input = main_dir + group_1 + "/leafcutter_0.2.9/" + group_1 + "_JSR_junction_counts.tsv"
    group_1_df = pd.read_csv(group_1_input, sep='\t', index_col=0)
    print(len(group_1_df), "junctions from ", group_1)
    group_1_df = group_1_df[(group_1_df >= minReads).any(axis=1)]
    group_1_df.index = group_1_df.index.map(modify_index)
    group_1_df.columns = [col + "_" + group_1 for col in group_1_df.columns] #add sample name to col
    print(len(group_1_df), "express junctions from ", group_1)
    #group_1_df = group_1_df.head(2000)
    plot_histogram(group_1_df, group_1, 'immune', ensembl_dir) #plot to choose the min JSR of all the cells
    
    #group_2_input = fibroblast_dir + group_2 + "/leafcutter_0.2.9/" + group_2 + "_JSR_junction_counts.tsv"
    #group_2_input = fibroblast_dir + group_2[0:11] + "/leafcutter_0.2.9/" + group_2 + "_JSR_junction_counts.tsv"
    group_2_input = main_dir + group_2 + "/leafcutter_0.2.9/" + group_2 + "_JSR_junction_counts.tsv"
    group_2_df = pd.read_csv(group_2_input, sep='\t', index_col=0)
    print(len(group_2_df), "junctions from ", group_2)
    group_2_df = group_2_df[(group_2_df >= minReads).any(axis=1)]
    group_2_df.index = group_2_df.index.map(modify_index)
    group_2_df.columns = [col + "_" + group_2 for col in group_2_df.columns] #add sample name to col
    print(len(group_2_df), "express junctions from ", group_2)
    #group_2_df = group_2_df.head(2000)
    plot_histogram(group_2_df, group_2, 'immune', ensembl_dir) 
    if species_1 == "h":
    #get the rows that have ortholgs
        orthologs_group_1, junction_dict_h = get_rows(points_df, group_1_df, 'position_h')
        output_file_dict = ensembl_dir + "/junction_dict_h.json"        
        # Save to a file
        with open(output_file_dict, "w") as f:
            json.dump(junction_dict_h, f)
        orthologs_group_1.index.name = 'h_junction'
    #else:
    #    orthologs_group_1 = get_rows(points_df, group_1_df, 'position_m')
    print(len(orthologs_group_1), "orthologs junctions from ", group_1)
    output_file = ensembl_dir + "/orthologous_junctions_" + group_1 + ".txt"
    #orthologs_group_1 = pd.read_csv(output_file, sep='\t', index_col=0)
    orthologs_group_1.to_csv(output_file, sep="\t")
    # if species_2 == "h":
    #     orthologs_group_2 = get_rows(points_df, group_2_df, 'position_h') 
    #     orthologs_group_2_human_index = orthologs_group_2
    #     orthologs_group_2_human_index = orthologs_group_2_human_index.rename_axis('h_junction')
    if species_2 == "m": #mouse
        orthologs_group_2, junction_dict_m = get_rows(points_df, group_2_df, 'position_m')
        output_file_dict = ensembl_dir + "junction_dict_m.json"        
        # Save to a file
        with open(output_file_dict, "w") as f:
            json.dump(junction_dict_m, f)
        #human positions are the index        
        orthologs_group_2['m_junction'] = orthologs_group_2.index
        orthologs_group_2.set_index('h_junction', inplace=True)        
        # Reorder columns so 'mouse_junction' is last
        #cols = [col for col in orthologs_group_2.columns if col != 'mouse_junction'] + ['mouse_junction']
        #orthologs_group_2 = orthologs_group_2[cols]
       
    print(len(orthologs_group_2), "orthologs junctions from ", group_2)
    output_file = ensembl_dir + "orthologous_junctions_" + group_2 + ".txt"
    orthologs_group_2.to_csv(output_file, sep="\t")
       
    # Select only the relevant columns from each DataFrame
    df1_selected = orthologs_group_1.reset_index()[['h_junction', 'm_junction']]
    df2_selected = orthologs_group_2.reset_index()[['h_junction', 'm_junction']]    
    # Concatenate the two DataFrames
    junction_sum = pd.concat([df1_selected, df2_selected])  
    junction_sum = junction_sum.drop_duplicates()
    print(len(junction_sum), "merge junctions - union")    
    #junction_sum = junction_sum.head(2000)
    junction_sum = get_exon_data(junction_sum, points_df) #add exons data
    junction_sum.reset_index(drop=True, inplace=True)
    junction_sum_human = duplicated_junctions(junction_sum, junction_dict_h, "h") #deal with same junctions in human that are different in mouse
    junction_sum_human_mouse = duplicated_junctions(junction_sum_human, junction_dict_m, "m") #deal with same junctions in human that are different in mouse
    
    merged_df = junction_sum_human_mouse.merge(group_1_df, how='left', left_on='h_junction', right_index=True)
    merged_df = merged_df.merge(group_2_df, how='left', left_on='m_junction', right_index=True)
    merged_df = merged_df.fillna(0)
    output_file = ensembl_dir + "junctions_sum_" + group_1 + "_" + group_2 + ".txt"
    merged_df.to_csv(output_file, sep="\t")
    
    #cluster the juncctions from human. i don't need this step here
    #junction_sum_final = assign_clusters_to_groups(junction_sum_human_mouse)
    #output_file = main_dir + "ensembl113/junctions_B/clusters_sum_" + group_1 + "_" + group_2 + ".txt"
    #junction_sum_final.to_csv(output_file, sep="\t")
    return
          
if __name__ == "__main__":
      main() 