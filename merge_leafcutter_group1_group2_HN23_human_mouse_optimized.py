#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import json
import sys
import networkx as nx
from general_def import assign_clusters_to_groups

def plot_histogram(df, species, cell, ensembl_dir):
     print("Plotting histogram for species:", species, "cell:", cell)
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
    print("get junctions with at least one point has ortholog", col_name)
    species_dict = {"h": "m", "m": "h"}
    species = col_name[-1]
    second_col = 'position_' + species_dict[species]
    junction_col = species_dict[species] + "_junction"
    # Preprocess points_df into a dictionary for fast lookup
    points_dict = dict(zip(points_df[col_name], points_df[second_col]))
    # Prepare storage
    rows = []
    junction_dict = {}
    # Loop using itertuples for faster iteration
    for row_tuple in df.itertuples(index=True):
        i = row_tuple.Index
        # Handle index format
        if isinstance(i, str):
            chr_, start, end = i.split(':')
        elif isinstance(i, tuple) and len(i) == 3:
            chr_, start, end = i
        else:
            continue  # Skip invalid indices
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
        # Create row copy
        row = pd.DataFrame([row_tuple[1:]], columns=df.columns, index=[i])
        if int(junction_dict[point_a]) > int(junction_dict[point_b]) and int(junction_dict[point_a]) > 0 and int(junction_dict[point_b]) > 0:
            row[junction_col] = f"{chrName}:{junction_dict[point_b]}:{junction_dict[point_a]}"
            # Swap in dict
            temp = junction_dict[point_a]
            junction_dict[point_a] = junction_dict[point_b]
            junction_dict[point_b] = temp
        else:
            row[junction_col] = f"{chrName}:{junction_dict[point_a]}:{junction_dict[point_b]}"
        rows.append(row)
    # Combine rows
    if rows:
        orthologs_group = pd.concat(rows)
    else:
        orthologs_group = pd.DataFrame(columns=df.columns)
    # Modify index
    orthologs_group.index = orthologs_group.index.map(modify_index)
    return orthologs_group, junction_dict
    
# Function to extract chromosome, start, and end positions from junction
def extract_chr_start_end(junction):
    chr, start, end = junction.split(':')
    return chr, int(start), int(end)

def add_exon_data(df, points_df):
    rank_h_dict = dict(zip(points_df['position_h'], points_df['rank_h']))
    symbol_dict = dict(zip(points_df['position_h'], points_df['symbol']))
    ensembl_dict = dict(zip(points_df['position_h'], points_df['ensembl']))
    rank_h_list = []
    symbol_h_list = []
    ensembl_h_list = []
    for index, row_data in df.iterrows():
        h_chr, h_start, h_end = row_data['h_junction'].split(':')
        h_start_key = f"{h_chr}:{h_start}"
        h_end_key = f"{h_chr}:{h_end}"
        rank_h_start = rank_h_dict.get(h_start_key, '*')
        rank_h_end = rank_h_dict.get(h_end_key, '*')
        rank_h_list.append(f"{rank_h_start}_{rank_h_end}")
        symbol_h = symbol_dict.get(h_start_key) or symbol_dict.get(h_end_key) or ''
        symbol_h_list.append(symbol_h)
        ensembl_h = ensembl_dict.get(h_start_key) or ensembl_dict.get(h_end_key) or ''
        ensembl_h_list.append(ensembl_h)
    df['rank_h'] = rank_h_list
    df['symbol_h'] = symbol_h_list
    df['ensembl_h'] = ensembl_h_list
    return df

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
        if int(h_start) < 0 or int(h_end) < 0:
            # Modified h_junction, keep existing symbol_h and ensembl_h, set rank_h to *_*
            rank_h_list.append('*_*')
            symbol_h = row_data['symbol_h']
            ensembl_h = row_data['ensembl_h']
        else:
            rank_h_start = rank_h_dict.get(h_start_key, '*')
            rank_h_end = rank_h_dict.get(h_end_key, '*')
            rank_h_list.append(f"{rank_h_start}_{rank_h_end}")
            symbol_h = symbol_dict.get(h_start_key) or symbol_dict.get(h_end_key) or ''
            ensembl_h = ensembl_dict.get(h_start_key) or ensembl_dict.get(h_end_key) or ''
        rank_m_start = rank_m_dict.get(m_start_key, '*')
        rank_m_end = rank_m_dict.get(m_end_key, '*')
        rank_m_list.append(f"{rank_m_start}_{rank_m_end}")
        symbol_h_list.append(symbol_h)
        ensembl_h_list.append(ensembl_h)
    # Assign new columns
    junction_sum['rank_h'] = rank_h_list
    junction_sum['rank_m'] = rank_m_list
    junction_sum['symbol_h'] = symbol_h_list
    junction_sum['ensembl_h'] = ensembl_h_list
    return junction_sum

def duplicated_junctions(df, junction_dict, speciesName):
    print("duplicated_junctions for species:", speciesName)
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
            #print(junction_value, group)
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
            rank_second_value = df.at[max_star_index, secondRank]
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

def removeJuntions(orthologs_group_2, orthologs_group_1, speciesName):
    """
    In `orthologs_group_2` detect junction pairs where a coordinate is the end of one
    junction and the start of another. Check `orthologs_group_1` for those junctions.
    If one of them is found in `orthologs_group_1`, remove the other one from 
    `orthologs_group_2`.
    Returns a modified copy of `orthologs_group_2` with problematic rows removed 
    (index preserved as `h_junction`).
    """
    print("removeJuntions for species:", speciesName)
    og2 = orthologs_group_2.reset_index().copy()
    og1 = orthologs_group_1.reset_index().copy()
    # determine which column to inspect based on species
    if speciesName == 'm':
        inspect_col = 'h_junction'
    else:
        inspect_col = 'm_junction'
    if inspect_col not in og2.columns or inspect_col not in og1.columns:
        return orthologs_group_2
    # Create a set of all junctions in og1 for fast lookup
    og1_junctions = set(og1[inspect_col].astype(str).dropna().unique())
    # Build maps for og2: start_coords and end_coords
    start_coords = {}  # map: (chr, start) -> [junctions]
    end_coords = {}    # map: (chr, end) -> [junctions]
    for j in og2[inspect_col].astype(str).dropna().unique():
        try:
            ch, s, e = j.split(":")
            s_key = (ch, int(s))
            e_key = (ch, int(e))
            start_coords.setdefault(s_key, []).append(j)
            end_coords.setdefault(e_key, []).append(j)
        except Exception:
            continue
    # Find rows to remove
    rows_to_remove = set()
    # Find pairs where one junction's end == another's start
    for coord in set(start_coords.keys()) & set(end_coords.keys()):
        end_junctions = end_coords.get(coord, [])
        start_junctions = start_coords.get(coord, [])
        for end_j in end_junctions:
            for start_j in start_junctions:
                if end_j == start_j:
                    continue
                # Check if either junction is in og1
                end_j_in_og1 = end_j in og1_junctions
                start_j_in_og1 = start_j in og1_junctions
                # If one is found in og1, remove the other from og2
                if end_j_in_og1 and not start_j_in_og1:
                    # Remove start_j from og2
                    mask = og2[inspect_col].astype(str) == start_j
                    for ridx in og2[mask].index:
                        rows_to_remove.add(ridx)
                elif start_j_in_og1 and not end_j_in_og1:
                    # Remove end_j from og2
                    mask = og2[inspect_col].astype(str) == end_j
                    for ridx in og2[mask].index:
                        rows_to_remove.add(ridx)
    # Remove marked rows
    if rows_to_remove:
        og2 = og2.drop(index=list(rows_to_remove)).reset_index(drop=True)
    # restore index and return
    if 'h_junction' in og2.columns:
        og2.set_index('h_junction', inplace=True)
    return og2

def readDF(group, species, main_dir, ensembl_dir, points_df, version):
    print("Reading data for group:", group, "species:", species)
    group_input = main_dir + group + "/leafcutter_0.2.9/" + group + "_JSR_junction_counts.tsv"
    #group_1_input = main_dir + group_1 + "/leafcutter_0.2.9/" + group_1 + "_JSR_junction_counts_cuta.tsv"
    group_df = pd.read_csv(group_input, sep='\t', index_col=0)
    print(len(group_df), "junctions from ", group)
    group_df.index = group_df.index.map(modify_index)
    group_df.columns = [col + "_" + group for col in group_df.columns] #add sample name to col
    #group_df = group_df.head(2000)
    plot_histogram(group_df, group, 'immune', ensembl_dir) #plot to choose the min JSR of all the cells

    if species == "h":
    #get the rows that have ortholgs
        orthologs_group, junction_dict_h = get_rows(points_df, group_df, 'position_h')
        output_file_dict = ensembl_dir + "/junction_dict_h_" + version + "b.json"        
        # Save to a file
        with open(output_file_dict, "w") as f:
            json.dump(junction_dict_h, f)
        orthologs_group.index.name = 'h_junction'
        junction_dict = junction_dict_h
    if species == "m": #mouse
        orthologs_group, junction_dict_m = get_rows(points_df, group_df, 'position_m')
        output_file_dict = ensembl_dir + "/junction_dict_m_" + version + "b.json"        
        # Save to a file
        with open(output_file_dict, "w") as f:
            json.dump(junction_dict_m, f)
        #human positions are the index        
        orthologs_group['m_junction'] = orthologs_group.index
        orthologs_group.set_index('h_junction', inplace=True)        
        junction_dict = junction_dict_m

    print(len(orthologs_group), "orthologs junctions from ", group)
    output_file = ensembl_dir + "/orthologous_junctions_" + group + "_" + version + "b.txt"
    orthologs_group.to_csv(output_file, sep="\t")
    return group_df, orthologs_group, junction_dict

def main():
    minReads = 10 #remove junctions with less than minReads reads in all min samples
    minSamples = 2 #remove junctions with less than minSamples samples with minReads reads
    min_psi = 0.0 #remove junctions that are less than 5% from all the junctions in a clusters, in all samples
    main_dir = '/gpfs0/tals/projects/Analysis/human_mouse_exons/'
    #if len(sys.argv) != 2:
    #    print("Usage: python get_input.py <value>")
    #    sys.exit(1)
    #version = sys.argv[1] #version of the results
    version = "HN6"
    ensembl_dir = main_dir + 'ensembl115/'
    points_input = ensembl_dir + 'unique_points_' + version + '.txt' #all the ortholgs points
    points_df = pd.read_csv(points_input, sep='\t',index_col=0)
    #human_mouse
    groups = ['GSE115736', 'GSE116177']
    species = ['h', 'm']
    #human_human
    #groups = ['GSE115736', 'GSE60424']
    #species = ['h', 'h']
    #mouse_mouse
    #groups = ['GSE116177', 'GSE180020']
    #species = ['m', 'm']
    for i in range(len(groups)):
        oneGroup = groups[i]
        oneSpecies = species[i]
        print("Processing group:", oneGroup, "species:", oneSpecies, "i:", i)
        if i == 0:
           group_1_df, orthologs_group_1, junction_dict_h = readDF(oneGroup, oneSpecies, main_dir, ensembl_dir, points_df, version)   
        if i == 1:
           group_2_df, orthologs_group_2, junction_dict_m = readDF(oneGroup, oneSpecies, main_dir, ensembl_dir, points_df, version)

    # Add exon data before modifications
    orthologs_group_1 = add_exon_data(orthologs_group_1.reset_index(), points_df).set_index('h_junction')
    orthologs_group_2 = add_exon_data(orthologs_group_2.reset_index(), points_df).set_index('h_junction')

    if species[0] != species[1]: #if the species are different, remove junctions that have problems in human side
        orthologs_group_2 = removeJuntions(orthologs_group_2, orthologs_group_1, speciesName="m") #remove junctions that have problems in human side
        orthologs_group_1 = removeJuntions(orthologs_group_1, orthologs_group_2, speciesName="h")
    # Select only the relevant columns from each DataFrame
    df1_selected = orthologs_group_1.reset_index()[['h_junction', 'm_junction', 'symbol_h', 'ensembl_h', 'rank_h']]
    df2_selected = orthologs_group_2.reset_index()[['h_junction', 'm_junction', 'symbol_h', 'ensembl_h', 'rank_h']]    
    # Concatenate the two DataFrames
    junction_sum = pd.concat([df1_selected, df2_selected])  
    junction_sum = junction_sum.drop_duplicates().reset_index(drop=True)

    print(len(junction_sum), "merge junctions - union")    
    #junction_sum = junction_sum.head(2000)
    junction_sum = get_exon_data(junction_sum, points_df) #add exons data
    junction_sum.reset_index(drop=True, inplace=True)
    if species[0] != species[1]:
        junction_sum_human = duplicated_junctions(junction_sum, junction_dict_h, "h") #deal with same junctions in human that are different in mouse
        junction_sum_human_mouse = duplicated_junctions(junction_sum_human, junction_dict_m, "m") #deal with same junctions in human that are different in mouse
        merged_df = junction_sum_human_mouse.merge(group_1_df, how='left', left_on='h_junction', right_index=True)
        merged_df = merged_df.merge(group_2_df, how='left', left_on='m_junction', right_index=True)
    else:
        merged_df = group_1_df.merge(group_2_df, how='left', left_index=True, right_index=True)
    #remove junctions with less than minSamples samples with minReads reads
    numeric_df = merged_df.select_dtypes(include='number')
    merged_df = merged_df[(numeric_df >= minReads).sum(axis=1) >= minSamples] 
    print(len(merged_df), "merge junctions - after filtering with minReads and minSamples")
    merged_df = merged_df.fillna(0)
    output_file = ensembl_dir + "junctions_merge_" + groups[0] + "_" + groups[1] + "_" + version + ".txt"
    merged_df.to_csv(output_file, sep="\t")
    return
          
if __name__ == "__main__":
      main()