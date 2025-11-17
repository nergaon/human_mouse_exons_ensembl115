#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 10:53:00 2024

@author: nergaon
"""

import requests
from biomart import BiomartServer
import pandas as pd

def main():
    working_dir = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/all/'
    # Define the URL for the Ensembl REST API
    url = "https://rest.ensembl.org/info/software?content-type=application/json"
    
    # Make a GET request to the API
    response = requests.get(url)
    
    # Check if the request was successful
    if response.status_code == 200:
        data = response.json()
        ensembl_version = data.get("release")
        if ensembl_version:
            print(f"Latest Ensembl version: {ensembl_version}")
        else:
            print("Could not retrieve Ensembl version.")
    else:
        print(f"Failed to connect to Ensembl API. Status code: {response.status_code}")
    
    # Connect to Ensembl Biomart server using the latest version
    server = BiomartServer("http://www.ensembl.org/biomart")
    #Choose the Ensembl Genes version
    #server = BiomartServer(f"http://www.ensembl.org/biomart?release={ensembl_version}")
    
    #Choose Mouse genes dataset
    mart = server.datasets["mmusculus_gene_ensembl"]
    
    # Set the attributes you want to retrieve
    attributes = [
        "ensembl_gene_id",
        "hsapiens_homolog_ensembl_gene",
        "hsapiens_homolog_associated_gene_name",
        "hsapiens_homolog_orthology_type",
        "hsapiens_homolog_perc_id",
        "hsapiens_homolog_perc_id_r1",
        "hsapiens_homolog_orthology_confidence"
    ]
    
    # Perform the query
    response = mart.search({
        'attributes': attributes
    })
    
    # Process the results into a list of dictionaries
    results = []
    for line in response.iter_lines():
        row = line.decode('utf-8').split("\t")
        results.append({
            "ensembl_gene_id": row[0],
            "hsapiens_homolog_ensembl_gene": row[1],
            "hsapiens_homolog_associated_gene_name": row[2],
            "hsapiens_homolog_orthology_type": row[3],
            "hsapiens_homolog_perc_id": row[4],
            "hsapiens_homolog_perc_id_r1": row[5],
            "hsapiens_homolog_orthology_confidence": row[6]
        })
    
    # Convert the results to a DataFrame
    mouse_ortologs = pd.DataFrame(results)
    
    # Save the DataFrame to a CSV file
    df_output = working_dir + 'mouse_orthologs_v' + str(ensembl_version) + '.txt'
    mouse_ortologs.to_csv(df_output, index=False)
    unique_gene_ids = mouse_ortologs['ensembl_gene_id'].unique()
    print("there are", len(unique_gene_ids), "genes in mouse")
    filtered_df = mouse_ortologs[mouse_ortologs['hsapiens_homolog_orthology_confidence'] == "1"]
    filtered_df = filtered_df[filtered_df['hsapiens_homolog_orthology_type'] == 'ortholog_one2one']
    unique_gene_ids = filtered_df['ensembl_gene_id'].unique()
    print("there are", len(unique_gene_ids), "good orthologs genes in mouse")
    filtered_df = filtered_df.reset_index(drop=True)
    df_output = working_dir + 'good_ortologs_df.txt'
    filtered_df.to_csv(df_output, sep="\t")
    print("Data saved")
    return
          
if __name__ == "__main__":
      main() 