#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 09:54:17 2025

@author: nergaon
"""

import requests
from Bio.Align import PairwiseAligner

def seq(species, chromosome, start, end):
    # Ensembl REST API URL
    url = f"https://rest.ensembl.org/sequence/region/{species}/{chromosome}:{start}..{end}:1?"
    # Request headers
    headers = {"Content-Type": "text/plain"}
    # Make the request
    response = requests.get(url, headers=headers)
    # Output only the sequence
    if response.ok:
        sequence = response.text.strip()
    else:
        print(f"Error: {response.status_code}")   
    return(sequence)

def alignment(seq_h, seq_m):
    aligner = PairwiseAligner() 
    aligner.mode = 'global' #local or global
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -12
    aligner.extend_gap_score = -4
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    alignments = aligner.align(seq_h, seq_m)
    alignment = alignments[0]
    identities = sum(aa1 == aa2 for aa1, aa2 in zip(alignment[0], alignment[1]))
    percentage_identity_h = round((identities / len(seq_h)) * 100)
    percentage_identity_m = round((identities / len(seq_m)) * 100)
    return(percentage_identity_h, percentage_identity_m)

def main():
    species = "human"
    chromosome = "10"
    start = 119576611
    end = 119577203
    seq_h = seq(species, chromosome, start, end)
    species = "mus_musculus"   
    chromosome = "7"
    start = 128045531
    end = 128046128
    seq_m = seq(species, chromosome, start, end)
    return

if __name__ == "__main__":
      main() 