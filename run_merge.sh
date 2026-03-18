#!/bin/sh -x
conda activate ../../../../software/Anaconda3-2025.06/envs/spyder_HN1/
#python /gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl115/merge_leafcutter_group1_group2_HN23.py 'HN6'
python /gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl115/merge_leafcutter_group1_group2_HN23_human_mouse_optimized.py 'HN6'
