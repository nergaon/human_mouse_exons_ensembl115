#!/bin/sh
#qsub -cwd -V -q tals.q /gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl115/sharedJunctions_2026/leafcutter_0.2.9_HN5.sh 'HN6'
#activate conda env
source /gpfs0/tals/projects/software/anaconda3/etc/profile.d/conda.sh
conda activate leafcutter

version=$1
AS_value="/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/AS_clusters_value_${version}.txt"

#human mouse
leafcutter_folders=(
  "/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/leafcutter_GSE115736_GSE116177/CD4T"
	"/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/leafcutter_GSE115736_GSE116177/CD8T"
	"/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/leafcutter_GSE115736_GSE116177/NveB"
	"/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/leafcutter_GSE115736_GSE116177/NK"
	"/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/leafcutter_GSE115736_GSE116177/Mono"
	"/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/leafcutter_GSE115736_GSE116177/Neut")

#human human
#leafcutter_folders=(
# "/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/leafcutter_GSE115736_GSE60424/CD4T"
#	"/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/leafcutter_GSE115736_GSE60424/CD8T"
#	"/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/leafcutter_GSE115736_GSE60424/NveB"
#	"/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/leafcutter_GSE115736_GSE60424/NK"
#	"/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/leafcutter_GSE115736_GSE60424/Mono"
#	"/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/leafcutter_GSE115736_GSE60424/Neut")

#mouse mouse
#leafcutter_folders=(
#    "/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/leafcutter_GSE116177_GSE180020/CD4T"
#    "/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/leafcutter_GSE116177_GSE180020/Cd8T"
#    "/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/leafcutter_GSE116177_GSE180020/BCell"
#    "/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/leafcutter_GSE116177_GSE180020/NK"
#    "/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/sharedJunctions_2026/leafcutter_GSE116177_GSE180020/Mono")

leafcutter_scripts="/gpfs0/tals/projects/software/leafcutter_0.2.9"

#Differential intron excision analysis
#python creates the group file
#change -i and -g
# Loop through each folder and perform the analysis
for leafcutter_folder in "${leafcutter_folders[@]}"; do
	cd $leafcutter_folder
  # Set the path to your groups file
  GROUPS_FILE="groups_file.txt"
  # Extract the smallest group size
  min_group_size=$(cut -f2 "$GROUPS_FILE" | sort | uniq -c | awk '{print $1}' | sort -n | head -n1)
# Calculate 60% of the smallest group size (rounded down)
  threshold=$(echo "$min_group_size * 60/100" | bc)
# Enforce minimum threshold of 2
  if [ "$threshold" -lt 2 ]; then
    threshold=2
  fi
	$leafcutter_scripts/leafcutter/scripts/leafcutter_ds.R --num_threads 4 -i "$threshold"  -g "$threshold" --exon_file=/gpfs0/tals/projects/data/Transcriptomes/human_hg38/GCF_000001405.39_GRCh38.p13_genomic_leafcutter.txt.gz "$AS_value" groups_file.txt
done

# Program:        regtools
# Version:        0.5.2
# Usage:          regtools junctions extract [options] indexed_alignments.bam
# Options:
# -a INT  Minimum anchor length. Junctions which satisfy a minimum
# anchor length on both sides are reported. [8]
# -m INT  Minimum intron length. [70]
# -M INT  Maximum intron length. [500000]
# -o FILE The file to write output to. [STDOUT]
# -r STR  The region to identify junctions
# in "chr:start-end" format. Entire BAM by default.
# -s INT  Strand specificity of RNA library preparation
# (0 = unstranded, 1 = first-strand/RF, 2, = second-strand/FR). REQUIRED
# -t STR  Tag used in bam to label strand. [XS]

#leafcutter_cluster_regtools.py
#-j, --juncfiles, text file with all junction files to be processed
#-o, --outprefix, output prefix (default leafcutter))
#-q, --quiet,don't print status messages to stdout (default=True)
#-r, --rundir, write to directory (default ./))
#-l, --maxintronlen, maximum intron length in bp (default 100,000bp)
#-m, --minclureads, minimum reads in a cluster (default 30 reads)
#-p, --mincluratio, minimum fraction of reads in a cluster that support a junction (default 0.001)
#-c, --cluster, refined cluster file when clusters are already made, (default = None)
#-k, --nochromcheck, Don't check that the chromosomes are well formated e.g. chr1, chr2, ..., or 1, 2, ..., (default = False)
#-C, --includeconst, also include constitutive introns, (default = False)

#leafcutter_ds.R
#-o,"--output_prefix", default = "leafcutter_ds", The prefix for the two output files, <prefix>_cluster_significance.txt (containing test status, log likelihood ratio, degree of freedom, and p-value for each cluster) and <prefix>_effect_sizes.txt (containing the effect sizes for each intron
#-s,"--max_cluster_size", default=Inf, Don't test clusters with more introns than this
#-i,"--min_samples_per_intron", default=5, Ignore introns used (i.e. at least one supporting read) in fewer than n samples 
#-g,"--min_samples_per_group", default=3, Require this many samples in each group to have at least min_coverage reads 
#-c,"--min_coverage", default=20, Require min_samples_per_group samples in each group to have at least this many reads 
#-t,"--timeout", default=30, Maximum time (in seconds) allowed for a single optimization run
#-p,"--num_threads", default=1, Number of threads to use
#-e,"--exon_file", default=NULL, File defining known exons, example in data/gencode19_exons.txt.gz. Columns should be chr, start, end, strand, gene_name. Optional, only just to label the clusters
#--init, One of 'smart' (default) or 'random'. 
#--seed, default=12345, Random seed if using random initialization.
