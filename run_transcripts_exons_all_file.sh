#!/bin/bash

# Read the file in the last line and get the gene number
while IFS= read -r line; do
  r=$(qstat -u nergaon|wc -l) 
  while [ $r -ge 15 ] #if i have more than 10 jobs the program will run only 10 and wait. ge-greater or equal 
	do
		r=$(qstat -u nergaon|wc -l) 
		sleep 3 
	done  
	echo "$line"
	qsub -cwd -q tals.q -V run_transcripts_exons_single.sh $line
done < "/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/missing_genes.txt"

