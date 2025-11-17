#!/bin/sh
#get orthologs data
input_file="/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/all/good_ortologs_df.txt"
line_count=$(wc -l < "$input_file")
echo $line_count
for i in $(seq 0 $line_count);
        do
                r=$(qstat -u nergaon|wc -l)
                while [ $r -ge 15 ] #if i have more than 15 jobs the program will run only 15 and wait. ge-greater or equal
                do
                        r=$(qstat -u nergaon|wc -l)
                        sleep 3
                done
                echo $i
                qsub -cwd -q tals.q -V run_transcripts_exons_single.sh $i
done


