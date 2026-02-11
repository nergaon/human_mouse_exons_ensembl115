#!/bin/tcsh -x
source /gpfs0/tals/projects/software/Anaconda3-2025.06/etc/profile.d/conda.sh
conda activate /gpfs0/tals/projects/software/Anaconda3-2025.06/envs/spyder_HN1

echo "run_sum_data"
#python /gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl115/sum_data_HN12.py "HN6"

echo "run_merge"
python /gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl115/merge_leafcutter_group1_group2_HN23.py "HN6"

#echo "mkdir_HN1"
#qsub -cwd -V -q tals.q mkdir_HN1.sh
#sleep 10

echo "runImmune"
python /gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl115/immune_HN9.py "HN6"

echo "leafcutter"
qsub -cwd -V -q tals.q leafcutter_0.2.9_HN5.sh "HN6"
r=$(qstat | grep leafcutter_0.2.9_HN4 | wc -l)
while [ $r -g 0 ] #if i have more than 15 jobs the program will run only 15 and wait. ge-greater or equal
do
    r=$(qstat | grep leafcutter_0.2.9_HN4 | wc -l)
    sleep 5
done

echo "merge_leafcutter_cells"
python /gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl115/merge_leafcutter_cells_HN12.py "HN6"

echo "new junctions"
python /gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl115/new_exons_HN2.py "HN6"

echo "runPlot"
rm -f /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/leafcutter_GSE115736_GSE116177/genes_figs/*png
rm -f /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/leafcutter_GSE115736_GSE116177/genes_figs/*svg

python /gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl115/plot_bar_HN11.py "HN6"
