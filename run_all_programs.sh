echo "run_merge"
qsub -cwd -V -q tals.q run_merge.sh
while ('qstat | grep run_merge | wc -l' > 0)
    sleep 10
end
#qsub -cwd -V -q tals.q mkdir_HN1.sh
while ('qstat | grep mkdir_HN1 | wc -l' > 0)
    sleep 10
end
echo "runImmune"
qsub -cwd -V -q tals.q runImmune.sh
while ('qstat | grep runImmune | wc -l' > 0)
    sleep 10
end
echo "leafcutter"
qsub -cwd -V -q tals.q leafcutter_0.2.9_HN4.sh
while ('qstat | grep leafcutter_0.2.9_HN4 | wc -l' > 0)
    sleep 10
end
echo "/merge_leafcutter_cells"
python /gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl115/merge_leafcutter_cells_HN12.py 'HN4'
rm -f /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/leafcutter_GSE115736_GSE116177/genes_figs/*png
rm -f /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl115/leafcutter_GSE115736_GSE116177/genes_figs/*svg
echo "runPlot"
qsub -cwd -V -q tals.q runPlot.sh