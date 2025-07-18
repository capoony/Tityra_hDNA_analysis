
# ############### HTML output #####################
# # run the following commands in terminal to open Firefox and view the HTML output files
# 
## Illumina Data - FASTQC of raw reads
firefox --new-tab /media/inter/mkapun/projects/Tityra/results/denovo/results/rawQC/Tityra_Illumina_fastqc/Tityra_1_fastqc.html
firefox --new-tab /media/inter/mkapun/projects/Tityra/results/denovo/results/rawQC/Tityra_Illumina_fastqc/Tityra_2_fastqc.html

## Illumina Data - FASTQC after trimming
firefox --new-tab /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1_fastqc.html
firefox --new-tab /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2_fastqc.html

## QUAST
firefox --new-tab /media/inter/mkapun/projects/Tityra/results/denovo/results/AssemblyQC/Quast/report.html

## Blobtools

conda activate /media/inter/pipelines/AutDeNovo/envs/blobtools

blobtools view   --out /media/inter/mkapun/projects/Tityra/results/denovo/results/AssemblyQC/blobtools/out   --interactive   /media/inter/mkapun/projects/Tityra/results/denovo/results/AssemblyQC/blobtools 
### now copy the URL that is printed in the commandline and paste it in Firefox

