#!/usr/bin/env bash

## name of Job
#PBS -N blobtools_Tityra

## Redirect output stream to this file.
#PBS -o /media/inter/mkapun/projects/Tityra/results/denovo/log/blobtools_Tityra_log.txt

## Stream Standard Output AND Standard Error to outputfile (see above)
#PBS -j oe

## Select 150 cores and 200gb of RAM
#PBS -l select=1:ncpus=150:mem=200g

## Go to pwd
cd /media/inter/pipelines/AutDeNovo

######## load dependencies #######

eval "$(conda shell.bash hook)"
conda activate envs/blobtools

mkdir /media/inter/mkapun/projects/Tityra/results/denovo/results/AssemblyQC/blobtools

## create a genome BlobDir
blobtools add --fasta /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa /media/inter/mkapun/projects/Tityra/results/denovo/results/AssemblyQC/blobtools

## add BLAST results
blobtools add --hits /media/inter/mkapun/projects/Tityra/results/denovo/results/BLAST/blastn_Tityra.txt --taxdump /media/scratch/NCBI_taxdump --threads 150 /media/inter/mkapun/projects/Tityra/results/denovo/results/AssemblyQC/blobtools

## add BUSCO results
#blobtools add --busco /media/inter/mkapun/projects/Tityra/results/denovo/results/AssemblyQC/Busco/Tityra/run_vertebrata_odb10/full_table.tsv --threads 150 /media/inter/mkapun/projects/Tityra/results/denovo/results/AssemblyQC/blobtools
blobtools add --busco /media/inter/mkapun/projects/Tityra/results/denovo/results/AssemblyQC/Busco_aves/Tityra/run_aves_odb10/full_table.tsv --threads 150 /media/inter/mkapun/projects/Tityra/results/denovo/results/AssemblyQC/blobtools

## add coverage data
blobtools add --cov /media/inter/mkapun/projects/Tityra/results/denovo/results/mapping/Tityra.bam --threads 150 /media/inter/mkapun/projects/Tityra/results/denovo/results/AssemblyQC/blobtools
