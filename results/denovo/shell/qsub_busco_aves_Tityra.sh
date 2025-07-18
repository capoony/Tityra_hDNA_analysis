#!/usr/bin/env bash

## name of Job
#PBS -N BUSCO_Tityra

## Redirect output stream to this file.
#PBS -o /media/inter/mkapun/projects/Tityra/results/denovo/log/Busco_Tityra_log.txt

## Stream Standard Output AND Standard Error to outputfile (see above)
#PBS -j oe

## Select 150 cores and 200gb of RAM
#PBS -l select=1:ncpus=150:mem=200g

## Go to pwd
cd /media/inter/pipelines/AutDeNovo

######## load dependencies #######

eval "$(conda shell.bash hook)"
conda activate envs/busco

######## run analyses #######

mkdir -p /media/inter/mkapun/projects/Tityra/results/denovo/results/AssemblyQC/Busco_aves

cd /media/inter/mkapun/projects/Tityra/results/denovo/results/AssemblyQC/Busco_aves

busco -i ../../../output/Tityra_ILL.fa -o Tityra -m genome -c 150 -f -l aves_odb10
