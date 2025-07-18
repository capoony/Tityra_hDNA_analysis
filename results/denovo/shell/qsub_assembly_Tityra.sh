
  #!/usr/bin/env bash

  ## name of Job
  #PBS -N denovo_Tityra

  ## Redirect output stream to this file.
  #PBS -o /media/inter/mkapun/projects/Tityra/results/denovo/log/assembly_Tityra_log.txt

  ## Stream Standard Output AND Standard Error to outputfile (see above)
  #PBS -j oe

  ## Select a maximum walltime of 2h
  #PBS -l walltime=100:00:00

  ## Select 8 cores and 1000gb of RAM
  #PBS -l select=1:ncpus=8:mem=1000g

  ## Go to pwd
  cd /media/inter/pipelines/AutDeNovo

  ######## load dependencies #######

  eval "$(conda shell.bash hook)"

  if [[ ( ILL == 'ILL' ) ]]
  then

    conda activate envs/spades

    if [[ -f /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2.fq.gz ]]; then
      spades.py       -1 /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1.fq.gz       -2 /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2.fq.gz       -t 8       -m 1000       -o /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra
    else
      spades.py       -s /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1.fq.gz       -t 8       -m 1000       -o /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra
    fi

    if [[ -f /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra/scaffolds.fasta ]]
    then
      mv /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra/scaffolds.fasta /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa
    else
      mv /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra/contigs.fasta /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa
    fi

  elif [[ ( ILL == 'ILL_ONT' ) ]]
  then

    conda activate envs/spades

    if [[ -f /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2.fq.gz ]]; then
      spades.py       -1 /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1.fq.gz       -2 /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2.fq.gz       --nanopore /media/inter/mkapun/projects/Tityra/results/denovo/data/ONT/.fq.gz       -t 8       -m 1000       -o /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra
    else
      spades.py       -s /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1.fq.gz       --nanopore /media/inter/mkapun/projects/Tityra/results/denovo/data/ONT/.fq.gz       -t 8       -m 1000       -o /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra
    fi

    if [[ -f /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra/scaffolds.fasta ]]
    then
      mv /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra/scaffolds.fasta /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa
    else
      mv /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra/contigs.fasta /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa
    fi

  elif [[ ( ILL == 'ILL_PB' ) ]]
  then

    conda activate envs/spades

    if [[ -f /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2.fq.gz ]]; then
      spades.py       -1 /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1.fq.gz       -2 /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2.fq.gz       --pacbio /media/inter/mkapun/projects/Tityra/results/denovo/data/PB/.fq.gz       -t 8       -m 1000       -o /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra
    else
      spades.py       -s /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1.fq.gz       --pacbio /media/inter/mkapun/projects/Tityra/results/denovo/data/PB/.fq.gz       -t 8       -m 1000       -o /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra
    fi

    if [[ -f /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra/scaffolds.fasta ]]
    then
      mv /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra/scaffolds.fasta /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa
    else
      mv /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra/contigs.fasta /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa
    fi

  elif [[ ( ILL == 'ILL_ONT_PB' ) ]]
  then

    conda activate envs/spades

    if [[ -f /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2.fq.gz ]]; then
      spades.py       -1 /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1.fq.gz       -2 /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2.fq.gz       --pacbio /media/inter/mkapun/projects/Tityra/results/denovo/data/PB/.fq.gz       --nanopore /media/inter/mkapun/projects/Tityra/results/denovo/data/ONT/.fq.gz       -t 8       -m 1000       -o /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra
    else
      spades.py       -s /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1.fq.gz       --pacbio /media/inter/mkapun/projects/Tityra/results/denovo/data/PB/.fq.gz       --nanopore /media/inter/mkapun/projects/Tityra/results/denovo/data/ONT/.fq.gz       -t 8       -m 1000       -o /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra
    fi

    if [[ -f /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra/scaffolds.fasta ]]
    then
      mv /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra/scaffolds.fasta /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa
    else
      mv /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra/contigs.fasta /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa
    fi

  elif [[ ( ILL == 'ONT' ) ]]
  then

    conda activate envs/flye

    ## Flye only accepts 128 threads max
    if [[ 8 -gt 128 ]]
    then
      threads=128
    fi

    flye     --nano-raw /media/inter/mkapun/projects/Tityra/results/denovo/data/ONT/.fq.gz     --out-dir /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra     --threads 8     --scaffold

    mv /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra/assembly.fasta /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa

  elif [[ ( ILL == 'PB' ) ]]
  then

    conda activate envs/flye

    ## Flye only accepts 128 threads max
    if [[ 8 -gt 128 ]]
    then
      threads=128
    fi

    flye     --pacbio-raw /media/inter/mkapun/projects/Tityra/results/denovo/data/PB/.fq.gz     --out-dir /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra     --threads 8     --scaffold

    mv /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra/assembly.fasta /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa

  elif [[ ( ILL == 'ONT_PB' ) ]]
  then

    conda activate envs/flye

    ## see here: https://github.com/fenderglass/Flye/blob/flye/docs/FAQ.md#can-i-use-both-pacbio-and-ont-reads-for-assembly

    ## Flye only accepts 128 threads max
    if [[ 8 -gt 128 ]]
    then
      threads=128
    fi

    flye     --pacbio-raw /media/inter/mkapun/projects/Tityra/results/denovo/data/ONT/.fq.gz     /media/inter/mkapun/projects/Tityra/results/denovo/data/PB/.fq.gz     --iterations 0     --out-dir /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra     --threads 8

    flye     --nano-raw /media/inter/mkapun/projects/Tityra/results/denovo/data/ONT/.fq.gz      --resume-from polishing     --out-dir /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra     --threads 8     --scaffold

    mv /media/inter/mkapun/projects/Tityra/results/denovo/results/assembly/Tityra/assembly.fasta /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa
  fi

