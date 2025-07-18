
  #!/usr/bin/env bash

  ## name of Job
  #PBS -N fastqc_Tityra

  ## Redirect output stream to this file.
  #PBS -o /media/inter/mkapun/projects/Tityra/results/denovo/log/raw_fastqc_Tityra_log.txt

  ## Stream Standard Output AND Standard Error to outputfile (see above)
  #PBS -j oe

  ## Select 150 cores and 200gb of RAM
  #PBS -l select=1:ncpus=2:mem=200g

  ## Go to pwd
  cd /media/inter/pipelines/AutDeNovo

  ######## load dependencies #######

  eval "$(conda shell.bash hook)"
  conda activate envs/fastqc

  mkdir -p /media/inter/mkapun/projects/Tityra/results/denovo/results/rawQC/Tityra_Illumina_fastqc

  ## Go to output folder
  cd /media/inter/mkapun/projects/Tityra/results/denovo/data

  ## loop through all raw FASTQ and test quality

  if [[ -f Illumina/Tityra_2.fq.gz ]]; then
    fastqc       --outdir ../results/rawQC/Tityra_Illumina_fastqc       --threads 2       Illumina/Tityra_1.fq.gz       Illumina/Tityra_2.fq.gz
  else
    fastqc       --outdir ../results/rawQC/Tityra_Illumina_fastqc       --threads 2       Illumina/Tityra_1.fq.gz
  fi


