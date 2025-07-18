

  #!/usr/bin/env bash

  ## name of Job
  #PBS -N QUAST_Tityra

  ## Redirect output stream to this file.
  #PBS -o /media/inter/mkapun/projects/Tityra/results/denovo/log/Quast_Tityra_log.txt

  ## Stream Standard Output AND Standard Error to outputfile (see above)
  #PBS -j oe

  ## Select 64 cores and 64gb of RAM
  #PBS -l select=1:ncpus=64:mem=64g

  ## Go to pwd
  cd /media/inter/pipelines/AutDeNovo

  ######## load dependencies #######

  eval "$(conda shell.bash hook)"
  conda activate envs/quast

  ######## run analyses #######

  quast.py   --output-dir /media/inter/mkapun/projects/Tityra/results/denovo/results/AssemblyQC/Quast   --threads 64   --eukaryote   -f   /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa


