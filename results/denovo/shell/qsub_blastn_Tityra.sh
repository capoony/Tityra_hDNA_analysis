
  #!/usr/bin/env bash

  ## name of Job
  #PBS -N BLASTN_Tityra

  ## Redirect output stream to this file.
  #PBS -o /media/inter/mkapun/projects/Tityra/results/denovo/log/BLAST_Tityra_log.txt

  ## Stream Standard Output AND Standard Error to outputfile (see above)
  #PBS -j oe

  ## Select 150 cores and 200gb of RAM
  #PBS -l select=1:ncpus=150:mem=200g

  ## Go to pwd
  cd /media/inter/pipelines/AutDeNovo

  ######## load dependencies #######

  eval "$(conda shell.bash hook)"
  conda activate envs/blast

  ######## run analyses #######

  blastn     -num_threads 150     -outfmt "6 qseqid staxids bitscore std"     -max_target_seqs 10     -max_hsps 1     -evalue 1e-25     -db /media/scratch/NCBI_nt_DB_210714/nt     -query /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa     > /media/inter/mkapun/projects/Tityra/results/denovo/results/BLAST/blastn_Tityra.txt


