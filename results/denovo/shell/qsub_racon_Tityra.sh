

  #!/usr/bin/env bash

  ## name of Job
  #PBS -N Racon_Tityra

  ## Redirect output stream to this file.
  #PBS -o /media/inter/mkapun/projects/Tityra/results/denovo/log/Racon_Tityra_log.txt

  ## Stream Standard Output AND Standard Error to outputfile (see above)
  #PBS -j oe

  ## Select 150 cores and 1000gb of RAM
  #PBS -l select=1:ncpus=150:mem=1000g

  ## Go to pwd
  cd /media/inter/pipelines/AutDeNovo

  ######## load dependencies #######

  eval "$(conda shell.bash hook)"
  conda activate envs/pigz

  ## concatenate Illumina data (if needed)

  if [[ ILL == *'ILL' ]]
  then
    pigz -p 150 -dc /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1.fq.gz       | sed 's/ 1:.*/\/1/g'       | pigz -p 150 > /media/inter/mkapun/projects/Tityra/results/denovo/results/Racon/Ill.fq.gz

    if [[ -f /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2.fq.gz ]]; then
      pigz -p 150 -dc /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2.fq.gz         | sed 's/ 2:.*/\/2/g'         | pigz -p 150 >> /media/inter/mkapun/projects/Tityra/results/denovo/results/Racon/Ill.fq.gz
    fi
  fi

  ## make copy of unpolished contigs

  cp /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL_unpolished.fa
  pigz -p 150 /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL_unpolished.fa

  ## do Racon polishing
  for (( i=1; i<=4; i++ ))

  do

    if [[ ILL == 'ILL' ]]
    then

      conda deactivate
      conda activate envs/minimap

      minimap2         -x sr         -t 150         /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa         /media/inter/mkapun/projects/Tityra/results/denovo/results/Racon/Ill.fq.gz         > /media/inter/mkapun/projects/Tityra/results/denovo/results/Racon/temp_reads_to_draft.paf

      conda deactivate
      conda activate envs/racon

      racon         -t 150         /media/inter/mkapun/projects/Tityra/results/denovo/results/Racon/Ill.fq.gz         /media/inter/mkapun/projects/Tityra/results/denovo/results/Racon/temp_reads_to_draft.paf         /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa         > /media/inter/mkapun/projects/Tityra/results/denovo/results/Racon/temp_draft_new.fa

    elif [[ ILL == *'ONT'* ]]
    then

      conda deactivate
      conda activate envs/minimap
      
      minimap2         -x map-ont         -t 150         /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa         /media/inter/mkapun/projects/Tityra/results/denovo/data/ONT/.fq.gz         > /media/inter/mkapun/projects/Tityra/results/denovo/results/Racon/temp_reads_to_draft.paf

      conda deactivate
      conda activate envs/racon
      
      racon         -t 150         /media/inter/mkapun/projects/Tityra/results/denovo/data/ONT/.fq.gz         /media/inter/mkapun/projects/Tityra/results/denovo/results/Racon/temp_reads_to_draft.paf         /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa         > /media/inter/mkapun/projects/Tityra/results/denovo/results/Racon/temp_draft_new.fa

    else

      conda deactivate
      conda activate envs/minimap
      
      minimap2         -x map-pb         -t 150         /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa         /media/inter/mkapun/projects/Tityra/results/denovo/data/PB/.fq.gz         > /media/inter/mkapun/projects/Tityra/results/denovo/results/Racon/temp_reads_to_draft.paf

      conda deactivate
      conda activate envs/racon
      
      racon         -t 150         /media/inter/mkapun/projects/Tityra/results/denovo/data/PB/.fq.gz         /media/inter/mkapun/projects/Tityra/results/denovo/results/Racon/temp_reads_to_draft.paf         /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa         > /media/inter/mkapun/projects/Tityra/results/denovo/results/Racon/temp_draft_new.fa

    fi

    mv /media/inter/mkapun/projects/Tityra/results/denovo/results/Racon/temp_draft_new.fa /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa

  done

  if [[ ILL == 'ILL' ]]
  then
    rm -f /media/inter/mkapun/projects/Tityra/results/denovo/results/Racon/Ill.fq.gz
  fi

  rm -rf /media/inter/mkapun/projects/Tityra/results/denovo/results/Racon/temp_reads_to_draft.paf


