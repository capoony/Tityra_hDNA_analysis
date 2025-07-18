
  #!/usr/bin/env bash

  ## name of Job
  #PBS -N mapping_Tityra

  ## Redirect output stream to this file.
  #PBS -o /media/inter/mkapun/projects/Tityra/results/denovo/log/mapping_Tityra_log.txt

  ## Stream Standard Output AND Standard Error to outputfile (see above)
  #PBS -j oe

  ## Select 150 cores and 200gb of RAM
  #PBS -l select=1:ncpus=150:mem=200g

  ## Go to pwd
  cd /media/inter/pipelines/AutDeNovo

  ######## load dependencies #######

  eval "$(conda shell.bash hook)"

  ######## run analyses #######

  if [[ ( ILL == 'ILL' ) || ( ILL == 'ILL_ONT' ) || ( ILL == 'ILL_PB' ) || ( ILL == 'ILL_ONT_PB' ) ]]
  then

    conda activate envs/bwa

    ## index reference
    bwa index /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa

    if [[ -f /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2.fq.gz ]]; then
      bwa mem       -t 150       /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa       /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1.fq.gz       /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2.fq.gz       | samtools view -bh | samtools sort       > /media/inter/mkapun/projects/Tityra/results/denovo/results/mapping/Tityra.bam
    else
      bwa mem       -t 150       /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa       /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1.fq.gz       | samtools view -bh | samtools sort       > /media/inter/mkapun/projects/Tityra/results/denovo/results/mapping/Tityra.bam
    fi

  elif [[ ( ILL == 'ONT' ) || ( ILL == 'ONT_PB' ) ]]
  then

    conda activate envs/minimap

    ## index reference
    minimap2 -d /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.mmi       /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa

    minimap2 -ax map-ont     -t 150     /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa     /media/inter/mkapun/projects/Tityra/results/denovo/data/ONT/.fq.gz     | samtools view -bh | samtools sort     > /media/inter/mkapun/projects/Tityra/results/denovo/results/mapping/Tityra.bam

  else

    conda activate envs/minimap

    ## index reference
    minimap2 -d /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.mmi       /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa

    minimap2 -ax map-pb     -t 150     /media/inter/mkapun/projects/Tityra/results/denovo/output/Tityra_ILL.fa     /media/inter/mkapun/projects/Tityra/results/denovo/data/PB/.fq.gz     | samtools view -bh | samtools sort     > /media/inter/mkapun/projects/Tityra/results/denovo/results/mapping/Tityra.bam
  fi


