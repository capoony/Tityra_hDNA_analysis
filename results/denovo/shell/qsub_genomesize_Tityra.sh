
  #!/usr/bin/env bash

  ## name of Job
  #PBS -N jellyfish_Tityra

  ## Redirect output stream to this file.
  #PBS -o /media/inter/mkapun/projects/Tityra/results/denovo/log/GenomeSize_Tityra_log.txt

  ## Stream Standard Output AND Standard Error to outputfile (see above)
  #PBS -j oe

  ## Select a maximum walltime of 2h
  #PBS -l walltime=48:00:00

  ## Select 150 cores and 200gb of RAM
  #PBS -l select=1:ncpus=150:mem=1000g

  ## Go to pwd
  cd /media/inter/pipelines/AutDeNovo

  ######## load dependencies #######

  eval "$(conda shell.bash hook)"
  conda activate envs/jellyfish

  if [[ ( ILL == 'ILL' ) || ( ILL == 'ILL_ONT' ) || ( ILL == 'ILL_PB' ) || ( ILL == 'ILL_ONT_PB' ) ]]
  then

    ## unzip files
    gunzip -c /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1.fq.gz > /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1.fq

    if [[ -f /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2.fq.gz ]]; then
      gunzip -c /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2.fq.gz > /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2.fq
      wait
      jellyfish count       -C       -m 31       -s 100M       -t 150       -F 2       -o /media/inter/mkapun/projects/Tityra/results/denovo/results/GenomeSize/Tityra_reads.jf       /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1.fq       /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2.fq
      rm -f /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1.fq /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_2_val_2.fq
    else
      wait
      jellyfish count       -C       -m 31       -s 100M       -t 150       -F 2       -o /media/inter/mkapun/projects/Tityra/results/denovo/results/GenomeSize/Tityra_reads.jf       /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1.fq
      rm -f /media/inter/mkapun/projects/Tityra/results/denovo/data/Illumina/Tityra_1_val_1.fq
    fi
  fi

  if [[ ( ILL == 'ONT' ) || ( ILL == 'ONT_PB' ) ]]
  then

    ## unzip ONT file

    gunzip -c /media/inter/mkapun/projects/Tityra/results/denovo/data/ONT/.fq.gz > /media/inter/mkapun/projects/Tityra/results/denovo/data/ONT/.fq

    ## run Jellyfish

    ## parameters
    # -C canononical; count both strands
    # -m 31 Length of mers
    # -s initial hash size
    jellyfish count       -C       -m 31       -s 100M       -t 150       -F 2       -o /media/inter/mkapun/projects/Tityra/results/denovo/results/GenomeSize/Tityra_reads.jf       /media/inter/mkapun/projects/Tityra/results/denovo/data/ONT/.fq

    rm -f  /media/inter/mkapun/projects/Tityra/results/denovo/data/ONT/.fq

  fi

  if [[ ( ILL == 'PB' ) ]]
  then

    ## unzip file

    gunzip -c /media/inter/mkapun/projects/Tityra/results/denovo/data/PB/.fq.gz > /media/inter/mkapun/projects/Tityra/results/denovo/data/PB/.fq

    ## run Jellyfish

    ## parameters
    # -C canononical; count both strands
    # -m 31 Length of mers
    # -s initial hash size
    jellyfish count       -C       -m 31       -s 100M       -t 150       -F 2       -o /media/inter/mkapun/projects/Tityra/results/denovo/results/GenomeSize/Tityra_reads.jf       /media/inter/mkapun/projects/Tityra/results/denovo/data/PB/.fq

    rm -f  /media/inter/mkapun/projects/Tityra/results/denovo/data/PB/.fq

  fi

  jellyfish histo     -t 150     /media/inter/mkapun/projects/Tityra/results/denovo/results/GenomeSize/Tityra_reads.jf     > /media/inter/mkapun/projects/Tityra/results/denovo/results/GenomeSize/Tityra_reads.histo

  ## run GenomeScope
  conda deactivate
  conda activate envs/genomescope

  genomescope2   -i /media/inter/mkapun/projects/Tityra/results/denovo/results/GenomeSize/Tityra_reads.histo   -k 31   -p 2   -o /media/inter/mkapun/projects/Tityra/results/denovo/results/GenomeSize/Tityra

