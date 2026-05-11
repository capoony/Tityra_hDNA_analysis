
#!/bin/bash
#PBS -S /bin/bash
#PBS -N freebayes
#PBS -o /media/inter/mkapun/projects/Tityra/data/log_cluster.txt
#PBS -j oe
#PBS -l select=1:ncpus=100:mem=200gb

# Load Conda and FreeBayes environment
source /opt/anaconda3/etc/profile.d/conda.sh

module load NGSmapper/minimap2-2.17
module load Tools/samtools-1.18

minimap2 -ax sr --secondary=no -t 100     /media/inter/mkapun/projects/Tityra/data/contaminants/joint_reference.fna.gz     /media/inter/mkapun/projects/Tityra/data/trimmed2/Tityra_leucura_1_trimmed.fastq.gz     /media/inter/mkapun/projects/Tityra/data/trimmed2/Tityra_leucura_2_trimmed.fastq.gz     | samtools view -bS - | samtools sort >/media/inter/mkapun/projects/Tityra/results/contaminants/mappings/Tityra_leucura_contaminants.bam
samtools index /media/inter/mkapun/projects/Tityra/results/contaminants/mappings/Tityra_leucura_contaminants.bam


