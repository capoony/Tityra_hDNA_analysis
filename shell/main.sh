#!/bin/bash

###############################################################################
# Tityra Analysis Pipeline
# Author: [Your Name]
# Description: End-to-end pipeline for processing Tityra sequencing data.
###############################################################################

# Set working directory
WD="/media/inter/mkapun/projects/Tityra"

###############################################################################
# 1. Copy Raw Data
###############################################################################
echo "Step 1: Copying raw data..."
mkdir -p "${WD}/data/raw"
cp /media/inter/SeqData/raw/Macrogen/Illumina/Tityra_250430/data/*.gz "${WD}/data/raw"

###############################################################################
# 2. Trim Reads with fastp
###############################################################################
echo "Step 2: Trimming reads with fastp..."
mkdir -p "${WD}/data/trimmed"
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs

fastp \
    -i "${WD}/data/raw/Tityra_leucura_1.fastq.gz" \
    -I "${WD}/data/raw/Tityra_leucura_2.fastq.gz" \
    -o "${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz" \
    -O "${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz" \
    --merge \
    --merged_out "${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz" \
    --length_required 25 \
    --dedup \
    --trim_poly_g \
    --html "${WD}/data/trimmed/Tityra_leucura.html" \
    --json "${WD}/data/trimmed/Tityra_leucura.json" \
    --detect_adapter_for_pe

conda deactivate

###############################################################################
# 3. Run ECMSD Pipeline
###############################################################################
echo "Step 3: Running ECMSD pipeline..."
bash /media/inter/pipelines/ECMSD/shell/ECMSD.sh \
    --fwd "${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz" \
    --rev "${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz" \
    --merged "${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz" \
    --out "${WD}/results/ECMSD" \
    --threads 200 \
    --Binsize 1000 \
    --RMUS-threshold 0.15 \
    --mapping_quality 20 \
    --taxonomic-hierarchy genus \
    --force

###############################################################################
# 4. Kraken2 Taxonomic Classification
###############################################################################
echo "Step 4: Running Kraken2 analysis..."
mkdir -p "${WD}/results/kraken2"
module load Assembly/kraken-2.1.2

# Paired-end reads
kraken2 \
    --db /media/scratch/kraken-2.1.2/db/pluspfp_20240904 \
    --threads 150 \
    --report "${WD}/results/kraken2/report_PE.txt" \
    --paired \
    "${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz" \
    "${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz" \
    >/dev/null

# Merged reads
kraken2 \
    --db /media/scratch/kraken-2.1.2/db/pluspfp_20240904 \
    --threads 150 \
    --report "${WD}/results/kraken2/report_merged.txt" \
    "${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz" \
    >/dev/null

# summarize results
python ${WD}/scripts/summarizeKraken.py \
    --input ${WD}/results/kraken2/ \
    --output ${WD}/results/kraken2/kraken_summary.csv

###############################################################################
# 5. mapDamage Analysis
###############################################################################
echo "Step 5: Running mapDamage analysis..."

# Download reference genome of Pachyramphus minor
REF_DIR="${WD}/data/ref"
REF_FASTA="GCA_013397135.1_ASM1339713v1_genomic.fna.gz"
REF_URL="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/397/135/GCA_013397135.1_ASM1339713v1/${REF_FASTA}"

mkdir -p "${REF_DIR}"
cd "${REF_DIR}"
if [ ! -f "${REF_FASTA}" ]; then
    wget -q "${REF_URL}"
fi

# Align reads with minimap2 and process with samtools
mkdir -p "${WD}/results/minimap2"
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs

# Align paired-end reads
minimap2 -ax sr --secondary=no -t 200 \
    "${REF_DIR}/${REF_FASTA}" \
    "${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz" \
    "${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz" |
    samtools view -bS -F 4 - |
    samtools sort -o "${WD}/results/minimap2/Tityra_leucura_PE.bam"
samtools index "${WD}/results/minimap2/Tityra_leucura_PE.bam"

# Align merged reads
minimap2 -ax sr --secondary=no -t 200 \
    "${REF_DIR}/${REF_FASTA}" \
    "${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz" |
    samtools view -bS -F 4 - |
    samtools sort -o "${WD}/results/minimap2/Tityra_leucura_merged.bam"
samtools index "${WD}/results/minimap2/Tityra_leucura_merged.bam"

# Merge BAM files
samtools merge -f "${WD}/results/minimap2/Tityra_leucura.bam" \
    "${WD}/results/minimap2/Tityra_leucura_PE.bam" \
    "${WD}/results/minimap2/Tityra_leucura_merged.bam"
samtools index "${WD}/results/minimap2/Tityra_leucura.bam"

samtools coverage --reference "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz" \
    "${WD}/results/minimap2/Tityra_leucura.bam" \
    >"${WD}/results/minimap2/Tityra_leucura.coverage.txt"

## Plot read depth sorted  descending

Rscript -e "
library(tidyverse)
## header starts With #
coverage_data <- read.table('${WD}/results/minimap2/Tityra_leucura.coverage.txt', header=TRUE, comment.char='-')
coverage_data <- coverage_data %>%
    arrange(desc(meandepth)) 

# plot without x-axis labels and add median depth as red vertical line and write median depth as text on top and restrict to first 1000 longest contigs based on covbases
coverage_data<- coverage_data %>%
    arrange(desc(covbases))

coverage_data.new<- coverage_data[1:1000, ]
ggplot(coverage_data.new, aes(x = reorder(rname, -meandepth), y = meandepth)) +
    geom_bar(stat='identity', fill='steelblue') +
    labs(x='Contig', y='Mean Depth', title='Median Read Depth for Tityra leucura for longest 1000 contigs') +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept=median(coverage_data.new\$meandepth), color='red', linetype='dashed') +
    geom_text(aes(x=100, y=median(coverage_data.new\$meandepth) + 0.5, label=paste('Median Depth:', round(median(coverage_data.new\$meandepth), 2))), color='red', size=4) +
    theme_bw()
ggsave('${WD}/results/minimap2/Tityra_leucura.coverage_plot.png', width=10, height=6)
"

# Run mapDamage

gunzip ${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz
bgzip -f "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna"
samtools faidx "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz"

conda deactivate
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/mapdamage2
mapDamage -i "${WD}/results/minimap2/Tityra_leucura.bam" \
    -r "${REF_DIR}/${REF_FASTA}" \
    --rescale \
    --folder="${WD}/results/mapDamage/Tityra_leucura"
conda deactivate

## convert PDFs to PNGs (only works on Linux systems)
for pdf in ${WD}/results/mapDamage/Tityra_leucura/*.pdf; do
    png=${pdf%.pdf}.png
    convert -density 300 $pdf -quality 90 $png
done

echo "Pipeline completed successfully."
