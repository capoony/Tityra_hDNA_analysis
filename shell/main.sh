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

###############################################################################
# 5. mapDamage Analysis
###############################################################################
echo "Step 5: Running mapDamage analysis..."

# Download reference genome if not present
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

echo "Pipeline completed successfully."
