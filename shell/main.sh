#!/bin/bash

###############################################################################
# Tityra Analysis Pipeline
# Author: Martin Kapun
# Description: End-to-end pipeline for processing Tityra sequencing data.
###############################################################################

set -euo pipefail

# Set working directory
WD=/media/inter/mkapun/projects/Tityra

###############################################################################
# 1. Copy Raw Data
###############################################################################
echo "Step 1: Copying raw data..."
mkdir -p "${WD}/data/raw" # Create directory for raw data if it doesn't exist
cp /media/inter/SeqData/raw/Macrogen/Illumina/Tityra_250430/data/*.gz "${WD}/data/raw" # Copy raw .gz files

###############################################################################
# 2. Trim Reads with fastp
###############################################################################
echo "Step 2: Trimming reads with fastp..."
mkdir -p "${WD}/data/trimmed" # Create directory for trimmed reads
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs # Activate conda environment with fastp

# Run fastp for paired-end trimming, merging, deduplication, and QC reporting
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

mkdir -p "${WD}/data/trimmed2" # Create directory for trimmed reads if it doesn't exist
# Run fastp for paired-end trimming, merging, deduplication, and QC reporting and additionally trim away the first three bases of each read at the 5'end 
fastp \
    -i "${WD}/data/raw/Tityra_leucura_1.fastq.gz" \
    -I "${WD}/data/raw/Tityra_leucura_2.fastq.gz" \
    -o "${WD}/data/trimmed2/Tityra_leucura_1_trimmed.fastq.gz" \
    -O "${WD}/data/trimmed2/Tityra_leucura_2_trimmed.fastq.gz" \
    --merge \
    --merged_out "${WD}/data/trimmed2/Tityra_leucura_merged.fastq.gz" \
    --length_required 25 \
    --cut_front --cut_front_window_size 3 --cut_front_mean_quality 20 \
    --dedup \
    --trim_poly_g \
    --html "${WD}/data/trimmed2/Tityra_leucura.html" \
    --json "${WD}/data/trimmed2/Tityra_leucura.json" \
    --detect_adapter_for_pe

conda deactivate # Deactivate conda environment

###############################################################################
# 3. Run ECMSD Pipeline
###############################################################################
echo "Step 3: Running ECMSD pipeline..."
# Run ECMSD pipeline for metagenomic analysis
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
mkdir -p "${WD}/results/kraken2" # Create output directory for Kraken2
module load Assembly/kraken-2.1.2 # Load Kraken2 module

# Run Kraken2 on paired-end reads for taxonomic classification
kraken2 \
    --db /media/scratch/kraken-2.1.2/db/pluspfp_20240904 \
    --threads 150 \
    --report "${WD}/results/kraken2/report_PE.txt" \
    --paired \
    "${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz" \
    "${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz" \
    >/dev/null

# Run Kraken2 on merged reads
kraken2 \
    --db /media/scratch/kraken-2.1.2/db/pluspfp_20240904 \
    --threads 150 \
    --report "${WD}/results/kraken2/report_merged.txt" \
    "${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz" \
    >/dev/null

# Summarize Kraken2 results using custom Python script
python "${WD}/scripts/summarizeKraken.py" \
    --input "${WD}/results/kraken2/" \
    --output "${WD}/results/kraken2/kraken_summary.csv"

###############################################################################
# 5. mapDamage Analysis
###############################################################################
echo "Step 5: Running mapDamage analysis..."

# Set reference genome variables
REF_DIR="${WD}/data/ref"
REF_FASTA="GCA_013397135.1_ASM1339713v1_genomic.fna.gz"
REF_URL="ftp:/ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/397/135/GCA_013397135.1_ASM1339713v1/${REF_FASTA}"

mkdir -p "${REF_DIR}" # Create reference directory
cd "${REF_DIR}" # Change to reference directory
if [ ! -f "${REF_FASTA}" ]; then
    wget -q "${REF_URL}" # Download reference genome if not present
fi

mkdir -p "${WD}/results/minimap2" # Create minimap2 output directory
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs # Activate conda env

# Align paired-end reads to reference genome with minimap2, filter, and sort
minimap2 -ax sr --secondary=no -t 200 \
    "${REF_DIR}/${REF_FASTA}" \
    "${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz" \
    "${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz" |
    samtools view -bS -F 4 - |
    samtools sort -o "${WD}/results/minimap2/Tityra_leucura_PE.bam"
samtools index "${WD}/results/minimap2/Tityra_leucura_PE.bam" # Index BAM

# Align merged reads to reference genome
minimap2 -ax sr --secondary=no -t 200 \
    "${REF_DIR}/${REF_FASTA}" \
    "${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz" |
    samtools view -bS -F 4 - |
    samtools sort -o "${WD}/results/minimap2/Tityra_leucura_merged.bam"
samtools index "${WD}/results/minimap2/Tityra_leucura_merged.bam"

# Merge paired-end and merged BAM files
samtools merge -f "${WD}/results/minimap2/Tityra_leucura.bam" \
    "${WD}/results/minimap2/Tityra_leucura_PE.bam" \
    "${WD}/results/minimap2/Tityra_leucura_merged.bam"
samtools index "${WD}/results/minimap2/Tityra_leucura.bam"

# Calculate coverage statistics
samtools coverage --reference "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz" \
    "${WD}/results/minimap2/Tityra_leucura.bam" \
    >"${WD}/results/minimap2/Tityra_leucura.coverage.txt"

# Plot read depth for the 1000 longest contigs using R
Rscript -e "
library(tidyverse)
coverage_data <- read.table('${WD}/results/minimap2/Tityra_leucura.coverage.txt', header=TRUE, comment.char='-')
coverage_data <- coverage_data %>% arrange(desc(meandepth), desc(covbases))
coverage_data.new <- coverage_data[1:1000, ]
ggplot(coverage_data.new, aes(x = reorder(rname, -meandepth), y = meandepth)) +
    geom_bar(stat='identity', fill='steelblue') +
    labs(x='Contig', y='Mean Depth', title='Median Read Depth for Tityra leucura for longest 1000 contigs') +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept=median(coverage_data.new\$meandepth), color='red', linetype='dashed') +
    geom_text(aes(x=100, y=median(coverage_data.new\$meandepth) + 0.5, label=paste('Median Depth:', round(median(coverage_data.new\$meandepth), 2))), color='red', size=4) +
    theme_bw()
ggsave('${WD}/results/minimap2/Tityra_leucura.coverage_plot.png', width=10, height=6)
"

# Prepare reference genome for mapDamage (bgzip and index)
gunzip -f "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz"
bgzip -f "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna"
samtools faidx "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz"

conda deactivate
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/mapdamage2 # Activate mapDamage env

# Run mapDamage to assess DNA damage patterns
mapDamage -i "${WD}/results/minimap2/Tityra_leucura.bam" \
    -r "${REF_DIR}/${REF_FASTA}" \
    --rescale \
    --folder="${WD}/results/mapDamage/Tityra_leucura"

conda deactivate

# Convert mapDamage PDF plots to PNG images
for pdf in "${WD}/results/mapDamage/Tityra_leucura/"*.pdf; do
    png="${pdf%.pdf}.png"
    convert -density 300 "$pdf" -quality 90 "$png"
done

echo "Pipeline completed successfully!"

###############################################################################
# 6. Map trimmed reads to contaminant references and retain unmapped reads
###############################################################################
echo "Step 7: Mapping trimmed reads to contaminant references..."

# Declare associative array of contaminant references and URLs
declare -A contaminants=(
    ["Vanrija_pseudolonga"]="https:/ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/906/515/GCF_020906515.1_ASM2090651v1/GCF_020906515.1_ASM2090651v1_genomic.fna.gz"
    ["Penicillium_coprophilum"]="https:/ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/826/855/GCF_028826855.1_ASM2882685v1/GCF_028826855.1_ASM2882685v1_genomic.fna.gz"
    ["Homo_sapiens"]="https:/ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
    ["Malassezia_restricta"]="https:/ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/290/485/GCF_003290485.1_ASM329048v1/GCF_003290485.1_ASM329048v1_genomic.fna.gz"
    ["Aspergillus_cristatus"]="https:/ftp.ncbi.nlm.nih.gov/genomes/all/GCA/044/706/195/GCA_044706195.1_ASM4470619v1/GCA_044706195.1_ASM4470619v1_genomic.fna.gz"
)

mkdir -p "${WD}/results/contaminants"
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs

mkdir -p "${WD}/results/contaminants/mappings"
mkdir -p "${WD}/data/contaminants"

# Download contaminant reference genomes if not present
for name in "${!contaminants[@]}"; do
    url="${contaminants[$name]}"
    ref_file="${WD}/data/contaminants/${name}.fna.gz"
    if [ ! -f "${ref_file}" ]; then
        echo "Downloading reference for ${name} from ${url}..."
        wget -q -O "${ref_file}" "${url}"
    fi
done

# Create joint reference file by concatenating all contaminant references
joint_ref="${WD}/data/contaminants/joint_reference.fna.gz"
if [ ! -f "${joint_ref}" ]; then
    echo "Creating joint reference file..."
    cat "${WD}/data/contaminants/"*.fna.gz >"${joint_ref}"
fi

# Map merged reads to joint contaminant reference and retain unmapped reads
minimap2 -ax sr --secondary=no -t 200 \
    "${joint_ref}" \
    "${WD}/data/trimmed2/Tityra_leucura_merged.fastq.gz" | awk '($5 < 20 || and($2,4)) || $1 ~ /^@/' | samtools view -bS - >"${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_merged.bam"
samtools index "${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_merged.bam"

# Extract unmapped merged reads
samtools fastq "${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_merged.bam" \
    | pigz >"${WD}/data/trimmed2/Tityra_leucura_contaminants_merged_unmapped.fastq.gz"

###############################################################################
# 7. Assemble Mitochondrial Genome
###############################################################################
echo "Step 6: Assembling mitochondrial genome..."
mkdir -p "${WD}/results/mitogenome" # Create output directory

sh ${WD}/shell/mitofinder_assembly.sh


###############################################################################
# 8. Get more COI and other mitochondrial sequences from other Tyrranidae
###############################################################################

mkdir ${WD}/data/Mito
for j in COI CO2 ND1 ND2 ND3 ND4 ND4L ND5 N6 CYTB ATP8 COX1 COX2 COX3 ATP6 rrnL rrnS; do

    for i in "Onychorhynchus coronatus" Tityra Schiffornis Laniocera Iodopleura Laniisoma Xenopsaris Pachyramphus; do
        
        python /media/inter/mkapun/projects/EfetchThePython/EfetchThePython.py \
            --email capopony@gmail.com \
            --api_key dfbe3972fc58b8478eba0929f14847f9cd08 \
            --Term "\"$i\"[Organism], $j" \
            --Output ${WD}/data/Mito/${j}_$i \
            --FASTA
    done
done

###############################################################################
# 7. Combine mitochondrial gene sequences
###############################################################################
# Combine downloaded gene sequences with Tityra sequences from MitoFinder

echo ""
echo "==================================================================="
echo "Combining mitochondrial gene sequences (assembly-based)..."
echo "==================================================================="

bash ${WD}/shell/combine_mito_genes_assembly.sh

###############################################################################
# 8. Phylogenetic analysis of species-level combined genes
###############################################################################

echo ""
echo "==================================================================="
echo "Running phylogenetic analysis (assembly-based)..."
echo "==================================================================="

bash ${WD}/shell/phylogenetic_analysis_assembly.sh

###############################################################################
# 8C. Concatenate genes by species (Assembly-based)
###############################################################################
echo ""
echo "==================================================================="
echo "Concatenating mitochondrial genes by species (assembly-based)..."
echo "==================================================================="

bash ${WD}/shell/concatenate_genes_assembly.sh

###############################################################################
# 10. Phylogenetic analysis of concatenated genes (Assembly-based)
###############################################################################
echo ""
echo "==================================================================="
echo "Running phylogenetic analysis on concatenated genes (assembly-based)..."
echo "==================================================================="

bash ${WD}/shell/phylogenetic_analysis_concatenated_assembly.sh
