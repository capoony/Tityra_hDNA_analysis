#!/bin/bash

###############################################################################
# Tityra Analysis Pipeline
# Author: Martin Kapun
# Description: End-to-end pipeline for processing Tityra sequencing data.
###############################################################################

set -euo pipefail

# Set working directory
WD="/media/inter/mkapun/projects/Tityra_hDNA_analysis"

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

# Summarize results
python "${WD}/scripts/summarizeKraken.py" \
    --input "${WD}/results/kraken2/" \
    --output "${WD}/results/kraken2/kraken_summary.csv"

###############################################################################
# 5. mapDamage Analysis
###############################################################################
echo "Step 5: Running mapDamage analysis..."

REF_DIR="${WD}/data/ref"
REF_FASTA="GCA_013397135.1_ASM1339713v1_genomic.fna.gz"
REF_URL="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/397/135/GCA_013397135.1_ASM1339713v1/${REF_FASTA}"

mkdir -p "${REF_DIR}"
cd "${REF_DIR}"
if [ ! -f "${REF_FASTA}" ]; then
    wget -q "${REF_URL}"
fi

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

# Plot read depth sorted descending
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

# Prepare reference for mapDamage
gunzip -f "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz"
bgzip -f "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna"
samtools faidx "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz"

conda deactivate
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/mapdamage2

mapDamage -i "${WD}/results/minimap2/Tityra_leucura.bam" \
    -r "${REF_DIR}/${REF_FASTA}" \
    --rescale \
    --folder="${WD}/results/mapDamage/Tityra_leucura"

conda deactivate

# Convert PDFs to PNGs
for pdf in "${WD}/results/mapDamage/Tityra_leucura/"*.pdf; do
    png="${pdf%.pdf}.png"
    convert -density 300 "$pdf" -quality 90 "$png"
done

echo "Pipeline completed successfully!"

###############################################################################
# 6. Assemble Mitochondrial Genome
###############################################################################
echo "Step 6: Assembling mitochondrial genome..."
mkdir -p "${WD}/results/mitogenome"
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs

# Download mitochondrial genome
wget -q -O "${WD}/results/mitogenome/Pachyramphus_minor_mito.fasta" \
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_051035.1&db=nuccore&report=fasta&extrafeat=on&conwithfeat=on&hide-cdd=on&withparts=on"

# Prepare reference for BLAST
gunzip -c "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz" >"${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna"

makeblastdb -in "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna" \
    -dbtype nucl -out "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic"

blastn -query "${WD}/results/mitogenome/Pachyramphus_minor_mito.fasta" \
    -db "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -out "${WD}/results/mitogenome/mito_blast_results.txt" \
    -num_threads 200

# Map reads to mitochondrial genome
minimap2 -ax sr --secondary=no -t 200 \
    "${WD}/results/mitogenome/Pachyramphus_minor_mito.fasta" \
    "${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz" \
    "${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz" |
    samtools view -bS -F 4 - |
    samtools sort -o "${WD}/results/mitogenome/Tityra_leucura_PE_mito.bam"
samtools index "${WD}/results/mitogenome/Tityra_leucura_PE_mito.bam"

minimap2 -ax sr --secondary=no -t 200 \
    "${WD}/results/mitogenome/Pachyramphus_minor_mito.fasta" \
    "${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz" |
    samtools view -bS -F 4 - |
    samtools sort -o "${WD}/results/mitogenome/Tityra_leucura_merged_mito.bam"
samtools index "${WD}/results/mitogenome/Tityra_leucura_merged_mito.bam"

samtools merge -f "${WD}/results/mitogenome/Tityra_leucura_mito.bam" \
    "${WD}/results/mitogenome/Tityra_leucura_PE_mito.bam" \
    "${WD}/results/mitogenome/Tityra_leucura_merged_mito.bam"
samtools index "${WD}/results/mitogenome/Tityra_leucura_mito.bam"

samtools coverage --reference "${WD}/results/mitogenome/Pachyramphus_minor_mito.fasta" \
    "${WD}/results/mitogenome/Tityra_leucura_mito.bam" \
    >"${WD}/results/mitogenome/Tityra_leucura_mito.coverage.txt"

# Get consensus sequence
bgzip -f "${WD}/results/mitogenome/Pachyramphus_minor_mito.fasta"
samtools faidx "${WD}/results/mitogenome/Pachyramphus_minor_mito.fasta.gz"
samtools consensus -T "${WD}/results/mitogenome/Pachyramphus_minor_mito.fasta.gz" \
    "${WD}/results/mitogenome/Tityra_leucura_mito.bam" \
    >"${WD}/results/mitogenome/Tityra_leucura_mito_consensus.fasta"

bgzip -f "${WD}/results/mitogenome/Tityra_leucura_mito_consensus.fasta"
samtools faidx "${WD}/results/mitogenome/Tityra_leucura_mito_consensus.fasta.gz"

###############################################################################
# 7. Map trimmed reads to contaminant references and retain unmapped reads
###############################################################################
echo "Step 7: Mapping trimmed reads to contaminant references..."

declare -A contaminants=(
    ["Vanrija_pseudolonga"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/906/515/GCF_020906515.1_ASM2090651v1/GCF_020906515.1_ASM2090651v1_genomic.fna.gz"
    ["Penicillium_coprophilum"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/826/855/GCF_028826855.1_ASM2882685v1/GCF_028826855.1_ASM2882685v1_genomic.fna.gz"
    ["Homo_sapiens"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
    ["Malassezia_restricta"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/290/485/GCF_003290485.1_ASM329048v1/GCF_003290485.1_ASM329048v1_genomic.fna.gz"
    ["Aspergillus_cristatus"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/044/706/195/GCA_044706195.1_ASM4470619v1/GCA_044706195.1_ASM4470619v1_genomic.fna.gz"
)

mkdir -p "${WD}/results/contaminants"
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs

mkdir -p "${WD}/results/contaminants/mappings"
mkdir -p "${WD}/data/contaminants"

# Download contaminant references
for name in "${!contaminants[@]}"; do
    url="${contaminants[$name]}"
    ref_file="${WD}/data/contaminants/${name}.fna.gz"
    if [ ! -f "${ref_file}" ]; then
        echo "Downloading reference for ${name} from ${url}..."
        wget -q -O "${ref_file}" "${url}"
    fi
done

# Create joint reference
joint_ref="${WD}/data/contaminants/joint_reference.fna.gz"
if [ ! -f "${joint_ref}" ]; then
    echo "Creating joint reference file..."
    cat "${WD}/data/contaminants/"*.fna.gz >"${joint_ref}"
fi

# Map reads to joint reference and retain unmapped reads
echo "Mapping reads to the joint reference..."
minimap2 -ax sr --secondary=no -t 200 \
    "${joint_ref}" \
    "${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz" \
    "${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz" |
    samtools view -bS -f 4 - |
    samtools sort -n -o "${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_PE.bam"
samtools index "${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_PE.bam"

samtools fastq \
    -1 "${WD}/data/trimmed/Tityra_leucura_contaminants_PE_unmapped_1.fastq.gz" \
    -2 "${WD}/data/trimmed/Tityra_leucura_contaminants_PE_unmapped_2.fastq.gz" \
    -0 "${WD}/data/trimmed/Tityra_leucura_contaminants_PE_unmapped_single.fastq.gz" \
    -s "${WD}/data/trimmed/Tityra_leucura_contaminants_PE_unmapped_single2.fastq.gz" \
    "${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_PE.bam"

minimap2 -ax sr --secondary=no -t 200 \
    "${joint_ref}" \
    "${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz" |
    samtools view -bS -f 4 - >"${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_merged.bam"
samtools index "${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_merged.bam"

samtools fastq "${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_merged.bam" \
    >"${WD}/data/trimmed/Tityra_leucura_contaminants_merged_unmapped.fastq.gz"

###############################################################################
# 8. Run AutDeNovo Pipeline
# for documentation, see: https://github.com/nhmvienna/AutDeNovo
###############################################################################
/media/inter/pipelines/AutDeNovo/AutDeNovo_exp.sh \
    Name=Tityra \
    OutputFolder="${WD}/results/denovo" \
    Fwd="${WD}/data/trimmed/Tityra_leucura_contaminants_merged_unmapped.fastq.gz" \
    threads=150 \
    RAM=200 \
    RAMAssembly=1000 \
    decont=no \
    SmudgePlot=no \
    BLASTdb=/media/scratch/NCBI_nt_DB_210714/nt \
    BuscoDB=vertebrata_odb10 \
    Taxdump=/media/scratch/NCBI_taxdump/ \
    Racon=4

echo "All pipeline steps completed successfully!"
