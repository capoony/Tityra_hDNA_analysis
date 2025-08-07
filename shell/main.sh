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
# 6. Assemble Mitochondrial Genome
###############################################################################
echo "Step 6: Assembling mitochondrial genome..."
mkdir -p "${WD}/results/mitogenome" # Create output directory
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs

# Download mitochondrial genome reference
wget -q -O "${WD}/results/mitogenome/Pachyramphus_minor_mito.fasta" \
    "https:/www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_051035.1&db=nuccore&report=fasta&extrafeat=on&conwithfeat=on&hide-cdd=on&withparts=on"

# Prepare reference for BLAST by decompressing and indexing
gunzip -c "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz" >"${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna"

# Create BLAST database from reference genome
makeblastdb -in "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna" \
    -dbtype nucl -out "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic"

# BLAST mitochondrial genome against reference to find location
blastn -query "${WD}/results/mitogenome/Pachyramphus_minor_mito.fasta" \
    -db "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -out "${WD}/results/mitogenome/mito_blast_results.txt" \
    -num_threads 200

# Map reads to mitochondrial genome reference
minimap2 -ax sr --secondary=no -t 200 \
    "${WD}/results/mitogenome/Pachyramphus_minor_mito.fasta" \
    "${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz" \
    "${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz" |
    samtools view -bS -F 4 - |
    samtools sort -o "${WD}/results/mitogenome/Tityra_leucura_PE_mito.bam"
samtools index "${WD}/results/mitogenome/Tityra_leucura_PE_mito.bam"

# Map merged reads to mitochondrial genome
minimap2 -ax sr --secondary=no -t 200 \
    "${WD}/results/mitogenome/Pachyramphus_minor_mito.fasta" \
    "${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz" |
    samtools view -bS -F 4 - |
    samtools sort -o "${WD}/results/mitogenome/Tityra_leucura_merged_mito.bam"
samtools index "${WD}/results/mitogenome/Tityra_leucura_merged_mito.bam"

# Merge BAM files for mitochondrial mapping
samtools merge -f "${WD}/results/mitogenome/Tityra_leucura_mito.bam" \
    "${WD}/results/mitogenome/Tityra_leucura_PE_mito.bam" \
    "${WD}/results/mitogenome/Tityra_leucura_merged_mito.bam"
samtools index "${WD}/results/mitogenome/Tityra_leucura_mito.bam"

# Calculate coverage for mitochondrial genome
samtools coverage --reference "${WD}/results/mitogenome/Pachyramphus_minor_mito.fasta" \
    "${WD}/results/mitogenome/Tityra_leucura_mito.bam" \
    >"${WD}/results/mitogenome/Tityra_leucura_mito.coverage.txt"

# Generate consensus mitochondrial sequence
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

# Map paired-end reads to joint contaminant reference and retain unmapped reads
echo "Mapping reads to the joint reference..."
minimap2 -ax sr --secondary=no -t 200 \
    "${joint_ref}" \
    "${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz" \
    "${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz" |
    samtools view -bS -f 4 - |
    samtools sort -n -o "${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_PE.bam"
samtools index "${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_PE.bam"

# Extract unmapped paired and single reads to separate fastq files
samtools fastq \
    -1 "${WD}/data/trimmed/Tityra_leucura_contaminants_PE_unmapped_1.fastq.gz" \
    -2 "${WD}/data/trimmed/Tityra_leucura_contaminants_PE_unmapped_2.fastq.gz" \
    -0 "${WD}/data/trimmed/Tityra_leucura_contaminants_PE_unmapped_single.fastq.gz" \
    -s "${WD}/data/trimmed/Tityra_leucura_contaminants_PE_unmapped_single2.fastq.gz" \
    "${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_PE.bam"

# Map merged reads to joint contaminant reference and retain unmapped reads
minimap2 -ax sr --secondary=no -t 200 \
    "${joint_ref}" \
    "${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz" |
    samtools view -bS -f 4 - >"${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_merged.bam"
samtools index "${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_merged.bam"

# Extract unmapped merged reads
samtools fastq "${WD}/results/contaminants/mappings/Tityra_leucura_contaminants_merged.bam" \
    >"${WD}/data/trimmed/Tityra_leucura_contaminants_merged_unmapped.fastq.gz"

###############################################################################
# 8. Run AutDeNovo Pipeline
# for documentation, see: https:/github.com/nhmvienna/AutDeNovo
###############################################################################
# Run AutDeNovo pipeline for de novo assembly and annotation
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

###############################################################################
# 9. Make Phylogeny with BUSCO genes
###############################################################################

echo "Step 9: Making phylogeny with BUSCO genes..."

mkdir -p ${WD}/results/phylogeny/data # Create directory for phylogeny data

## download reference genomes (commented out, see genomes.txt for URLs)
while IFS=$"," read -r Name ID URL; do
    #curl -L -o "${WD}/results/phylogeny/data/${Name}.fna.gz" "${URL}"
    echo ${Name} >>${WD}/results/phylogeny/data/genomes.names
done <${WD}/results/phylogeny/data/genomes.txt

## copy Pachyrhampus genome to phylogeny data directory
cp ${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz \
    ${WD}/results/phylogeny/data/Pachyramphus_minor.fna.gz
echo "Pachyramphus_minor" >>${WD}/results/phylogeny/data/genomes.names

## copy Tityra genome (de novo assembly) to phylogeny data directory
cp ${WD}/results/denovo/output/Tityra_ILL.fa.gz \
    ${WD}/results/phylogeny/data/Tityra_leucura.fna.gz
echo "Tityra_leucura" >>${WD}/results/phylogeny/data/genomes.names

mkdir -p ${WD}/results/phylogeny/BUSCO # Create BUSCO output directory

# Run BUSCO for each genome to identify single-copy orthologs
while IFS=$"," read -r Name; do

    echo """

  #!/bin/sh

  ## name of Job
  #PBS -N BUSCO_${Name}

  ## Redirect output stream to this file.
  #PBS -o ${WD}/results/phylogeny/BUSCO/${Name}_log.txt

  ## Stream Standard Output AND Standard Error to outputfile (see above)
  #PBS -j oe

  ## Select a maximum of 200 cores and 1000gb of RAM
  #PBS -l select=1:ncpus=50:mem=500g

  ######## load dependencies #######

  source /opt/anaconda3/etc/profile.d/conda.sh
  conda activate busco_5.4.3

  ######## run analyses #######

  cd ${WD}/results/phylogeny/BUSCO

  pigz -d ${WD}/results/phylogeny/data/${Name}.fna.gz

  busco -i ${WD}/results/phylogeny/data/${Name}.fna \
    -o ${Name} \
    -m genome \
    -c 200 \
    -f \
    -l aves_odb10

  pigz -f ${WD}/results/phylogeny/data/${Name}.fna

  """ >${WD}/results/phylogeny/BUSCO/${Name}.sh

    sh ${WD}/results/phylogeny/BUSCO/${Name}.sh

done <${WD}/results/phylogeny/data/genomes.names

## concatenate BUSCO results
mkdir -p ${WD}/results/phylogeny/concatenated

cd ${WD}/results/phylogeny/BUSCO

## identify all "complete" BUSCO genes across all genomes and concatenate their IDs
for file in $(find -iname "full_table*.tsv"); do
    grep -v "^#" ${file} | awk '$2=="Complete" {print $1}' >>${WD}/results/phylogeny/concatenated/complete_busco_ids.txt
done

## sort all BUSCO IDs and count occurrences
sort ${WD}/results/phylogeny/concatenated/complete_busco_ids.txt |
    uniq -c \
        >${WD}/results/phylogeny/concatenated/complete_busco_ids_with_counts.txt

## filter for BUSCO genes that are present in all 14 genomes
awk '$1 == 14 {print $2}' ${WD}/results/phylogeny/concatenated/complete_busco_ids_with_counts.txt \
    >${WD}/results/phylogeny/concatenated/final_busco_ids.txt

mkdir -p ${WD}/results/phylogeny/concatenated/busco_aa

# Copy and rename BUSCO protein sequences for each genome and gene
while IFS=$"," read -r Name; do

    while IFS=$"," read -r gene; do

        cp ${WD}/results/phylogeny/BUSCO/${Name}/run_aves_odb10/busco_sequences/single_copy_busco_sequences/${gene}.faa \
            ${WD}/results/phylogeny/concatenated/busco_aa/${Name}_${gene}
        sed -i 's/^>/>'${Name}'|/g' ${WD}/results/phylogeny/concatenated/busco_aa/${Name}_${gene}
    done <${WD}/results/phylogeny/concatenated/final_busco_ids.txt

done <${WD}/results/phylogeny/data/genomes.names

### concatenate and reduce to shared genes

mkdir -p ${WD}/results/phylogeny/prealigned

# Concatenate protein sequences for each BUSCO gene
while read gene; do
    cat ${WD}/results/phylogeny/concatenated/busco_aa/*_${gene} \
        >>${WD}/results/phylogeny/prealigned/${gene}_aa.fasta
done <${WD}/results/phylogeny/concatenated/final_busco_ids.txt

mkdir ${WD}/results/phylogeny/mafft

conda activate mafft-7.487
# Align each BUSCO gene with MAFFT
for i in ${WD}/results/phylogeny/prealigned/*_aa.fasta; do

    TMP=${i##*/}
    ID=${TMP%_*}

    mafft \
        --thread 50 \
        --auto \
        ${i} \
        >${WD}/results/phylogeny/mafft/${ID}_aln.fasta

done

# Fix sequence IDs after MAFFT alignment
for i in ${WD}/results/phylogeny/prealigned/*_aa.fasta; do

    TMP=${i##*/}
    ID=${TMP%_*}

    python /media/inter/mkapun/projects/TorpedoDeNovo/scripts/fixIDAfterMafft.py \
        --Alignment ${WD}/results/phylogeny/mafft/${ID}_aln.fasta \
        --input ${i} \
        >${WD}/results/phylogeny/mafft/${ID}_aln_fixed.fasta

    rm -f ${WD}/results/phylogeny/mafft/${ID}_aln.fasta

done

## make phylogeny

mkdir ${WD}/results/phylogeny/phylogeny

# Concatenate all alignments for phylogenetic analysis
python /media/inter/mkapun/projects/TorpedoDeNovo/scripts/proteins2genome.py \
    --input ${WD}/results/phylogeny/mafft \
    >${WD}/results/phylogeny/phylogeny/alignment.fa

module load Phylogeny/RAxML-2.8.10

## make new directory

cd ${WD}/results/phylogeny/phylogeny

## run ML tree reconstruction with RAxML
raxmlHPC-PTHREADS-SSE3 \
    -m PROTGAMMAWAG \
    -N 20 \
    -p 772374015 \
    -n Tityra \
    -s ${WD}/results/phylogeny/phylogeny/alignment.fa \
    -o Acanthisitta_chloris \
    -T 200

# Run RAxML bootstrapping
raxmlHPC-PTHREADS-SSE3 \
    -m PROTGAMMAWAG \
    -N autoMRE \
    -p 772374015 \
    -b 444353738 \
    -n bootrep \
    -s ${WD}/results/phylogeny/phylogeny/alignment.fa \
    -o Acanthisitta_chloris \
    -T 200

# Reconcile best ML tree with bootstrap replicates
raxmlHPC-SSE3 -f b \
    -m GTRGAMMA \
    -t RAxML_bestTree.Tityra \
    -z RAxML_bootstrap.bootrep \
    -n FINAL \
    -o Acanthisitta_chloris

# Plot phylogenetic tree using R and ggtree
Rscript -e """
# load necessary R libraries
library('ggtree')
library('gridExtra')
library('ggrepel')
library('ape')
library('ggplot2')
library('phangorn')
library('dplyr')
library(phytools) # to determine the maximum tree height and add midpoint root

## load tree file and root with outgroup taxa
tree<-read.tree('${WD}/results/phylogeny/phylogeny/RAxML_bipartitions.FINAL')


## plot tree

## caluculate tree height (on x-axis)
Xmax<-max(nodeHeights(tree))

tree<-root(tree,outgroup='Acanthisitta_chloris')

## plot tree
PLOT.tree<-ggtree(tree)+
  ggtitle('Tityra_BUSCO')+
  theme_tree2()+
  theme_bw()+
  ggplot2::xlim(0,
    Xmax+0.25)+
  xlab('av. subst./site') +
  geom_nodelab()+
  theme(axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  geom_tiplab()

## export tree
ggsave(filename='${WD}/results/phylogeny/phylogeny/Tityra_BUSCO.pdf',
  PLOT.tree)
ggsave(filename='${WD}/results/phylogeny/phylogeny/Tityra_BUSCO.png',
  PLOT.tree)
"""

## OK, now repeat with mapping reads to Pachyrhamphus minor BUSCO genes
sh ${WD}/shell/BUSCO_Pminor.sh ${WD}

## OK, now repeat with mapping Tityra reads to Ochrsit BUSCO genes
sh ${WD}/shell/BUSCO_Ocrist.sh ${WD}

## OK, now repeat with mapping Tityra reads to Ochrsit BUSCO genes
sh ${WD}/shell/BUSCO_tsav.sh ${WD}

## Make Phylogenies for all genes based on the denovo assembly 
sh ${WD}/shell/BUSCO_phylogeny_denovo_indgenes.sh ${WD}
