#!/bin/bash

###############################################################################
# Parallel variant calling for all graph-based genes
# Usage: bash parallel_variant_calling.sh <num_jobs>
###############################################################################

set -euo pipefail

WD=/media/inter/mkapun/projects/Tityra
VG=/opt/bioinformatics/cactus-2.1.1/bin/vg
JOBS=${1:-10}  # Number of parallel jobs
THREADS_PER_JOB=4  # Threads per vg map job

echo "=================================================================="
echo "Variant calling for all graph-based BUSCO genes (parallel)"
echo "Parallel jobs: ${JOBS}"
echo "Threads per job: ${THREADS_PER_JOB}"
echo "=================================================================="

mkdir -p ${WD}/results/graph_alignment/mapped
mkdir -p ${WD}/results/graph_alignment/consensus

READS_PE1="${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz"
READS_PE2="${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz"
READS_MERGED="${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz"

export WD VG THREADS_PER_JOB READS_PE1 READS_PE2 READS_MERGED

# Function to process one gene
process_gene() {
    local graph_xg=$1
    local gene=$(basename ${graph_xg} .xg)
    
    # Skip if consensus already exists
    if [ -f "${WD}/results/graph_alignment/consensus/${gene}_Tityra.fasta" ]; then
        return 0
    fi
    
    # Map paired-end reads
    ${VG} map \
        -x ${graph_xg} \
        -g ${WD}/results/graph_alignment/graphs/${gene}.gcsa \
        -f ${READS_PE1} -f ${READS_PE2} \
        -t ${THREADS_PER_JOB} \
        > ${WD}/results/graph_alignment/mapped/${gene}.gam 2>/dev/null || return 1
    
    # Map merged reads (CRITICAL!)
    ${VG} map \
        -x ${graph_xg} \
        -g ${WD}/results/graph_alignment/graphs/${gene}.gcsa \
        -f ${READS_MERGED} \
        -t ${THREADS_PER_JOB} \
        >> ${WD}/results/graph_alignment/mapped/${gene}.gam 2>/dev/null || return 1
    
    # Pack alignments
    ${VG} pack \
        -x ${graph_xg} \
        -g ${WD}/results/graph_alignment/mapped/${gene}.gam \
        -o ${WD}/results/graph_alignment/mapped/${gene}.pack \
        2>/dev/null || return 1
    
    # Call variants with relaxed filters for ancient DNA
    ${VG} call \
        ${graph_xg} \
        -k ${WD}/results/graph_alignment/mapped/${gene}.pack \
        -d 1 \
        -s 1 \
        > ${WD}/results/graph_alignment/mapped/${gene}_raw.vcf \
        2>/dev/null || return 1
    
    # Get first reference
    first_ref=$(grep "^>" ${WD}/results/graph_alignment/busco_regions/${gene}.fasta | head -n 1 | sed 's/^>//')
    
    # Extract reference sequence
    awk -v ref="^>${first_ref}" 'BEGIN{p=0} $0 ~ ref {p=1; print; next} /^>/ {p=0} p' \
        ${WD}/results/graph_alignment/busco_regions/${gene}.fasta \
        > ${WD}/results/graph_alignment/mapped/${gene}_ref.fasta
    
    # Filter variants: only keep those with depth > 5
    bcftools view \
        -i 'FORMAT/DP>5' \
        ${WD}/results/graph_alignment/mapped/${gene}_raw.vcf \
        -o ${WD}/results/graph_alignment/mapped/${gene}_filtered.vcf \
        2>/dev/null || return 1
    
    # Compress and index filtered VCF
    bgzip -c ${WD}/results/graph_alignment/mapped/${gene}_filtered.vcf > ${WD}/results/graph_alignment/mapped/${gene}.vcf.gz 2>/dev/null
    tabix -p vcf ${WD}/results/graph_alignment/mapped/${gene}.vcf.gz 2>/dev/null
    samtools faidx ${WD}/results/graph_alignment/mapped/${gene}_ref.fasta 2>/dev/null
    
    # Create consensus
    bcftools consensus \
        -f ${WD}/results/graph_alignment/mapped/${gene}_ref.fasta \
        ${WD}/results/graph_alignment/mapped/${gene}.vcf.gz \
        > ${WD}/results/graph_alignment/consensus/${gene}_Tityra.fasta \
        2>/dev/null || return 1
    
    # Rename sequence
    sed -i "1s/.*/>${gene}|Tityra_leucura_graph/" ${WD}/results/graph_alignment/consensus/${gene}_Tityra.fasta
    
    # Clean up large files
    rm -f ${WD}/results/graph_alignment/mapped/${gene}.gam
    rm -f ${WD}/results/graph_alignment/mapped/${gene}_ref.fasta
    
    return 0
}

export -f process_gene

# Process all genes in parallel
find ${WD}/results/graph_alignment/graphs/ -name "*.xg" | \
    parallel -j ${JOBS} --progress process_gene {}

echo ""
echo "=================================================================="
echo "Variant calling complete!"
echo "Consensus sequences: $(ls ${WD}/results/graph_alignment/consensus/*_Tityra.fasta 2>/dev/null | wc -l)"
echo "=================================================================="
