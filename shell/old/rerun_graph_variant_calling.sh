#!/bin/bash

###############################################################################
# Re-run variant calling for all graph-based genes with corrected parameters
# Fixes: 
# 1. Ensure merged reads are mapped
# 2. Use relaxed depth filter (-d 1 -s 1)
# 3. Set VG path explicitly
###############################################################################

set -euo pipefail

WD=/media/inter/mkapun/projects/Tityra
VG=/opt/bioinformatics/cactus-2.1.1/bin/vg
THREADS=4

echo "=================================================================="
echo "Re-running variant calling for all 20 graph-based BUSCO genes"
echo "Working directory: ${WD}"
echo "VG binary: ${VG}"
echo "=================================================================="

# Clear all existing results to force complete re-run
echo "Clearing old consensus and variant files..."
rm -f ${WD}/results/graph_alignment/consensus/*_Tityra.fasta
rm -f ${WD}/results/graph_alignment/mapped/*.gam
rm -f ${WD}/results/graph_alignment/mapped/*.pack
rm -f ${WD}/results/graph_alignment/mapped/*_raw.vcf*
rm -f ${WD}/results/graph_alignment/mapped/*.vcf.gz*

# Read files
READS_PE1="${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz"
READS_PE2="${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz"
READS_MERGED="${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz"

# Process each gene
for graph_xg in ${WD}/results/graph_alignment/graphs/*.xg; do
    gene=$(basename ${graph_xg} .xg)
    
    echo ""
    echo "========================================" echo "Processing ${gene}..."
    echo "========================================"
    
    # Map paired-end reads
    echo "  Mapping PE reads..."
    ${VG} map \
        -x ${graph_xg} \
        -g ${WD}/results/graph_alignment/graphs/${gene}.gcsa \
        -f ${READS_PE1} -f ${READS_PE2} \
        -t ${THREADS} \
        > ${WD}/results/graph_alignment/mapped/${gene}.gam 2>&1
    
    # Map merged reads (CRITICAL - this was missing before!)
    echo "  Mapping merged reads..."
    ${VG} map \
        -x ${graph_xg} \
        -g ${WD}/results/graph_alignment/graphs/${gene}.gcsa \
        -f ${READS_MERGED} \
        -t ${THREADS} \
        >> ${WD}/results/graph_alignment/mapped/${gene}.gam 2>&1
    
    # Check mapping stats
    total_aln=$(${VG} stats -a ${WD}/results/graph_alignment/mapped/${gene}.gam 2>&1 | grep "Total alignments" | awk '{print $3}')
    aligned=$(${VG} stats -a ${WD}/results/graph_alignment/mapped/${gene}.gam 2>&1 | grep "Total aligned" | awk '{print $3}')
    echo "  Total alignments: ${total_aln}, aligned: ${aligned}"
    
    # Pack alignments
    echo "  Packing alignments..."
    ${VG} pack \
        -x ${graph_xg} \
        -g ${WD}/results/graph_alignment/mapped/${gene}.gam \
        -o ${WD}/results/graph_alignment/mapped/${gene}.pack \
        2>&1
    
    # Call variants with relaxed filters for ancient DNA
    echo "  Calling variants (-d 1 -s 1)..."
    ${VG} call \
        ${graph_xg} \
        -k ${WD}/results/graph_alignment/mapped/${gene}.pack \
        -d 1 \
        -s 1 \
        > ${WD}/results/graph_alignment/mapped/${gene}_raw.vcf \
        2>&1
    
    # Get first reference for consensus
    first_ref=$(grep "^>" ${WD}/results/graph_alignment/busco_regions/${gene}.fasta | head -n 1 | sed 's/^>//')
    
    # Extract first reference sequence
    awk -v ref="^>${first_ref}" 'BEGIN{p=0} $0 ~ ref {p=1; print; next} /^>/ {p=0} p' \
        ${WD}/results/graph_alignment/busco_regions/${gene}.fasta \
        > ${WD}/results/graph_alignment/mapped/${gene}_ref.fasta
    
    # Filter variants: only keep those with depth > 5
    bcftools view \
        -i 'FORMAT/DP>5' \
        ${WD}/results/graph_alignment/mapped/${gene}_raw.vcf \
        -o ${WD}/results/graph_alignment/mapped/${gene}_filtered.vcf \
        2>&1
    
    # Compress and index filtered VCF
    bgzip -c ${WD}/results/graph_alignment/mapped/${gene}_filtered.vcf > ${WD}/results/graph_alignment/mapped/${gene}.vcf.gz 2>&1
    tabix -p vcf ${WD}/results/graph_alignment/mapped/${gene}.vcf.gz 2>&1
    
    # Index reference
    samtools faidx ${WD}/results/graph_alignment/mapped/${gene}_ref.fasta 2>&1
    
    # Count variants (raw and filtered)
    raw_count=$(grep -v "^#" ${WD}/results/graph_alignment/mapped/${gene}_raw.vcf | wc -l)
    filtered_count=$(zcat ${WD}/results/graph_alignment/mapped/${gene}.vcf.gz | grep -v "^#" | wc -l)
    echo "  Variants: ${raw_count} raw, ${filtered_count} after DP>5 filter"
    
    # Apply variants to create consensus
    echo "  Creating consensus..."
    bcftools consensus \
        -f ${WD}/results/graph_alignment/mapped/${gene}_ref.fasta \
        ${WD}/results/graph_alignment/mapped/${gene}.vcf.gz \
        > ${WD}/results/graph_alignment/consensus/${gene}_Tityra.fasta \
        2>&1
    
    # Rename sequence header
    sed -i "1s/.*/>${gene}|Tityra_leucura_graph/" ${WD}/results/graph_alignment/consensus/${gene}_Tityra.fasta
    
    # Clean up large GAM file
    rm -f ${WD}/results/graph_alignment/mapped/${gene}.gam
    rm -f ${WD}/results/graph_alignment/mapped/${gene}_ref.fasta
    
    echo "  ✓ Complete: ${variant_count} variants, consensus created"
done

echo ""
echo "=================================================================="
echo "Variant calling complete for all genes"
echo "Results in: ${WD}/results/graph_alignment/consensus/"
echo "=================================================================="
