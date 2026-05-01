#!/bin/bash

###############################################################################
# Parallel graph construction for all BUSCO genes
# Usage: bash parallel_graph_construction.sh <num_jobs>
###############################################################################

set -euo pipefail

WD=/media/inter/mkapun/projects/Tityra
VG=/opt/bioinformatics/cactus-2.1.1/bin/vg
THREADS=${1:-20}  # Number of parallel jobs

echo "=================================================================="
echo "Building variation graphs for all BUSCO genes (parallel)"
echo "Jobs: ${THREADS}"
echo "=================================================================="

mkdir -p ${WD}/results/graph_alignment/graphs
mkdir -p ${WD}/results/graph_alignment/aligned

# Activate conda for MAFFT
source /opt/mamba/etc/profile.d/conda.sh
conda activate mafft-7.487

export WD VG

# Function to process one gene
process_gene() {
    local busco_fasta=$1
    local gene=$(basename ${busco_fasta} .fasta)
    
    # Skip if graph already exists
    if [ -f "${WD}/results/graph_alignment/graphs/${gene}.xg" ] && \
       [ -f "${WD}/results/graph_alignment/graphs/${gene}.gcsa" ]; then
        return 0
    fi
    
    # Align with MAFFT
    mafft --auto --quiet ${busco_fasta} > ${WD}/results/graph_alignment/aligned/${gene}_aligned.fasta 2>/dev/null || return 1
    
    # Build variation graph
    ${VG} construct \
        -r ${WD}/results/graph_alignment/aligned/${gene}_aligned.fasta \
        -M ${gene} \
        > ${WD}/results/graph_alignment/graphs/${gene}.vg 2>/dev/null || return 1
    
    # Prune complex regions
    ${VG} prune \
        ${WD}/results/graph_alignment/graphs/${gene}.vg \
        > ${WD}/results/graph_alignment/graphs/${gene}.pruned.vg 2>/dev/null || return 1
    
    # Index for mapping
    ${VG} index \
        -x ${WD}/results/graph_alignment/graphs/${gene}.xg \
        -g ${WD}/results/graph_alignment/graphs/${gene}.gcsa \
        -k 16 \
        ${WD}/results/graph_alignment/graphs/${gene}.pruned.vg 2>/dev/null || return 1
    
    # Clean up intermediate files
    rm -f ${WD}/results/graph_alignment/graphs/${gene}.vg
    rm -f ${WD}/results/graph_alignment/graphs/${gene}.pruned.vg
    
    return 0
}

export -f process_gene

# Process all genes in parallel
find ${WD}/results/graph_alignment/busco_regions/ -name "*.fasta" | \
    parallel -j ${THREADS} --progress process_gene {}

conda deactivate

echo ""
echo "=================================================================="
echo "Graph construction complete!"
echo "Graphs built: $(ls ${WD}/results/graph_alignment/graphs/*.xg 2>/dev/null | wc -l)"
echo "=================================================================="
