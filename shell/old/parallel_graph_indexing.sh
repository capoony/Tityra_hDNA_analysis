#!/bin/bash

###############################################################################
# Index all VG graphs in parallel
###############################################################################

set -euo pipefail

WD=/media/inter/mkapun/projects/Tityra
VG=/opt/bioinformatics/cactus-2.1.1/bin/vg
JOBS=${1:-20}

echo "=================================================================="
echo "Indexing variation graphs (parallel)"
echo "Jobs: ${JOBS}"
echo "=================================================================="

export WD VG

# Function to index one graph
index_graph() {
    local vg_file=$1
    local gene=$(basename ${vg_file} .vg)
    
    # Skip if indexes already exist
    if [ -f "${WD}/results/graph_alignment/graphs/${gene}.xg" ] && \
       [ -f "${WD}/results/graph_alignment/graphs/${gene}.gcsa" ]; then
        return 0
    fi
    
    # Prune complex regions
    ${VG} prune \
        ${vg_file} \
        > ${WD}/results/graph_alignment/graphs/${gene}.pruned.vg 2>/dev/null || return 1
    
    # Index for mapping
    ${VG} index \
        -x ${WD}/results/graph_alignment/graphs/${gene}.xg \
        -g ${WD}/results/graph_alignment/graphs/${gene}.gcsa \
        -k 16 \
        ${WD}/results/graph_alignment/graphs/${gene}.pruned.vg 2>/dev/null || return 1
    
    # Clean up pruned file
    rm -f ${WD}/results/graph_alignment/graphs/${gene}.pruned.vg
    
    return 0
}

export -f index_graph

# Index all graphs in parallel
find ${WD}/results/graph_alignment/graphs/ -name "*.vg" | \
    parallel -j ${JOBS} --progress index_graph {}

echo ""
echo "=================================================================="
echo "Indexing complete!"
echo "XG indexes: $(ls ${WD}/results/graph_alignment/graphs/*.xg 2>/dev/null | wc -l)"
echo "GCSA indexes: $(ls ${WD}/results/graph_alignment/graphs/*.gcsa 2>/dev/null | wc -l)"
echo "=================================================================="
