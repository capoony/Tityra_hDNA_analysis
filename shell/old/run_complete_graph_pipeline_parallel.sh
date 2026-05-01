#!/bin/bash

###############################################################################
# Complete parallelized graph-based phylogeny pipeline
# Runs all ~8100 BUSCO genes through graph alignment
###############################################################################

set -euo pipefail

WD=/media/inter/mkapun/projects/Tityra
GRAPH_JOBS=20  # Parallel jobs for graph construction
MAPPING_JOBS=10  # Parallel jobs for variant calling

echo "=================================================================="
echo "COMPLETE GRAPH-BASED PHYLOGENY PIPELINE (PARALLELIZED)"
echo "=================================================================="
echo "Working directory: ${WD}"
echo "Graph construction jobs: ${GRAPH_JOBS}"
echo "Variant calling jobs: ${MAPPING_JOBS}"
echo "=================================================================="
echo ""

# Step 1: Extract all BUSCO genes
echo "STEP 1: Extracting all BUSCO genes from references..."
bash ${WD}/shell/extract_all_busco_genes.sh

extracted=$(ls ${WD}/results/graph_alignment/busco_regions/*.fasta 2>/dev/null | wc -l)
echo "Extracted: ${extracted} genes"
echo ""

# Step 2: Build variation graphs in parallel
echo "STEP 2: Building variation graphs (parallel, ${GRAPH_JOBS} jobs)..."
bash ${WD}/shell/parallel_graph_construction.sh ${GRAPH_JOBS}

graphs=$(ls ${WD}/results/graph_alignment/graphs/*.xg 2>/dev/null | wc -l)
echo "Graphs built: ${graphs}"
echo ""

# Step 3: Map reads and call variants in parallel
echo "STEP 3: Mapping reads and calling variants (parallel, ${MAPPING_JOBS} jobs)..."
bash ${WD}/shell/parallel_variant_calling.sh ${MAPPING_JOBS}

consensus=$(ls ${WD}/results/graph_alignment/consensus/*_Tityra.fasta 2>/dev/null | wc -l)
echo "Consensus sequences: ${consensus}"
echo ""

# Step 4: Translate to proteins
echo "STEP 4: Translating consensus sequences to proteins..."
mkdir -p ${WD}/results/phylogeny_graph/busco_aa

python3 << 'PYTHON_TRANSLATE'
import os
from Bio import SeqIO
from Bio.Seq import Seq

wd = os.environ.get('WD', '/media/inter/mkapun/projects/Tityra')
consensus_dir = f"{wd}/results/graph_alignment/consensus"
output_dir = f"{wd}/results/phylogeny_graph/busco_aa"

translated = 0
for fname in os.listdir(consensus_dir):
    if not fname.endswith('_Tityra.fasta'):
        continue
    
    gene_id = fname.replace('_Tityra.fasta', '')
    
    # Read consensus nucleotide sequence
    consensus_file = f"{consensus_dir}/{fname}"
    for record in SeqIO.parse(consensus_file, "fasta"):
        # Translate to protein (frame 1)
        protein_seq = str(record.seq.translate())
        
        # Write protein sequence
        with open(f"{output_dir}/{gene_id}.faa", 'w') as out:
            out.write(f">{gene_id}|Tityra_leucura_graph\n{protein_seq}\n")
        
        translated += 1
        break

print(f"Translated {translated} sequences to protein")
PYTHON_TRANSLATE

echo ""
echo "=================================================================="
echo "PIPELINE COMPLETE!"
echo "=================================================================="
echo "BUSCO genes extracted: ${extracted}"
echo "Variation graphs built: ${graphs}"
echo "Consensus sequences created: ${consensus}"
echo ""
echo "Next steps:"
echo "1. Update final_busco_ids.txt to include graph-based genes"
echo "2. Copy graph-based proteins to phylogeny concatenation"
echo "3. Re-run phylogeny with all genomes"
echo "=================================================================="
