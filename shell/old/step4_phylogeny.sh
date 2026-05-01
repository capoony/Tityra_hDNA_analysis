#!/bin/bash
###############################################################################
# Step 4: Integrate Tityra Consensus into Phylogeny Pipeline
# Adds variant-called Tityra sequences to existing alignments
###############################################################################

WD=/media/inter/mkapun/projects/Tityra

CONSENSUS_DIR=${WD}/results/graph_alignment/consensus
ALIGNED_REF=${WD}/results/graph_alignment/aligned
PHYLO_NEW=${WD}/results/phylogeny_tityra

# Load RAxML module
module load Phylogeny/RAxML-2.8.10

mkdir -p ${PHYLO_NEW}/mafft
mkdir -p ${PHYLO_NEW}/phylogeny
mkdir -p ${PHYLO_NEW}/concatenated

echo "==================================================================="
echo "STEP 4: Integrating Tityra Consensus into Phylogeny"
echo "==================================================================="

# Count available sequences
NUM_TITYRA=$(ls ${CONSENSUS_DIR}/*_protein.fasta 2>/dev/null | wc -l)
NUM_REF=$(ls ${ALIGNED_REF}/*.fasta 2>/dev/null | wc -l)
echo "Tityra consensus sequences: ${NUM_TITYRA}"
echo "Reference alignments: ${NUM_REF}"

# Find overlapping genes
echo "Finding overlapping genes..."
cd ${CONSENSUS_DIR}
>  ${PHYLO_NEW}/gene_list.txt

for tityra_file in *_protein.fasta; do
    gene=$(basename ${tityra_file} _Tityra_protein.fasta)
    
    if [ -f "${ALIGNED_REF}/${gene}.fasta" ]; then
        echo "${gene}" >> ${PHYLO_NEW}/gene_list.txt
    fi
done

NUM_OVERLAP=$(wc -l < ${PHYLO_NEW}/gene_list.txt)
echo "Overlapping genes: ${NUM_OVERLAP}"

echo ""
echo "==================================================================="
echo "Adding Tityra to Alignments and Building Trees"
echo "==================================================================="

# Function to process one gene
process_gene() {
    gene=$1
    GENE_DIR=${PHYLO_NEW}/phylogeny/${gene}
    
    mkdir -p ${GENE_DIR}
    cd ${GENE_DIR}
    
    # The reference sequences are already aligned, and Tityra is the same length
    # Just combine them directly (both are nucleotide sequences)
    cat ${ALIGNED_REF}/${gene}.fasta > ${gene}_final.fasta
    cat ${CONSENSUS_DIR}/${gene}_Tityra.fasta >> ${gene}_final.fasta
    
    if [ ! -s "${gene}_final.fasta" ]; then
        return
    fi
    
    # Copy to mafft directory for reference
    cp ${gene}_final.fasta ${PHYLO_NEW}/mafft/${gene}_aln_with_tityra.fasta
    
    # Build tree with RAxML
    raxmlHPC-PTHREADS-SSE3 \
        -m PROTGAMMAWAG \
        -N 20 \
        -p 772374015 \
        -n tree \
        -s ${gene}_final.fasta \
        -T 1 \
        > /dev/null 2>&1
    
    if [ -f "RAxML_bestTree.tree" ]; then
        echo "Completed: ${gene}"
    fi
}

export -f process_gene
export PHYLO_NEW ALIGNED_REF CONSENSUS_DIR WD

# Process genes in parallel
cat ${PHYLO_NEW}/gene_list.txt | /bin/parallel -j 50 --will-cite process_gene {}

echo ""
echo "==================================================================="
echo "Building Concatenated Tree"
echo "==================================================================="

cd ${PHYLO_NEW}/concatenated

# Concatenate all aligned genes
python3 << 'PYTHON_SCRIPT'
from Bio import SeqIO
import os

phylo_dir = os.environ['PHYLO_NEW']
gene_list_file = f"{phylo_dir}/gene_list.txt"

# Read gene list
with open(gene_list_file) as f:
    genes = [line.strip() for line in f]

# Dictionary to store concatenated sequences
concat_seqs = {}

# Process each gene
for gene in genes:
    aln_file = f"{phylo_dir}/phylogeny/{gene}/{gene}_final.fasta"
    
    if not os.path.exists(aln_file):
        continue
    
    # Read alignment
    for record in SeqIO.parse(aln_file, "fasta"):
        species = record.id.split("|")[0]  # Remove everything after first |
        
        if species not in concat_seqs:
            concat_seqs[species] = []
        
        concat_seqs[species].append(str(record.seq))

# Write concatenated alignment
with open(f"{phylo_dir}/concatenated/concatenated.fasta", "w") as out:
    for species, seqs in concat_seqs.items():
        concat_seq = "".join(seqs)
        out.write(f">{species}\n{concat_seq}\n")

print(f"Concatenated {len(genes)} genes for {len(concat_seqs)} species")
PYTHON_SCRIPT

# Build concatenated tree with RAxML
raxmlHPC-PTHREADS-SSE3 \
    -m PROTGAMMAWAG \
    -N 20 \
    -p 772374015 \
    -n concatenated \
    -s concatenated.fasta \
    -T 50 \
    > /dev/null 2>&1

echo ""
echo "==================================================================="
echo "STEP 4 Complete"
echo "==================================================================="
echo "Individual gene trees: ${PHYLO_NEW}/phylogeny/"
echo "Concatenated tree: ${PHYLO_NEW}/concatenated/RAxML_bestTree.concatenated"
echo ""
echo "Key files:"
echo "  - ${PHYLO_NEW}/concatenated/RAxML_bestTree.concatenated (main tree)"
echo "  - ${PHYLO_NEW}/gene_list.txt (${NUM_OVERLAP} genes included)"
echo "  - ${PHYLO_NEW}/mafft/ (alignments with Tityra)"
echo "==================================================================="
