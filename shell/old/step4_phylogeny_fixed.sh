#!/bin/bash
###############################################################################
# Step 4: Integrate Tityra Consensus into Multi-Species Phylogeny
# - Uses ALL species from /results/phylogeny/BUSCO
# - Adds Tityra consensus sequences from graph alignment
# - Builds individual gene trees and concatenated species tree
###############################################################################

WD=/media/inter/mkapun/projects/Tityra

CONSENSUS_DIR=${WD}/results/graph_alignment/consensus
BUSCO_DIR=${WD}/results/phylogeny/BUSCO
PHYLO_NEW=${WD}/results/phylogeny_tityra

# Load RAxML module
module load Phylogeny/RAxML-2.8.10

mkdir -p ${PHYLO_NEW}/gene_alignments
mkdir -p ${PHYLO_NEW}/individual_trees
mkdir -p ${PHYLO_NEW}/concatenated

echo "==================================================================="
echo "STEP 4: Multi-Species Phylogeny with Tityra Consensus"
echo "==================================================================="

# Identify all species with completed BUSCO runs
echo "Finding species with BUSCO data..."
> ${PHYLO_NEW}/species_list.txt

for dir in ${BUSCO_DIR}/*/; do
    species=$(basename "$dir")
    # Skip if not a directory or is busco_downloads
    if [ ! -d "${dir}" ] || [ "${species}" == "busco_downloads" ]; then
        continue
    fi
    
    # Check if BUSCO run completed successfully
    if [ -f "${dir}run_aves_odb10/full_table.tsv" ]; then
        echo "${species}" >> ${PHYLO_NEW}/species_list.txt
    fi
done

NUM_SPECIES=$(wc -l < ${PHYLO_NEW}/species_list.txt)
echo "Found ${NUM_SPECIES} species with BUSCO data:"
cat ${PHYLO_NEW}/species_list.txt

# Identify BUSCO genes that are complete in ALL species
echo ""
echo "==================================================================="
echo "Finding BUSCO genes present in all ${NUM_SPECIES} species"
echo "==================================================================="

# Collect all complete BUSCO gene IDs across all species
> ${PHYLO_NEW}/all_complete_genes.txt

while read species; do
    busco_table="${BUSCO_DIR}/${species}/run_aves_odb10/full_table.tsv"
    if [ -f "${busco_table}" ]; then
        grep -v "^#" "${busco_table}" | awk '$2=="Complete" {print $1}' >> ${PHYLO_NEW}/all_complete_genes.txt
    fi
done < ${PHYLO_NEW}/species_list.txt

# Count occurrences and find genes present in all species
sort ${PHYLO_NEW}/all_complete_genes.txt | uniq -c | \
    awk -v ns="${NUM_SPECIES}" '$1 == ns {print $2}' > ${PHYLO_NEW}/universal_genes.txt

NUM_GENES=$(wc -l < ${PHYLO_NEW}/universal_genes.txt)
echo "Found ${NUM_GENES} BUSCO genes present in all ${NUM_SPECIES} species"

# Filter to only genes for which we have Tityra consensus
echo "Filtering to genes with Tityra consensus sequences..."
> ${PHYLO_NEW}/final_gene_list.txt

while read gene; do
    if [ -f "${CONSENSUS_DIR}/${gene}_Tityra.fasta" ]; then
        echo "${gene}" >> ${PHYLO_NEW}/final_gene_list.txt
    fi
done < ${PHYLO_NEW}/universal_genes.txt

NUM_FINAL=$(wc -l < ${PHYLO_NEW}/final_gene_list.txt)
echo "Final gene set: ${NUM_FINAL} genes with Tityra consensus"

echo ""
echo "==================================================================="
echo "Building Individual Gene Trees"
echo "==================================================================="

# Function to process one gene
process_gene() {
    gene=$1
    
    GENE_DIR=${PHYLO_NEW}/individual_trees/${gene}
    mkdir -p ${GENE_DIR}
    cd ${GENE_DIR}
    
    # Create alignment file with all species
    > ${gene}_alignment.fasta
    
    # Add each species' sequence
    while read species; do
        busco_dir="${BUSCO_DIR}/${species}/run_aves_odb10/busco_sequences/single_copy_busco_sequences"
        
        # Find the BUSCO sequence file (nucleotide)
        seq_file=$(find "${busco_dir}" -name "${gene}.fna" 2>/dev/null | head -1)
        
        if [ -f "${seq_file}" ]; then
            # Add species name as header and sequence
            echo ">${species}" >> ${gene}_alignment.fasta
            grep -v "^>" "${seq_file}" | tr -d '\n' >> ${gene}_alignment.fasta
            echo "" >> ${gene}_alignment.fasta
        fi
    done < ${PHYLO_NEW}/species_list.txt
    
    # Add Tityra consensus sequence
    if [ -f "${CONSENSUS_DIR}/${gene}_Tityra.fasta" ]; then
        cat "${CONSENSUS_DIR}/${gene}_Tityra.fasta" >> ${gene}_alignment.fasta
    fi
    
    # Check if alignment has sequences
    num_seqs=$(grep -c "^>" ${gene}_alignment.fasta)
    if [ ${num_seqs} -lt 4 ]; then
        echo "Skipped: ${gene} (only ${num_seqs} sequences)"
        return
    fi
    
    # Copy to gene_alignments directory
    cp ${gene}_alignment.fasta ${PHYLO_NEW}/gene_alignments/${gene}_aln.fasta
    
    # Build tree with RAxML using GTRGAMMA model for nucleotide data
    raxmlHPC-PTHREADS-SSE3 \
        -m GTRGAMMA \
        -N 20 \
        -p 772374015 \
        -n tree \
        -s ${gene}_alignment.fasta \
        -T 1 \
        > raxml.log 2>&1
    
    if [ -f "RAxML_bestTree.tree" ]; then
        echo "Completed: ${gene} (${num_seqs} species)"
    else
        echo "Failed: ${gene}"
    fi
}

export -f process_gene
export PHYLO_NEW BUSCO_DIR CONSENSUS_DIR WD

# Process genes in parallel
cat ${PHYLO_NEW}/final_gene_list.txt | /bin/parallel -j 50 --will-cite process_gene {}

echo ""
echo "==================================================================="
echo "Building Concatenated Species Tree"
echo "==================================================================="

cd ${PHYLO_NEW}/concatenated

# Python script to concatenate alignments with proper gap padding
python3 << 'PYTHON_SCRIPT'
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

phylo_dir = os.environ['PHYLO_NEW']
gene_list_file = f"{phylo_dir}/final_gene_list.txt"
species_list_file = f"{phylo_dir}/species_list.txt"

# Read gene list
with open(gene_list_file) as f:
    genes = [line.strip() for line in f if line.strip()]

# Read species list
with open(species_list_file) as f:
    all_species = [line.strip() for line in f if line.strip()]

# Add Tityra to species list
all_species.append("Tityra_leucura")

print(f"Concatenating {len(genes)} genes for {len(all_species)} species")

# Dictionary to store concatenated sequences for each species
concat_seqs = {species: [] for species in all_species}
gene_lengths = []

# Process each gene
for gene in genes:
    aln_file = f"{phylo_dir}/individual_trees/{gene}/{gene}_alignment.fasta"
    
    if not os.path.exists(aln_file):
        continue
    
    # Read alignment
    seqs_in_gene = {}
    for record in SeqIO.parse(aln_file, "fasta"):
        species = record.id
        seqs_in_gene[species] = str(record.seq)
    
    # Determine sequence length for this gene
    if seqs_in_gene:
        gene_length = len(next(iter(seqs_in_gene.values())))
        gene_lengths.append(gene_length)
        
        # Add sequence for each species (or gaps if missing)
        for species in all_species:
            if species in seqs_in_gene:
                concat_seqs[species].append(seqs_in_gene[species])
            else:
                # Add gap padding for missing genes
                concat_seqs[species].append("-" * gene_length)

# Write concatenated alignment
with open(f"{phylo_dir}/concatenated/concatenated.fasta", "w") as out:
    for species in sorted(all_species):
        concat_seq = "".join(concat_seqs[species])
        out.write(f">{species}\n{concat_seq}\n")

total_length = sum(gene_lengths)
print(f"Concatenated alignment length: {total_length} bp ({len(genes)} genes)")
print(f"Species: {len(all_species)}")

# Verify all sequences have same length
seq_lengths = {species: len("".join(seqs)) for species, seqs in concat_seqs.items()}
if len(set(seq_lengths.values())) == 1:
    print(f"✓ All sequences have same length: {list(seq_lengths.values())[0]} bp")
else:
    print("✗ ERROR: Sequences have different lengths!")
    for species, length in sorted(seq_lengths.items()):
        print(f"  {species}: {length} bp")
PYTHON_SCRIPT

# Check if concatenated alignment was created successfully
if [ ! -f "concatenated.fasta" ]; then
    echo "ERROR: Failed to create concatenated alignment"
    exit 1
fi

# Build concatenated tree with RAxML using GTRGAMMA model
echo "Building concatenated species tree with RAxML..."
raxmlHPC-PTHREADS-SSE3 \
    -m GTRGAMMA \
    -N 20 \
    -p 772374015 \
    -n concatenated \
    -s concatenated.fasta \
    -T 50

if [ -f "RAxML_bestTree.concatenated" ]; then
    echo "✓ Concatenated tree built successfully"
else
    echo "✗ ERROR: Concatenated tree building failed"
    echo "Check RAxML output above for errors"
fi

echo ""
echo "==================================================================="
echo "STEP 4 Complete"
echo "==================================================================="
echo "Species included: ${NUM_SPECIES} + Tityra_leucura"
echo "Genes analyzed: ${NUM_FINAL}"
echo ""
echo "Results:"
echo "  - Individual gene trees: ${PHYLO_NEW}/individual_trees/"
echo "  - Gene alignments: ${PHYLO_NEW}/gene_alignments/"
echo "  - Concatenated alignment: ${PHYLO_NEW}/concatenated/concatenated.fasta"
echo "  - Species tree: ${PHYLO_NEW}/concatenated/RAxML_bestTree.concatenated"
echo "==================================================================="
