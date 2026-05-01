#!/bin/bash
###############################################################################
# BUSCO-based Multi-Species Phylogeny Pipeline
# Uses complete and partial BUSCO genes from aves_odb10 database
# Builds individual gene trees with IQ-TREE and species tree with ASTRAL
###############################################################################

WD=/media/inter/mkapun/projects/Tityra

# Define directories - ALL species from phylogeny BUSCO folder
BUSCO_DIR="${WD}/results/phylogeny/BUSCO"
RESULTS="${WD}/results/BUSCO_aves_phylogeny"

# Create result subdirectories
mkdir -p ${RESULTS}/gene_sequences
mkdir -p ${RESULTS}/alignments
mkdir -p ${RESULTS}/gene_trees
mkdir -p ${RESULTS}/astral

echo "==================================================================="
echo "BUSCO-based Phylogenomic Analysis"
echo "==================================================================="
echo "Analysis Directory: ${RESULTS}"
echo ""

###############################################################################
# Step 1: Identify all species with BUSCO data
###############################################################################

echo "Step 1: Identifying species with BUSCO data..."
> ${RESULTS}/species_list.txt

# Add all birds from phylogeny BUSCO folder (including Tityra)
for dir in ${BUSCO_DIR}/*/; do
    species=$(basename "$dir")
    # Skip if not a directory or is busco_downloads
    if [ ! -d "${dir}" ] || [ "${species}" == "busco_downloads" ]; then
        continue
    fi
    
    # Check if BUSCO run completed successfully
    if [ -f "${dir}run_aves_odb10/full_table.tsv" ]; then
        echo "${species}" >> ${RESULTS}/species_list.txt
        echo "  ✓ ${species}"
    fi
done

NUM_SPECIES=$(wc -l < ${RESULTS}/species_list.txt)
echo ""
echo "Total species: ${NUM_SPECIES}"

###############################################################################
# Step 2: Extract Complete and Duplicated BUSCO genes from all species
###############################################################################

echo ""
echo "==================================================================="
echo "Step 2: Extracting BUSCO gene IDs"
echo "==================================================================="

# Collect all Complete and Duplicated BUSCO gene IDs across all species
> ${RESULTS}/all_busco_genes.txt

while read species; do
    busco_table="${BUSCO_DIR}/${species}/run_aves_odb10/full_table.tsv"
    
    if [ -f "${busco_table}" ]; then
        # Extract Complete and Duplicated genes
        grep -v "^#" "${busco_table}" | \
            awk '$2=="Complete" || $2=="Duplicated" {print $1}' >> ${RESULTS}/all_busco_genes.txt
    fi
done < ${RESULTS}/species_list.txt

# Count occurrences and find genes present in all species
sort ${RESULTS}/all_busco_genes.txt | uniq -c | \
    awk -v ns="${NUM_SPECIES}" '$1 == ns {print $2}' > ${RESULTS}/universal_genes.txt

NUM_GENES=$(wc -l < ${RESULTS}/universal_genes.txt)
echo "Genes present in all ${NUM_SPECIES} species: ${NUM_GENES}"

###############################################################################
# Step 3: Extract protein sequences for each gene from each species
###############################################################################

echo ""
echo "==================================================================="
echo "Step 3: Extracting protein sequences"
echo "==================================================================="

# Function to extract protein sequence for one gene from one species
extract_sequence() {
    species=$1
    gene=$2
    
    # All species (including Tityra) use standard BUSCO structure with .faa protein files
    busco_seqs="${BUSCO_DIR}/${species}/run_aves_odb10/busco_sequences/single_copy_busco_sequences"
    seq_file="${busco_seqs}/${gene}.faa"
    
    if [ -f "${seq_file}" ]; then
        echo ">${species}" >> ${RESULTS}/gene_sequences/${gene}.fasta
        grep -v "^>" "${seq_file}" | tr -d '\n' >> ${RESULTS}/gene_sequences/${gene}.fasta
        echo "" >> ${RESULTS}/gene_sequences/${gene}.fasta
    fi
}

export -f extract_sequence
export BUSCO_DIR RESULTS

# Extract sequences for all genes in parallel
echo "Extracting ${NUM_GENES} genes from ${NUM_SPECIES} species..."

while read gene; do
    > ${RESULTS}/gene_sequences/${gene}.fasta
    while read species; do
        extract_sequence "${species}" "${gene}"
    done < ${RESULTS}/species_list.txt
done < ${RESULTS}/universal_genes.txt

echo "Sequence extraction complete"

###############################################################################
# Step 4: Align each gene with MAFFT
###############################################################################

echo ""
echo "==================================================================="
echo "Step 4: Aligning genes with MAFFT"
echo "==================================================================="

# Load MAFFT
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate mafft-7.487

# Function to align one gene
align_gene() {
    gene=$1
    
    gene_fasta="${RESULTS}/gene_sequences/${gene}.fasta"
    alignment="${RESULTS}/alignments/${gene}_aln.fasta"
    
    # Check if sequence file exists and has sequences
    if [ ! -f "${gene_fasta}" ]; then
        return
    fi
    
    num_seqs=$(grep -c "^>" "${gene_fasta}")
    if [ ${num_seqs} -lt 4 ]; then
        return
    fi
    
    # Align with MAFFT
    mafft --auto --quiet "${gene_fasta}" > "${alignment}" 2>/dev/null
    
    if [ -f "${alignment}" ]; then
        echo "Aligned: ${gene} (${num_seqs} species)"
    fi
}

export -f align_gene
export RESULTS

# Align genes in parallel
cat ${RESULTS}/universal_genes.txt | /bin/parallel -j 50 --will-cite align_gene {}

NUM_ALIGNMENTS=$(ls ${RESULTS}/alignments/*.fasta 2>/dev/null | wc -l)
echo "Generated ${NUM_ALIGNMENTS} alignments"

###############################################################################
# Step 5: Build individual gene trees with IQ-TREE
###############################################################################

echo ""
echo "==================================================================="
echo "Step 5: Building gene trees with IQ-TREE"
echo "==================================================================="

# Load IQ-TREE (use latest version)
conda activate iq-tree-3.0.1

# Function to build tree for one gene
build_tree() {
    gene=$1
    
    alignment="${RESULTS}/alignments/${gene}_aln.fasta"
    
    if [ ! -f "${alignment}" ]; then
        return
    fi
    
    # Create gene-specific directory for IQ-TREE output
    gene_dir="${RESULTS}/gene_trees/${gene}"
    mkdir -p "${gene_dir}"
    cd "${gene_dir}"
    
    # Build tree with IQ-TREE
    # -s: alignment file
    # -m: model selection (automatic)
    # -bb: ultrafast bootstrap (1000 replicates)
    # -nt: number of threads
    # -quiet: suppress output
    iqtree -s "${alignment}" -m MFP -bb 1000 -nt 1 -quiet > /dev/null 2>&1
    
    if [ -f "${gene}.fasta.treefile" ]; then
        # Copy tree file to common location for ASTRAL
        cp "${gene}_aln.fasta.treefile" "${RESULTS}/astral/${gene}.tree"
        echo "Tree built: ${gene}"
    fi
}

export -f build_tree
export RESULTS

# Build trees in parallel
cat ${RESULTS}/universal_genes.txt | /bin/parallel -j 50 --will-cite build_tree {}

NUM_TREES=$(ls ${RESULTS}/astral/*.tree 2>/dev/null | wc -l)
echo "Generated ${NUM_TREES} gene trees"

###############################################################################
# Step 6: Build species tree with ASTRAL
###############################################################################

echo ""
echo "==================================================================="
echo "Step 6: Building species tree with ASTRAL"
echo "==================================================================="

# Concatenate all gene trees into single file
cat ${RESULTS}/astral/*.tree > ${RESULTS}/astral/all_gene_trees.trees

# Load ASTRAL
module load Phylogeny/Astral-5.7.7

# Run ASTRAL to build species tree
cd ${RESULTS}/astral

echo "Running ASTRAL on ${NUM_TREES} gene trees..."
java -jar ${ASTRAL} \
    -i all_gene_trees.trees \
    -o species_tree.tre \
    2> astral.log

if [ -f "species_tree.tre" ]; then
    echo "✓ Species tree built successfully"
    echo ""
    echo "Species tree:"
    cat species_tree.tre
else
    echo "✗ ERROR: ASTRAL failed to build species tree"
    echo "Check astral.log for details"
fi

###############################################################################
# Step 7: Generate summary statistics
###############################################################################

echo ""
echo "==================================================================="
echo "Summary Statistics"
echo "==================================================================="

echo "Species included: ${NUM_SPECIES}"
echo "Universal BUSCO genes: ${NUM_GENES}"
echo "Successful alignments: ${NUM_ALIGNMENTS}"
echo "Gene trees built: ${NUM_TREES}"
echo ""
echo "Results directory: ${RESULTS}"
echo "  - Gene sequences: ${RESULTS}/gene_sequences/"
echo "  - Alignments: ${RESULTS}/alignments/"
echo "  - Gene trees: ${RESULTS}/gene_trees/"
echo "  - Species tree: ${RESULTS}/astral/species_tree.tre"
echo "==================================================================="
