#!/bin/bash

###############################################################################
# Combine Mitochondrial Gene Sequences from Assembly-based MitoFinder
# Description: Join gene-specific sequences from downloaded data and add
#              Tityra sequences from MitoFinder assembly results
###############################################################################

set -euo pipefail

# Set working directory
WD=/media/inter/mkapun/projects/Tityra

# Input directories
MITO_DIR="${WD}/data/Mito"
TITYRA_GENES_FILE="${WD}/results/mitofinder_assembly/Tityra_leucura/Tityra_leucura/Tityra_leucura_MitoFinder_mitfi_Final_Results/Tityra_leucura_mtDNA_contig_1_genes_NT.fasta"

# Output directory
OUTPUT_DIR="${WD}/results/combined_mito_genes_assembly"
mkdir -p "${OUTPUT_DIR}"

# Log file
LOG="${WD}/logs/combine_mito_genes_assembly.log"
mkdir -p "${WD}/logs"

exec > >(tee -a "${LOG}") 2>&1

echo "==================================================================="
echo "Combining Mitochondrial Gene Sequences (Assembly-based)"
echo "Started: $(date)"
echo "==================================================================="

# Check if Tityra genes file exists
if [ ! -f "${TITYRA_GENES_FILE}" ]; then
    echo "Error: Tityra genes file not found: ${TITYRA_GENES_FILE}"
    exit 1
fi

echo "Using Tityra genes from: ${TITYRA_GENES_FILE}"

echo ""
echo "==================================================================="
echo "Combining sequences for each gene"
echo "==================================================================="

# Get list of unique genes from downloaded data
GENES=$(ls -1 "${MITO_DIR}"/*.fasta | xargs -n1 basename | cut -d'_' -f1 | sort -u)

# Map synonym genes to their primary names
declare -A GENE_MAP
GENE_MAP["COX1"]="COI"
GENE_MAP["COX2"]="CO2"
GENE_MAP["COX3"]="CO3"

# Track which genes we've already processed
declare -A PROCESSED_GENES

for GENE in ${GENES}; do
    # Check if this is a synonym that maps to another gene
    if [ -n "${GENE_MAP[$GENE]+x}" ]; then
        PRIMARY_GENE="${GENE_MAP[$GENE]}"
        # This is a synonym, use the primary gene name
        if [ -n "${PROCESSED_GENES[$PRIMARY_GENE]+x}" ]; then
            # Primary gene already processed, skip
            echo "Skipping ${GENE} (already processed as ${PRIMARY_GENE})"
            continue
        fi
        GENE="${PRIMARY_GENE}"
    fi
    
    # Check if we've already processed this gene
    if [ -n "${PROCESSED_GENES[$GENE]+x}" ]; then
        continue
    fi
    PROCESSED_GENES[$GENE]=1
    
    echo ""
    echo "-------------------------------------------------------------------"
    echo "Processing ${GENE}..."
    echo "-------------------------------------------------------------------"
    
    OUTPUT_FILE="${OUTPUT_DIR}/${GENE}_combined.fasta"
    > "${OUTPUT_FILE}"  # Clear file
    
    # Collect all genus files for this gene AND its synonyms
    GENUS_FILES=""
    for SEARCH_GENE in "${GENE}" $(echo "${!GENE_MAP[@]}" | tr ' ' '\n' | grep -v "^$" | while read key; do [ "${GENE_MAP[$key]}" = "${GENE}" ] && echo "$key"; done); do
        FILES=$(ls -1 "${MITO_DIR}/${SEARCH_GENE}_"*.fasta 2>/dev/null || true)
        if [ -n "${FILES}" ]; then
            GENUS_FILES="${GENUS_FILES} ${FILES}"
        fi
    done
    
    if [ -z "${GENUS_FILES}" ]; then
        echo "No sequences found for ${GENE}, skipping..."
        continue
    fi
    
    # Process each genus file
    for FILE in ${GENUS_FILES}; do
        GENUS=$(basename "${FILE}" | cut -d'_' -f2 | cut -d'.' -f1)
        GENE_NAME=$(basename "${FILE}" | cut -d'_' -f1)
        echo "  Processing ${GENUS} (${GENE_NAME})..."
        
        # Get unique species and randomly select up to 3 samples per species
        python3 << PYTHON_EOF
import sys
import random
from collections import defaultdict

# Set random seed for reproducibility
random.seed(42)

sequences = defaultdict(list)

# Read sequences
with open("${FILE}") as f:
    header = None
    seq = []
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if header:
                sequences[species].append((header, ''.join(seq)))
            header = line
            # Extract species name - first two words (genus + species)
            parts = line.replace(">", "").split()
            species = "_".join(parts[:2]) if len(parts) >= 2 else parts[0]
            seq = []
        else:
            seq.append(line)
    if header:
        sequences[species].append((header, ''.join(seq)))

# Write up to 3 randomly selected samples per species
output_file = "${OUTPUT_FILE}"
count = 0
with open(output_file, 'a') as out:
    for species, seqs in sequences.items():
        # Randomly select up to 3 samples
        selected = random.sample(seqs, min(3, len(seqs)))
        for header, seq in selected:
            out.write(f"{header}\\n{seq}\\n")
            count += 1

print(f"  Added {count} sequences from ${GENUS}")
PYTHON_EOF
    done
    
    # Deduplicate and ensure only 3 samples per species across all files
    echo "  Deduplicating and limiting to 3 samples per species..."
    python3 << DEDUP_EOF
import random
from collections import defaultdict

# Set random seed for reproducibility
random.seed(42)

input_file = "${OUTPUT_FILE}"
sequences = defaultdict(list)

# Read all sequences
with open(input_file, 'r') as f:
    header = None
    seq = []
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if header:
                sequences[species].append((header, ''.join(seq)))
            header = line
            # Extract species name - everything before the first space
            species = line.replace(">", "").split()[0] if line.strip() else ""
            seq = []
        else:
            seq.append(line)
    if header:
        sequences[species].append((header, ''.join(seq)))

# Write deduplicated sequences (max 3 per species)
with open(input_file, 'w') as out:
    total = 0
    for species, seqs in sequences.items():
        # Remove duplicates and randomly select up to 3
        unique_seqs = list({seq: (h, seq) for h, seq in seqs}.values())
        selected = random.sample(unique_seqs, min(3, len(unique_seqs)))
        for header, seq in selected:
            out.write(f"{header}\\n{seq}\\n")
            total += 1
    print(f"  Total unique sequences after deduplication: {total}")

DEDUP_EOF
    
    # Add Tityra sequence for this gene
    echo "  Extracting Tityra leucura sequence..."
    
    # Extract the specific gene from Tityra genes file
    # Gene names in the file might have variations (COX1/COI, etc.)
    python3 << TITYRA_EOF
import sys

gene = "${GENE}"
input_file = "${TITYRA_GENES_FILE}"
output_file = "${OUTPUT_FILE}"

# Gene name aliases
aliases = {
    "COI": ["COI", "CO1", "COX1"],
    "CO2": ["CO2", "COII", "COX2"],
    "CO3": ["CO3", "COIII", "COX3"],
    "CYTB": ["CYB", "CYTB", "CYTOCHROME_B"],
    "ATP6": ["ATP6", "ATPASE_6", "ATPASE6"],
    "ATP8": ["ATP8", "ATPASE_8", "ATPASE8"],
    "ND1": ["ND1", "NAD1"],
    "ND2": ["ND2", "NAD2"],
    "ND3": ["ND3", "NAD3"],
    "ND4": ["ND4", "NAD4"],
    "ND4L": ["ND4L", "NAD4L"],
    "ND5": ["ND5", "NAD5"],
    "ND6": ["ND6", "NAD6"]
}

# Get list of possible names for this gene
possible_names = aliases.get(gene, [gene])

# Read Tityra genes file and find matching gene
written = False
with open(input_file, 'r') as f:
    header = None
    seq = []
    found = False
    
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            # Write previous sequence if it matched
            if found and header and seq and not written:
                with open(output_file, 'a') as out:
                    out.write(f">Tityra_leucura\\n{''.join(seq)}\\n")
                print(f"  Added Tityra leucura {gene} sequence")
                written = True
                break
            
            # Check if this header matches our gene
            header = line
            gene_name = line.split()[0].replace(">", "").upper()
            
            found = False
            for possible in possible_names:
                if possible.upper() in gene_name or gene_name in possible.upper():
                    found = True
                    seq = []
                    break
        elif found:
            seq.append(line)
    
    # Check last sequence (if not already written)
    if found and header and seq and not written:
        with open(output_file, 'a') as out:
            out.write(f">Tityra_leucura\\n{''.join(seq)}\\n")
        print(f"  Added Tityra leucura {gene} sequence")

TITYRA_EOF
    
    # Count final sequences
    SEQ_COUNT=$(grep -c "^>" "${OUTPUT_FILE}" || echo 0)
    echo "  Total sequences in ${GENE}_combined.fasta: ${SEQ_COUNT}"
done

echo ""
echo "==================================================================="
echo "Pipeline completed: $(date)"
echo "==================================================================="
echo ""
echo "Combined gene files saved to: ${OUTPUT_DIR}"
echo ""

# Summary
echo "Summary:"
for FILE in "${OUTPUT_DIR}"/*_combined.fasta; do
    if [ -f "${FILE}" ]; then
        GENE=$(basename "${FILE}" | sed 's/_combined.fasta//')
        COUNT=$(grep -c "^>" "${FILE}" || echo 0)
        echo "  ${GENE}: ${COUNT} sequences"
    fi
done
