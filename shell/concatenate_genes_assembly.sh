#!/bin/bash

###############################################################################
# Concatenate Mitochondrial Genes by Species (Rawreads-based)
# Description: Create concatenated gene alignments per species, allowing
#              different vouchers/specimens for each gene to maximize
#              species representation with multiple sequences
###############################################################################

set -euo pipefail

# Set working directory
WD=/media/inter/mkapun/projects/Tityra

# Input directory (from aligned genes)
INPUT_DIR="${WD}/results/phylogenetic_analysis_assembly/alignments"

# Output directory
OUTPUT_DIR="${WD}/results/concatenated_genes_assembly"
mkdir -p "${OUTPUT_DIR}"

# Log file
LOG="${WD}/logs/concatenate_genes_assembly.log"
mkdir -p "${WD}/logs"

exec > >(tee -a "${LOG}") 2>&1

echo "==================================================================="
echo "Concatenating Mitochondrial Genes by Species (Assembly-based)"
echo "Started: $(date)"
echo "==================================================================="

# List of genes to concatenate (4 mitochondrial genes)
# Mitochondrial: COI, CYTB, ND2, CO2 (ATP6 excluded)
GENES=(COI CYTB ND2 CO2)

echo "Using maximum of ${#GENES[@]} genes: ${GENES[*]}"
echo "Will ensure at least one species per genus"
echo "Will filter out species with >50% missing data (within genus priority)"
echo ""

# Check that aligned gene files exist
echo "Checking for aligned gene files..."
AVAILABLE_GENES=()
for GENE in "${GENES[@]}"; do
    GENE_FILE="${INPUT_DIR}/${GENE}_aligned_gblocks.fasta"
    if [ -f "${GENE_FILE}" ]; then
        echo "  ✓ ${GENE}: ${GENE_FILE}"
        AVAILABLE_GENES+=("${GENE}")
    else
        echo "  ✗ ${GENE}: Not found"
    fi
done

if [ ${#AVAILABLE_GENES[@]} -eq 0 ]; then
    echo "Error: No aligned gene files found"
    exit 1
fi

echo ""
echo "Will concatenate ${#AVAILABLE_GENES[@]} genes: ${AVAILABLE_GENES[*]}"

###############################################################################
# Extract species information and create concatenated sequences
###############################################################################

echo ""
echo "==================================================================="
echo "Building species-level concatenated sequences"
echo "==================================================================="

# Export available genes for Python script
export AVAILABLE_GENES_STR="${AVAILABLE_GENES[*]}"
export INPUT_DIR
export OUTPUT_DIR

# Python script to concatenate sequences by species
python3 << 'PYTHON_SCRIPT'
import os
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Get directories from environment
input_dir = os.environ.get('INPUT_DIR')
output_dir = os.environ.get('OUTPUT_DIR')

# Get available genes from environment
genes_str = os.environ.get('AVAILABLE_GENES_STR', '')
genes = genes_str.split()

if not genes:
    print("Error: No genes specified")
    sys.exit(1)

print(f"Processing {len(genes)} genes: {', '.join(genes)}")
print()

# Dictionary to store sequences by species and gene
# species_data[species][gene] = [(voucher, sequence), ...]
species_data = defaultdict(lambda: defaultdict(list))

# Dictionary to store alignment lengths for each gene
gene_lengths = {}

# Gene name synonyms for matching downloaded files
gene_synonyms = {
    'RAG1': ['RAG1', 'RAG-1'],
    'RAG2': ['RAG2', 'RAG-2']
}

# Read all sequences and organize by species
for gene in genes:
    # Try primary gene name first
    gene_file = f"{input_dir}/{gene}_aligned_gblocks.fasta"
    
    # If gene not found and has synonyms, try those
    if not os.path.exists(gene_file) and gene in gene_synonyms:
        for syn in gene_synonyms[gene]:
            alt_file = f"{input_dir}/{syn}_aligned_gblocks.fasta"
            if os.path.exists(alt_file):
                gene_file = alt_file
                break
    
    if not os.path.exists(gene_file):
        print(f"Warning: {gene} file not found, skipping")
        continue
    
    print(f"Reading {gene}...")
    
    # Track alignment length for this gene
    alignment_length = None
    
    for record in SeqIO.parse(gene_file, "fasta"):
        # Extract species name (first two parts of sequence ID)
        parts = record.id.split('_')
        if len(parts) >= 2:
            species = f"{parts[0]}_{parts[1]}"
            
            # Extract accession from description if available
            # Format is typically: "species_voucher accession ..."
            accession = "NA"
            if record.description:
                desc_parts = record.description.split()
                # Look for accession pattern (e.g., KX902341, MK880909)
                for part in desc_parts:
                    # Check if it looks like an accession (letters followed by numbers)
                    if len(part) > 3 and any(c.isalpha() for c in part) and any(c.isdigit() for c in part):
                        if not part.startswith('>'):
                            accession = part
                            break
            
            voucher_info = f"{record.id}|{accession}"
            
            # Store sequence length
            if alignment_length is None:
                alignment_length = len(record.seq)
            elif len(record.seq) != alignment_length:
                print(f"  Warning: Inconsistent alignment length in {gene}: {len(record.seq)} vs {alignment_length}")
            
            species_data[species][gene].append((voucher_info, str(record.seq)))
    
    gene_lengths[gene] = alignment_length
    print(f"  {gene}: {alignment_length} bp, {len(species_data)} species")

print()
print("===================================================================")
print("Species summary:")
print("===================================================================")

# Count how many genes each species has
species_gene_counts = {}
for species in sorted(species_data.keys()):
    gene_count = len(species_data[species])
    species_gene_counts[species] = gene_count
    print(f"{species:30s}: {gene_count:2d}/{len(genes)} genes")

print()
print(f"Total species: {len(species_data)}")
print(f"Species with all genes: {sum(1 for c in species_gene_counts.values() if c == len(genes))}")
print(f"Species with ≥50% genes: {sum(1 for c in species_gene_counts.values() if c >= len(genes)/2)}")

# Filter to ensure at least one species per genus, prioritizing data completeness
min_genes_required = len(genes) / 2

# Group species by genus
genus_species = defaultdict(list)
for species in species_data.keys():
    genus = species.split('_')[0]
    gene_count = len(species_data[species])
    genus_species[genus].append((species, gene_count))

# Select best representative per genus (most complete data)
filtered_species = {}
for genus, species_list in genus_species.items():
    # Sort by gene count (descending), then by species name
    species_list.sort(key=lambda x: (-x[1], x[0]))
    
    # Always include the best species from each genus
    best_species, best_count = species_list[0]
    filtered_species[best_species] = species_data[best_species]
    
    # Also include other species from genus if they meet the min_genes threshold
    for species, gene_count in species_list[1:]:
        if gene_count >= min_genes_required:
            filtered_species[species] = species_data[species]

print()
print(f"Filtering strategy: At least one species per genus + species with ≥{min_genes_required:.1f} genes")
print(f"Number of genera: {len(genus_species)}")
print(f"Species before filtering: {len(species_data)}")
print(f"Species after filtering: {len(filtered_species)}")
print(f"Excluded {len(species_data) - len(filtered_species)} species")

# Create concatenated sequences
print()
print("===================================================================")
print("Creating concatenated sequences")
print("===================================================================")

concatenated_records = []

for species in sorted(filtered_species.keys()):
    # Build concatenated sequence for this species
    concat_seq = ""
    concat_parts = []
    
    for gene in genes:
        if gene in filtered_species[species]:
            # Select the first (longest or best) sequence for this gene
            # You could implement more sophisticated selection here
            vouchers = filtered_species[species][gene]
            if len(vouchers) > 1:
                # If multiple sequences, pick the one with fewest gaps
                best_voucher, best_seq = min(vouchers, key=lambda x: x[1].count('-') + x[1].count('N'))
                concat_parts.append(f"{gene}:{best_voucher}")
            else:
                best_voucher, best_seq = vouchers[0]
                concat_parts.append(f"{gene}:{best_voucher}")
            
            concat_seq += best_seq
        else:
            # Gene missing for this species, add gaps
            concat_seq += "-" * gene_lengths[gene]
            concat_parts.append(f"{gene}:missing")
    
    # Create record
    record_id = species
    # Create detailed header with accession information for each gene
    gene_accessions = ' '.join(concat_parts)
    record_description = f"genes={len([p for p in concat_parts if 'missing' not in p])}/{len(genes)} {gene_accessions}"
    
    record = SeqRecord(
        Seq(concat_seq),
        id=record_id,
        description=record_description
    )
    concatenated_records.append(record)
    
    print(f"{species:30s}: {len(concat_seq)} bp ({len([p for p in concat_parts if 'missing' not in p])}/{len(genes)} genes)")
    print(f"  Accessions: {gene_accessions}")

# Write concatenated sequences
output_file = f"{output_dir}/concatenated_all_genes.fasta"
SeqIO.write(concatenated_records, output_file, "fasta")

print()
print(f"Concatenated sequences written to: {output_file}")
print(f"Total sequences: {len(concatenated_records)}")
print(f"Concatenated length: {len(concatenated_records[0].seq)} bp")

# Write partition file for IQ-TREE
partition_file = f"{output_dir}/partitions.txt"
with open(partition_file, 'w') as f:
    start = 1
    for gene in genes:
        end = start + gene_lengths[gene] - 1
        f.write(f"DNA, {gene} = {start}-{end}\n")
        start = end + 1

print(f"Partition file written to: {partition_file}")

PYTHON_SCRIPT

echo ""
echo "==================================================================="
echo "Concatenation completed"
echo "Completed: $(date)"
echo "==================================================================="
