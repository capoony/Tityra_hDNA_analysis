#!/bin/bash
###############################################################################
# Align complete mitochondrial genomes
# Author: Generated for Tityra project
# Description: Two alignment approaches:
#              1) Annotation-aware: Extract and align genes separately
#              2) Whole-genome: Align complete sequences with MAFFT
###############################################################################

WD="/media/inter/mkapun/projects/Tityra"

# Create output directories
mkdir -p "${WD}/results/complete_mitogenome_alignment"
mkdir -p "${WD}/results/complete_mitogenome_alignment/genes"
mkdir -p "${WD}/results/complete_mitogenome_alignment/whole_genome"

###############################################################################
# APPROACH 1: Annotation-aware alignment
###############################################################################
echo "======================================================================="
echo "APPROACH 1: Annotation-aware alignment (gene-by-gene)"
echo "======================================================================="

echo "Extracting genes from GenBank files..."

# Python script to extract all genes from GenBank files
python3 << 'PYEOF'
from Bio import SeqIO
import os
from collections import defaultdict

wd = "/media/inter/mkapun/projects/Tityra"
gb_dir = f"{wd}/data/refseq_genbank"
out_dir = f"{wd}/results/complete_mitogenome_alignment/genes"

# Dictionary to store genes: gene_name -> {species: sequence}
genes = defaultdict(dict)

# Standard mitochondrial gene names to look for
target_genes = [
    'COX1', 'COI', 'COX2', 'CO2', 'COII', 'COX3', 'CO3', 'COIII',
    'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6',
    'CYTB', 'CYB', 'ATP6', 'ATP8', 'NADH1', 'NADH2', 'NADH3',
    'NADH4', 'NADH4L', 'NADH5', 'NADH6'
]

# Normalize gene names
gene_aliases = {
    'COI': 'COX1', 'CO1': 'COX1',
    'COII': 'COX2', 'CO2': 'COX2',
    'COIII': 'COX3', 'CO3': 'COX3',
    'CYB': 'CYTB',
    'NADH1': 'ND1', 'NADH2': 'ND2', 'NADH3': 'ND3',
    'NADH4': 'ND4', 'NADH4L': 'ND4L', 'NADH5': 'ND5', 'NADH6': 'ND6'
}

print("Processing GenBank files...")
for gb_file in sorted(os.listdir(gb_dir)):
    if not gb_file.endswith('.gb'):
        continue
    
    species = gb_file.replace('.gb', '')
    print(f"  {species}")
    
    record = SeqIO.read(f"{gb_dir}/{gb_file}", "genbank")
    
    # Extract protein-coding genes
    for feature in record.features:
        if feature.type == 'CDS' or feature.type == 'gene':
            # Get gene name
            gene_name = None
            if 'gene' in feature.qualifiers:
                gene_name = feature.qualifiers['gene'][0].upper()
            elif 'product' in feature.qualifiers:
                product = feature.qualifiers['product'][0].upper()
                # Try to extract gene name from product
                for tg in target_genes:
                    if tg in product:
                        gene_name = tg
                        break
            
            if gene_name:
                # Normalize gene name
                gene_name = gene_aliases.get(gene_name, gene_name)
                
                # Extract sequence
                try:
                    if feature.type == 'CDS':
                        seq = feature.extract(record.seq)
                    else:
                        seq = feature.location.extract(record.seq)
                    
                    if len(seq) > 0:
                        genes[gene_name][species] = str(seq)
                except Exception as e:
                    pass

print(f"\nTotal genes found: {len(genes)}")

# Write each gene to a separate FASTA file
for gene_name, species_seqs in genes.items():
    if len(species_seqs) >= 2:  # Only keep genes present in at least 2 species
        fasta_file = f"{out_dir}/{gene_name}.fasta"
        with open(fasta_file, 'w') as f:
            for species, seq in sorted(species_seqs.items()):
                f.write(f">{species}\n")
                # Write sequence in 80 character lines
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + '\n')
        print(f"Wrote {gene_name}.fasta ({len(species_seqs)} species)")

print("\nGene extraction complete!")
PYEOF

if [ $? -ne 0 ]; then
    echo "Error: Gene extraction failed"
    exit 1
fi

# Count extracted genes
n_genes=$(ls -1 "${WD}/results/complete_mitogenome_alignment/genes"/*.fasta 2>/dev/null | wc -l)
echo ""
echo "Extracted ${n_genes} genes"

if [ "$n_genes" -gt 0 ]; then
    echo ""
    echo "Aligning each gene with MAFFT..."

    MAFFT="/opt/anaconda3/envs/mafft-7.487/bin/mafft"

    # Align each gene separately
    for gene_file in "${WD}/results/complete_mitogenome_alignment/genes"/*.fasta; do
        gene_name=$(basename "$gene_file" .fasta)
        aligned_file="${gene_file%.fasta}_aligned.fasta"
        
        echo "  Aligning ${gene_name}..."
        "$MAFFT" --auto --thread 4 "$gene_file" > "$aligned_file" 2>/dev/null
    done

    echo ""
    echo "Concatenating aligned genes..."

    # Python script to concatenate alignments and create partitions file
    python3 << 'PYEOF'
from Bio import SeqIO
from collections import defaultdict
import os

wd = "/media/inter/mkapun/projects/Tityra"
gene_dir = f"{wd}/results/complete_mitogenome_alignment/genes"
out_file = f"{wd}/results/complete_mitogenome_alignment/annotation_aware_aligned.fasta"
partitions_file = f"{wd}/results/complete_mitogenome_alignment/partitions.txt"
partitions_nexus = f"{wd}/results/complete_mitogenome_alignment/partitions.nex"

# Dictionary to store concatenated sequences
concat_seqs = defaultdict(str)
gene_order = []
gene_lengths = []
species_list = set()

# Read all aligned gene files
for gene_file in sorted(os.listdir(gene_dir)):
    if not gene_file.endswith('_aligned.fasta'):
        continue
    
    gene_name = gene_file.replace('_aligned.fasta', '')
    gene_order.append(gene_name)
    
    alignment_file = f"{gene_dir}/{gene_file}"
    
    # Read alignment
    gene_seqs = {}
    for record in SeqIO.parse(alignment_file, "fasta"):
        gene_seqs[record.id] = str(record.seq)
        species_list.add(record.id)
    
    # Add sequences to concatenation (with gaps for missing species)
    alignment_length = len(list(gene_seqs.values())[0]) if gene_seqs else 0
    gene_lengths.append(alignment_length)
    
    for species in species_list:
        if species in gene_seqs:
            concat_seqs[species] += gene_seqs[species]
        else:
            # Add gaps for missing gene
            concat_seqs[species] += '-' * alignment_length
    
    print(f"Added {gene_name}: {alignment_length} bp")

# Write concatenated alignment
with open(out_file, 'w') as f:
    for species in sorted(concat_seqs.keys()):
        f.write(f">{species}\n")
        seq = concat_seqs[species]
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + '\n')

# Create partitions file for IQ-TREE (RAxML format)
# Format: DNA, GENE_NAME = start-end
position = 1
with open(partitions_file, 'w') as f:
    for gene_name, length in zip(gene_order, gene_lengths):
        if length > 0:
            end_pos = position + length - 1
            f.write(f"DNA, {gene_name} = {position}-{end_pos}\n")
            position = end_pos + 1

# Create partitions file in NEXUS format (alternative)
with open(partitions_nexus, 'w') as f:
    f.write("#nexus\n")
    f.write("begin sets;\n")
    position = 1
    for gene_name, length in zip(gene_order, gene_lengths):
        if length > 0:
            end_pos = position + length - 1
            f.write(f"  charset {gene_name} = {position}-{end_pos};\n")
            position = end_pos + 1
    f.write("end;\n")

print(f"\nAnnotation-aware alignment: {out_file}")
print(f"  Species: {len(concat_seqs)}")
print(f"  Genes: {len(gene_order)}")
print(f"  Total length: {len(list(concat_seqs.values())[0])} bp")
print(f"  Gene order: {', '.join(gene_order)}")
print(f"\nPartitions file created:")
print(f"  RAxML format: {partitions_file}")
print(f"  NEXUS format: {partitions_nexus}")

# Print partition summary
print(f"\nPartition summary:")
position = 1
for gene_name, length in zip(gene_order, gene_lengths):
    if length > 0:
        end_pos = position + length - 1
        print(f"  {gene_name:20s}: {position:6d}-{end_pos:6d} ({length:5d} bp)")
        position = end_pos + 1
PYEOF
fi

###############################################################################
# APPROACH 2: Whole-genome alignment
###############################################################################
echo ""
echo "======================================================================="
echo "APPROACH 2: Whole-genome alignment (MAFFT)"
echo "======================================================================="

echo "Extracting complete sequences from GenBank files..."

combined_fasta="${WD}/results/complete_mitogenome_alignment/whole_genome/all_mitogenomes.fasta"
> "${combined_fasta}"

# Extract complete sequences
for gb_file in "${WD}/data/refseq_genbank"/*.gb; do
    if [ -f "$gb_file" ]; then
        species=$(basename "$gb_file" .gb)
        echo "  Processing: ${species}"
        
        python3 << EOF
from Bio import SeqIO
try:
    record = SeqIO.read("${gb_file}", "genbank")
    with open("${combined_fasta}", "a") as f:
        f.write(f">${species}\n")
        seq_str = str(record.seq)
        for i in range(0, len(seq_str), 80):
            f.write(seq_str[i:i+80] + "\n")
    print(f"  Added: ${species} ({len(record.seq)} bp)")
except Exception as e:
    print(f"  Error: {e}")
EOF
    fi
done

n_sequences=$(grep -c "^>" "${combined_fasta}")
echo ""
echo "Total sequences: ${n_sequences}"

if [ "$n_sequences" -ge 2 ]; then
    echo ""
    echo "Running MAFFT alignment on complete genomes..."
    
    MAFFT="/opt/anaconda3/envs/mafft-7.487/bin/mafft"

    "$MAFFT" --auto --adjustdirection --thread 20 \
        "${combined_fasta}" \
        > "${WD}/results/complete_mitogenome_alignment/whole_genome_aligned.fasta"

    echo "Whole-genome alignment complete!"
fi

###############################################################################
# Summary
###############################################################################
echo ""
echo "======================================================================="
echo "Alignment Summary"
echo "======================================================================="

if [ -f "${WD}/results/complete_mitogenome_alignment/annotation_aware_aligned.fasta" ]; then
    n_seqs=$(grep -c "^>" "${WD}/results/complete_mitogenome_alignment/annotation_aware_aligned.fasta")
    length=$(head -2 "${WD}/results/complete_mitogenome_alignment/annotation_aware_aligned.fasta" | tail -1 | tr -d '\n' | wc -c)
    echo "Annotation-aware alignment:"
    echo "  File: annotation_aware_aligned.fasta"
    echo "  Species: ${n_seqs}"
    echo "  Length: ${length} bp"
fi

if [ -f "${WD}/results/complete_mitogenome_alignment/whole_genome_aligned.fasta" ]; then
    n_seqs=$(grep -c "^>" "${WD}/results/complete_mitogenome_alignment/whole_genome_aligned.fasta")
    length=$(head -2 "${WD}/results/complete_mitogenome_alignment/whole_genome_aligned.fasta" | tail -1 | tr -d '\n' | wc -c)
    echo ""
    echo "Whole-genome alignment:"
    echo "  File: whole_genome_aligned.fasta"
    echo "  Species: ${n_seqs}"
    echo "  Length: ${length} bp"
fi

echo ""
echo "Both alignments complete!"
