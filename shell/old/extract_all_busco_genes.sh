#!/bin/bash

###############################################################################
# Extract all ~8100 BUSCO genes from the 3 references
###############################################################################

set -euo pipefail

WD=/media/inter/mkapun/projects/Tityra

echo "=================================================================="
echo "Extracting all BUSCO genes present in 2+ references"
echo "=================================================================="

mkdir -p ${WD}/results/graph_alignment/busco_regions

python3 << 'PYTHON_EXTRACT'
import os
import sys
from Bio import SeqIO

wd = os.environ.get('WD', '/media/inter/mkapun/projects/Tityra')
busco_dir = f"{wd}/results/graph_alignment/busco_refs"
output_dir = f"{wd}/results/graph_alignment/busco_regions"

# Find all complete BUSCO genes across references
busco_genes = set()
busco_refs = {}

print("Finding BUSCO references...", file=sys.stderr)
for busco_result in os.listdir(busco_dir):
    if not busco_result.endswith('_busco'):
        continue
    
    ref_name = busco_result.replace('_busco', '')
    seq_dir = f"{busco_dir}/{busco_result}/run_aves_odb10/busco_sequences/single_copy_busco_sequences"
    
    if not os.path.exists(seq_dir):
        print(f"Warning: No sequences found for {ref_name}", file=sys.stderr)
        continue
    
    busco_refs[ref_name] = seq_dir
    print(f"  Found: {ref_name}", file=sys.stderr)
    
    # Find all .fna (nucleotide coding sequence) files
    for fname in os.listdir(seq_dir):
        if fname.endswith('.fna'):
            gene_id = fname.replace('.fna', '')
            busco_genes.add(gene_id)

print(f"\nFound {len(busco_genes)} unique BUSCO genes across {len(busco_refs)} references", file=sys.stderr)

# For each BUSCO gene, collect sequences from all references
extracted = 0
skipped = 0

for i, gene in enumerate(sorted(busco_genes), 1):
    if i % 500 == 0:
        print(f"Progress: {i}/{len(busco_genes)} genes processed ({extracted} extracted, {skipped} skipped)", file=sys.stderr)
    
    gene_seqs = []
    
    for ref_name, seq_dir in busco_refs.items():
        fna_file = f"{seq_dir}/{gene}.fna"
        
        if not os.path.exists(fna_file):
            continue
        
        # Read the BUSCO-extracted coding sequence
        try:
            for record in SeqIO.parse(fna_file, "fasta"):
                gene_seqs.append(f">{ref_name}\n{str(record.seq)}\n")
                break  # Only take first sequence (should be only one)
        except Exception as e:
            print(f"Warning: Error parsing {gene} from {ref_name}: {e}", file=sys.stderr)
            continue
    
    # Write multi-FASTA only if we have sequences from at least 2 references
    if len(gene_seqs) >= 2:
        with open(f"{output_dir}/{gene}.fasta", 'w') as out:
            out.writelines(gene_seqs)
        extracted += 1
    else:
        skipped += 1

print(f"\n{'='*70}", file=sys.stderr)
print(f"BUSCO extraction complete!", file=sys.stderr)
print(f"  Extracted: {extracted} genes with >=2 reference sequences", file=sys.stderr)
print(f"  Skipped: {skipped} genes (only 1 reference)", file=sys.stderr)
print(f"  Output: {output_dir}", file=sys.stderr)
print(f"{'='*70}", file=sys.stderr)
PYTHON_EXTRACT

echo ""
echo "Final count: $(ls ${WD}/results/graph_alignment/busco_regions/*.fasta 2>/dev/null | wc -l) genes"
echo "Done!"
