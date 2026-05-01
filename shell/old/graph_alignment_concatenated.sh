#!/bin/bash

###############################################################################
# Concatenated Graph-Based Variant Calling
# Maps reads once to a mega-graph of all BUSCO genes
# Much faster than individual mapping
###############################################################################

set -euo pipefail

WD=${1:-/media/inter/mkapun/projects/Tityra}
THREADS=${2:-50}
MIN_DEPTH=5  # Minimum average depth to include gene

echo "==================================================================="
echo "Concatenated Graph-Based Variant Calling"
echo "Working directory: ${WD}"
echo "Threads: ${THREADS}"
echo "Minimum depth threshold: ${MIN_DEPTH}x"
echo "==================================================================="

# Setup paths
VG="${WD}/results/graph_alignment/vg"
CONDA_ENV=/media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs

BCFTOOLS=${CONDA_ENV}/bin/bcftools
BGZIP=${CONDA_ENV}/bin/bgzip
TABIX=${CONDA_ENV}/bin/tabix
SAMTOOLS=${CONDA_ENV}/bin/samtools

PE1=${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz
PE2=${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz
MERGED=${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz

mkdir -p ${WD}/results/graph_alignment/concatenated

###############################################################################
# Step 1: Concatenate all aligned BUSCO genes into one MSA
###############################################################################

echo "Step 1: Concatenating all aligned BUSCO genes..."

if [ ! -f "${WD}/results/graph_alignment/concatenated/all_genes.fasta" ]; then
    cat ${WD}/results/graph_alignment/aligned/*.fasta > ${WD}/results/graph_alignment/concatenated/all_genes.fasta
    num_seqs=$(grep -c "^>" ${WD}/results/graph_alignment/concatenated/all_genes.fasta || true)
    echo "Concatenated MSA contains ${num_seqs} sequences"
else
    echo "Concatenated MSA already exists"
fi

###############################################################################
# Step 2: Build one large variation graph
###############################################################################

echo "Step 2: Building mega variation graph..."

if [ ! -f "${WD}/results/graph_alignment/concatenated/mega.xg" ]; then
    echo "  Constructing graph from concatenated MSA..."
    ${VG} construct \
        -M ${WD}/results/graph_alignment/concatenated/all_genes.fasta \
        -m 32 \
        > ${WD}/results/graph_alignment/concatenated/mega.vg 2>/dev/null
    
    echo "  Indexing graph..."
    ${VG} index \
        -x ${WD}/results/graph_alignment/concatenated/mega.xg \
        ${WD}/results/graph_alignment/concatenated/mega.vg 2>/dev/null
    
    echo "  Pruning complex regions..."
    ${VG} prune \
        ${WD}/results/graph_alignment/concatenated/mega.vg \
        > ${WD}/results/graph_alignment/concatenated/mega.pruned.vg 2>/dev/null
    
    echo "  Building GCSA index (this may take a while)..."
    ${VG} index \
        -g ${WD}/results/graph_alignment/concatenated/mega.gcsa \
        -k 16 \
        ${WD}/results/graph_alignment/concatenated/mega.pruned.vg 2>/dev/null
    
    rm -f ${WD}/results/graph_alignment/concatenated/mega.vg
    rm -f ${WD}/results/graph_alignment/concatenated/mega.pruned.vg
    
    echo "  Graph construction complete!"
else
    echo "  Mega graph already exists, skipping..."
fi

###############################################################################
# Step 3: Map reads ONCE to mega-graph
###############################################################################

echo "Step 3: Mapping reads to mega-graph..."

if [ ! -f "${WD}/results/graph_alignment/concatenated/mega.gam" ]; then
    echo "  Mapping PE reads..."
    ${VG} map \
        -x ${WD}/results/graph_alignment/concatenated/mega.xg \
        -g ${WD}/results/graph_alignment/concatenated/mega.gcsa \
        -f ${PE1} -f ${PE2} \
        -t ${THREADS} \
        > ${WD}/results/graph_alignment/concatenated/mega.gam 2>/dev/null
    
    echo "  Mapping merged reads (append)..."
    ${VG} map \
        -x ${WD}/results/graph_alignment/concatenated/mega.xg \
        -g ${WD}/results/graph_alignment/concatenated/mega.gcsa \
        -f ${MERGED} \
        -t ${THREADS} \
        >> ${WD}/results/graph_alignment/concatenated/mega.gam 2>/dev/null
    
    echo "  Mapping complete!"
    echo "  GAM file size: $(du -h ${WD}/results/graph_alignment/concatenated/mega.gam | cut -f1)"
else
    echo "  Mega GAM already exists, skipping..."
fi

###############################################################################
# Step 4: Pack and call variants on entire graph
###############################################################################

echo "Step 4: Calling variants on mega-graph..."

if [ ! -f "${WD}/results/graph_alignment/concatenated/mega.vcf.gz" ]; then
    echo "  Packing alignments..."
    ${VG} pack \
        -x ${WD}/results/graph_alignment/concatenated/mega.xg \
        -g ${WD}/results/graph_alignment/concatenated/mega.gam \
        -o ${WD}/results/graph_alignment/concatenated/mega.pack \
        2>/dev/null
    
    echo "  Pack file size: $(du -h ${WD}/results/graph_alignment/concatenated/mega.pack | cut -f1)"
    
    echo "  Calling variants (relaxed thresholds for ancient DNA)..."
    ${VG} call \
        ${WD}/results/graph_alignment/concatenated/mega.xg \
        -k ${WD}/results/graph_alignment/concatenated/mega.pack \
        -d 1 -s 1 \
        > ${WD}/results/graph_alignment/concatenated/mega_raw.vcf \
        2>/dev/null
    
    echo "  Filtering variants (DP > 5)..."
    ${BCFTOOLS} view \
        -i 'FORMAT/DP>5' \
        ${WD}/results/graph_alignment/concatenated/mega_raw.vcf \
        -o ${WD}/results/graph_alignment/concatenated/mega_filtered.vcf \
        2>/dev/null
    
    echo "  Compressing and indexing VCF..."
    ${BGZIP} -c ${WD}/results/graph_alignment/concatenated/mega_filtered.vcf \
        > ${WD}/results/graph_alignment/concatenated/mega.vcf.gz 2>/dev/null
    ${TABIX} -p vcf ${WD}/results/graph_alignment/concatenated/mega.vcf.gz 2>/dev/null
    
    rm -f ${WD}/results/graph_alignment/concatenated/mega_raw.vcf
    rm -f ${WD}/results/graph_alignment/concatenated/mega_filtered.vcf
    
    total_variants=$(${BCFTOOLS} view -H ${WD}/results/graph_alignment/concatenated/mega.vcf.gz | wc -l)
    echo "  Total variants called (DP>5): ${total_variants}"
else
    echo "  Mega VCF already exists, skipping..."
fi

###############################################################################
# Step 5: Calculate depth statistics per gene/path
###############################################################################

echo "Step 5: Calculating depth statistics per gene..."

# Get list of paths from VCF
${BCFTOOLS} view -h ${WD}/results/graph_alignment/concatenated/mega.vcf.gz | \
    grep "^##contig" | \
    sed 's/.*ID=\([^,]*\).*/\1/' \
    > ${WD}/results/graph_alignment/concatenated/path_list.txt

# Calculate average depth per path using VCF
echo "gene_id,path_name,variant_count,avg_depth,max_depth" > ${WD}/results/graph_alignment/concatenated/depth_stats.csv

while read -r path; do
    # Extract gene ID
    gene=$(echo "${path}" | sed 's/_ref[0-9]*_Pachyramphus_minor//;s/_ref[0-9]*_Oxyruncus_cristatus//;s/_ref[0-9]*_Tyrannus_savana//')
    
    # Extract variants for this path and calculate depth stats
    depth_stats=$(${BCFTOOLS} view -r "${path}" -H ${WD}/results/graph_alignment/concatenated/mega.vcf.gz | \
        awk -F'\t' '{
            split($10, a, ":");
            depth = a[3];
            sum += depth;
            count++;
            if (depth > max) max = depth;
        }
        END {
            if (count > 0) {
                printf "%d,%.2f,%d", count, sum/count, max;
            } else {
                printf "0,0,0";
            }
        }' 2>/dev/null)
    
    echo "${gene},${path},${depth_stats}" >> ${WD}/results/graph_alignment/concatenated/depth_stats.csv
done < ${WD}/results/graph_alignment/concatenated/path_list.txt

echo "  Depth statistics calculated"

# Create list of genes passing depth threshold
awk -F',' -v min_depth="${MIN_DEPTH}" '$4 >= min_depth {print $1}' \
    ${WD}/results/graph_alignment/concatenated/depth_stats.csv | \
    sort -u > ${WD}/results/graph_alignment/concatenated/genes_pass_depth.txt

passed=$(wc -l < ${WD}/results/graph_alignment/concatenated/genes_pass_depth.txt)
total=$(tail -n +2 ${WD}/results/graph_alignment/concatenated/depth_stats.csv | awk -F',' '{print $1}' | sort -u | wc -l)
echo "  Genes passing depth threshold (≥${MIN_DEPTH}x): ${passed}/${total}"

###############################################################################
# Step 6: Split VCF by gene and generate consensus sequences
###############################################################################

echo "Step 6: Generating consensus for high-depth genes only..."

processed=0
passed_genes=0
low_depth=0

while read -r gene; do
    # Skip if consensus already exists
    if [ -f "${WD}/results/graph_alignment/consensus/${gene}_Tityra.fasta" ]; then
        ((processed++))
        ((passed_genes++))
        continue
    fi
    
    if [ $((processed % 100)) -eq 0 ]; then
        echo "  Processing ${processed}: ${gene}"
    fi
    
    # Find path(s) for this gene - use primary reference (Pachyramphus_minor)
    path=$(grep "^${gene}_ref.*Pachyramphus_minor" ${WD}/results/graph_alignment/concatenated/path_list.txt | head -n 1)
    
    if [ -z "${path}" ]; then
        # Try without ref prefix
        path=$(grep "Pachyramphus_minor" ${WD}/results/graph_alignment/concatenated/path_list.txt | grep "${gene}" | head -n 1)
    fi
    
    if [ -z "${path}" ]; then
        echo "    Warning: No path found for ${gene}"
        ((processed++))
        continue
    fi
    
    # Extract variants for this path/gene
    ${BCFTOOLS} view \
        -r "${path}" \
        ${WD}/results/graph_alignment/concatenated/mega.vcf.gz \
        -o ${WD}/results/graph_alignment/concatenated/${gene}.vcf \
        2>/dev/null || {
        echo "    Warning: Failed to extract variants for ${gene}"
        ((processed++))
        continue
    }
    
    # Count variants
    variant_count=$(grep -v "^#" ${WD}/results/graph_alignment/concatenated/${gene}.vcf | wc -l || echo "0")
    
    if [ "${variant_count}" -eq 0 ]; then
        # No variants - just copy reference sequence
        ${SAMTOOLS} faidx \
            ${WD}/results/graph_alignment/busco_regions/${gene}.fasta \
            ref1_Pachyramphus_minor 2>/dev/null | \
            sed "s/>ref1_Pachyramphus_minor/>${gene}|Tityra_leucura_graph/" \
            > ${WD}/results/graph_alignment/consensus/${gene}_Tityra.fasta
    else
        # Apply variants to generate consensus
        ${BGZIP} -c ${WD}/results/graph_alignment/concatenated/${gene}.vcf \
            > ${WD}/results/graph_alignment/concatenated/${gene}.vcf.gz 2>/dev/null
        ${TABIX} -p vcf ${WD}/results/graph_alignment/concatenated/${gene}.vcf.gz 2>/dev/null
        
        # Extract reference sequence
        ${SAMTOOLS} faidx \
            ${WD}/results/graph_alignment/busco_regions/${gene}.fasta \
            ref1_Pachyramphus_minor \
            > ${WD}/results/graph_alignment/concatenated/${gene}_ref.fasta 2>/dev/null
        
        # Generate consensus
        ${BCFTOOLS} consensus \
            -f ${WD}/results/graph_alignment/concatenated/${gene}_ref.fasta \
            ${WD}/results/graph_alignment/concatenated/${gene}.vcf.gz \
            -o ${WD}/results/graph_alignment/concatenated/${gene}_consensus.fasta \
            2>/dev/null || {
            echo "    Warning: Consensus failed for ${gene}"
            ((processed++))
            rm -f ${WD}/results/graph_alignment/concatenated/${gene}*
            continue
        }
        
        # Rename header
        sed "s/>ref1_Pachyramphus_minor/>${gene}|Tityra_leucura_graph/" \
            ${WD}/results/graph_alignment/concatenated/${gene}_consensus.fasta \
            > ${WD}/results/graph_alignment/consensus/${gene}_Tityra.fasta
    fi
    
    # Cleanup
    rm -f ${WD}/results/graph_alignment/concatenated/${gene}*
    
    ((processed++))
    ((passed_genes++))
done < ${WD}/results/graph_alignment/concatenated/genes_pass_depth.txt

echo "Consensus generation complete:"
echo "  Genes passing depth filter: ${passed_genes}"
echo "  Genes filtered out (low depth): $((total - passed_genes))"

###############################################################################
# Step 7: Translate to proteins (only high-depth genes)
###############################################################################

echo "Step 7: Translating consensus sequences to proteins..."

mkdir -p ${WD}/results/phylogeny_graph/busco_aa

python3 << 'PYTHON_TRANSLATE'
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq

wd = os.environ['WD']
consensus_dir = f"{wd}/results/graph_alignment/consensus"
output_dir = f"{wd}/results/phylogeny_graph/busco_aa"

# Read list of genes passing depth filter
genes_pass_file = f"{wd}/results/graph_alignment/concatenated/genes_pass_depth.txt"
with open(genes_pass_file) as f:
    genes_pass = set(line.strip() for line in f)

translated = 0
skipped = 0

for consensus_file in os.listdir(consensus_dir):
    if not consensus_file.endswith('_Tityra.fasta'):
        continue
    
    gene = consensus_file.replace('_Tityra.fasta', '')
    
    # Skip genes that didn't pass depth filter
    if gene not in genes_pass:
        skipped += 1
        continue
    
    consensus_path = f"{consensus_dir}/{consensus_file}"
    
    try:
        for record in SeqIO.parse(consensus_path, "fasta"):
            # Translate in frame 1
            protein_seq = record.seq.translate()
            
            # Write protein sequence
            with open(f"{output_dir}/{gene}.faa", 'w') as out:
                out.write(f">{gene}|Tityra_leucura_graph\n")
                out.write(f"{str(protein_seq)}\n")
            
            translated += 1
            break
    except Exception as e:
        print(f"Warning: Could not translate {gene}: {e}", file=sys.stderr)
        continue

print(f"Translated {translated} high-depth genes to protein", file=sys.stderr)
print(f"Skipped {skipped} low-depth genes", file=sys.stderr)
PYTHON_TRANSLATE

echo "==================================================================="
echo "Concatenated graph-based variant calling complete!"
echo "High-depth consensus sequences (≥${MIN_DEPTH}x): $(wc -l < ${WD}/results/graph_alignment/concatenated/genes_pass_depth.txt)"
echo "Protein sequences: $(ls ${WD}/results/phylogeny_graph/busco_aa/*.faa 2>/dev/null | wc -l)"
echo "Depth statistics: ${WD}/results/graph_alignment/concatenated/depth_stats.csv"
echo "==================================================================="
