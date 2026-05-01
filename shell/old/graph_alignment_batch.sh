#!/bin/bash

###############################################################################
# Space-Efficient Batch Graph-Based Variant Calling
# Processes genes in batches of 100, cleaning up intermediate files immediately
###############################################################################

set -euo pipefail

WD=${1:-/media/inter/mkapun/projects/Tityra}
THREADS=${2:-50}
BATCH_SIZE=100  # Process 100 genes at a time
MIN_DEPTH=5     # Minimum average depth to include gene

echo "==================================================================="
echo "Space-Efficient Batch Graph-Based Variant Calling"
echo "Working directory: ${WD}"
echo "Threads: ${THREADS}"
echo "Batch size: ${BATCH_SIZE} genes"
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

BATCH_DIR=${WD}/results/graph_alignment/batch_temp
mkdir -p ${BATCH_DIR}
mkdir -p ${WD}/results/graph_alignment/consensus

# Get all genes
cd ${WD}/results/graph_alignment/aligned
GENES=(*.fasta)
TOTAL_GENES=${#GENES[@]}

echo "Total genes to process: ${TOTAL_GENES}"

###############################################################################
# Process in batches
###############################################################################

for ((i=0; i<${TOTAL_GENES}; i+=${BATCH_SIZE})); do
    BATCH_NUM=$((i / BATCH_SIZE + 1))
    END=$((i + BATCH_SIZE - 1))
    if [ ${END} -ge ${TOTAL_GENES} ]; then
        END=$((TOTAL_GENES - 1))
    fi
    
    BATCH_GENES=("${GENES[@]:${i}:${BATCH_SIZE}}")
    ACTUAL_BATCH_SIZE=${#BATCH_GENES[@]}
    
    echo ""
    echo "==================================================================="
    echo "Processing Batch ${BATCH_NUM}: Genes $((i+1)) to $((i+ACTUAL_BATCH_SIZE))"
    echo "==================================================================="
    
    # Clean batch directory
    rm -rf ${BATCH_DIR}/*
    
    ###########################################################################
    # Step 1: Concatenate this batch
    ###########################################################################
    
    echo "Concatenating ${ACTUAL_BATCH_SIZE} genes..."
    for gene in "${BATCH_GENES[@]}"; do
        cat ${WD}/results/graph_alignment/aligned/${gene} >> ${BATCH_DIR}/batch.fasta
    done
    
    ###########################################################################
    # Step 2: Build batch graph
    ###########################################################################
    
    echo "Building variation graph for batch..."
    ${VG} construct -M ${BATCH_DIR}/batch.fasta -m 32 -t ${THREADS} > ${BATCH_DIR}/batch.vg
    
    echo "Indexing XG..."
    ${VG} index -x ${BATCH_DIR}/batch.xg ${BATCH_DIR}/batch.vg -t ${THREADS}
    rm ${BATCH_DIR}/batch.vg  # Clean up immediately
    
    echo "Pruning graph..."
    ${VG} prune ${BATCH_DIR}/batch.xg -k 16 -t ${THREADS} > ${BATCH_DIR}/batch.pruned.vg
    
    echo "Building GCSA index..."
    ${VG} index -g ${BATCH_DIR}/batch.gcsa -k 16 ${BATCH_DIR}/batch.pruned.vg -t ${THREADS}
    rm ${BATCH_DIR}/batch.pruned.vg  # Clean up
    
    ###########################################################################
    # Step 3: Map reads to batch graph
    ###########################################################################
    
    echo "Mapping paired-end reads..."
    ${VG} map -x ${BATCH_DIR}/batch.xg -g ${BATCH_DIR}/batch.gcsa \
        -f ${PE1} -f ${PE2} -t ${THREADS} > ${BATCH_DIR}/batch.gam
    
    echo "Appending merged reads..."
    ${VG} map -x ${BATCH_DIR}/batch.xg -g ${BATCH_DIR}/batch.gcsa \
        -f ${MERGED} -t ${THREADS} >> ${BATCH_DIR}/batch.gam
    
    # Clean up indices
    rm -f ${BATCH_DIR}/batch.gcsa*
    
    ###########################################################################
    # Step 4: Call variants on batch
    ###########################################################################
    
    echo "Packing read coverage..."
    ${VG} pack -x ${BATCH_DIR}/batch.xg -g ${BATCH_DIR}/batch.gam \
        -o ${BATCH_DIR}/batch.pack -t ${THREADS}
    
    rm ${BATCH_DIR}/batch.gam  # Clean up GAM immediately
    
    echo "Calling variants..."
    ${VG} call ${BATCH_DIR}/batch.xg -k ${BATCH_DIR}/batch.pack \
        -d 1 -s 1 -t ${THREADS} > ${BATCH_DIR}/batch_raw.vcf
    
    rm ${BATCH_DIR}/batch.pack  # Clean up pack file
    
    echo "Filtering variants (DP>${MIN_DEPTH})..."
    ${BCFTOOLS} view -i "FORMAT/DP>${MIN_DEPTH}" ${BATCH_DIR}/batch_raw.vcf \
        -o ${BATCH_DIR}/batch_filtered.vcf
    
    rm ${BATCH_DIR}/batch_raw.vcf  # Clean up
    
    echo "Compressing and indexing VCF..."
    ${BGZIP} -c ${BATCH_DIR}/batch_filtered.vcf > ${BATCH_DIR}/batch.vcf.gz
    ${TABIX} -p vcf ${BATCH_DIR}/batch.vcf.gz
    
    rm ${BATCH_DIR}/batch_filtered.vcf  # Clean up
    
    ###########################################################################
    # Step 5: Split by gene and generate consensus
    ###########################################################################
    
    echo "Extracting paths from VCF..."
    PATHS=$(${BCFTOOLS} query -l ${BATCH_DIR}/batch.vcf.gz | grep "ref1_Pachyramphus_minor" || true)
    
    for path in ${PATHS}; do
        # Extract gene name from path
        GENE=$(echo ${path} | sed 's/^ref1_Pachyramphus_minor#0#//' | sed 's/:.*$//')
        GENE_FILE=$(echo ${GENE} | sed 's/:/_/g').fasta
        
        # Skip if already processed
        if [ -f "${WD}/results/graph_alignment/consensus/${GENE}_Tityra.fasta" ]; then
            continue
        fi
        
        echo "  Processing ${GENE}..."
        
        # Extract reference for this gene
        ${SAMTOOLS} faidx ${BATCH_DIR}/batch.xg "${path}" > ${BATCH_DIR}/${GENE}_ref.fasta 2>/dev/null || continue
        
        # Check if variants exist for this gene
        ${BCFTOOLS} view -r "${path}" ${BATCH_DIR}/batch.vcf.gz > ${BATCH_DIR}/${GENE}.vcf
        
        # Check average depth
        AVG_DEPTH=$(${BCFTOOLS} query -f '%DP\n' ${BATCH_DIR}/${GENE}.vcf | \
            awk '{sum+=$1; n++} END {if(n>0) print sum/n; else print 0}')
        
        if (( $(echo "${AVG_DEPTH} < ${MIN_DEPTH}" | bc -l) )); then
            echo "    Skipping ${GENE}: average depth ${AVG_DEPTH}x < ${MIN_DEPTH}x"
            rm -f ${BATCH_DIR}/${GENE}*
            continue
        fi
        
        # Compress and index gene VCF
        ${BGZIP} -c ${BATCH_DIR}/${GENE}.vcf > ${BATCH_DIR}/${GENE}.vcf.gz
        ${TABIX} -p vcf ${BATCH_DIR}/${GENE}.vcf.gz
        
        # Generate consensus
        ${BCFTOOLS} consensus -f ${BATCH_DIR}/${GENE}_ref.fasta \
            ${BATCH_DIR}/${GENE}.vcf.gz > ${BATCH_DIR}/${GENE}_consensus.fasta
        
        # Rename header
        sed "s/^>.*/>Tityra_leucura_graph/" ${BATCH_DIR}/${GENE}_consensus.fasta \
            > ${WD}/results/graph_alignment/consensus/${GENE}_Tityra.fasta
        
        echo "    Generated consensus for ${GENE} (depth: ${AVG_DEPTH}x)"
        
        # Clean up gene files
        rm -f ${BATCH_DIR}/${GENE}*
    done
    
    # Clean up batch files
    rm -rf ${BATCH_DIR}/*
    
    echo "Batch ${BATCH_NUM} complete"
done

###############################################################################
# Step 6: Translate to proteins
###############################################################################

echo ""
echo "==================================================================="
echo "Translating consensus sequences to proteins"
echo "==================================================================="

python3 << 'PYTHON_SCRIPT'
import os
from Bio import SeqIO
from Bio.Seq import Seq

wd = os.environ['WD']
consensus_dir = f"{wd}/results/graph_alignment/consensus"

if not os.path.exists(consensus_dir):
    print("No consensus sequences found")
    exit(0)

processed = 0
skipped = 0

for fasta_file in os.listdir(consensus_dir):
    if not fasta_file.endswith('_Tityra.fasta'):
        continue
    
    gene_name = fasta_file.replace('_Tityra.fasta', '')
    protein_file = f"{consensus_dir}/{gene_name}_Tityra_protein.fasta"
    
    if os.path.exists(protein_file):
        continue
    
    fasta_path = f"{consensus_dir}/{fasta_file}"
    
    try:
        record = SeqIO.read(fasta_path, "fasta")
        
        # Translate to protein
        protein_seq = record.seq.translate()
        
        # Check for stop codons in middle
        if '*' in str(protein_seq)[:-1]:
            print(f"Skipping {gene_name}: internal stop codon")
            skipped += 1
            continue
        
        # Write protein sequence
        with open(protein_file, 'w') as f:
            f.write(f">{record.id}\n")
            f.write(str(protein_seq) + '\n')
        
        processed += 1
        
    except Exception as e:
        print(f"Error processing {gene_name}: {e}")
        skipped += 1

print(f"Translation complete: {processed} genes translated, {skipped} skipped")
PYTHON_SCRIPT

###############################################################################
# Final summary
###############################################################################

echo ""
echo "==================================================================="
echo "Pipeline complete!"
echo "==================================================================="

CONSENSUS_COUNT=$(ls ${WD}/results/graph_alignment/consensus/*_Tityra.fasta 2>/dev/null | wc -l)
PROTEIN_COUNT=$(ls ${WD}/results/graph_alignment/consensus/*_protein.fasta 2>/dev/null | wc -l)

echo "Consensus sequences: ${CONSENSUS_COUNT}"
echo "Protein sequences: ${PROTEIN_COUNT}"
echo "Results in: ${WD}/results/graph_alignment/consensus/"
