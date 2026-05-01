#!/bin/bash

###############################################################################
# Extract RAG-1 and RAG-2 genes from Tityra leucura assembly using BLAST
# Description: Use reference RAG sequences to find and extract RAG genes
#              from the assembled genome
###############################################################################

set -euo pipefail

# Set working directory
WD=/media/inter/mkapun/projects/Tityra

# Input files
ASSEMBLY="${WD}/results/denovo/output/Tityra_ILL.fa.gz"
MITO_DIR="${WD}/data/Mito"

# Output directory
OUTPUT_DIR="${WD}/results/rag_genes_from_assembly"
mkdir -p "${OUTPUT_DIR}"

# Log file
LOG="${WD}/logs/extract_rag_genes.log"
mkdir -p "${WD}/logs"

exec > >(tee -a "${LOG}") 2>&1

echo "==================================================================="
echo "Extracting RAG genes from Tityra leucura assembly"
echo "Started: $(date)"
echo "==================================================================="

# Uncompress assembly if needed
ASSEMBLY_UNCOMPRESSED="${OUTPUT_DIR}/Tityra_ILL.fa"
if [ ! -f "${ASSEMBLY_UNCOMPRESSED}" ]; then
    echo "Uncompressing assembly..."
    gunzip -c "${ASSEMBLY}" > "${ASSEMBLY_UNCOMPRESSED}"
fi

# Process RAG-1 and RAG-2
for RAG_GENE in RAG-1 RAG-2; do
    echo ""
    echo "==================================================================="
    echo "Processing ${RAG_GENE}"
    echo "==================================================================="
    
    # Find a reference sequence (use a closely related species)
    # Try Schiffornis, Laniocera, or Pachyramphus as references
    REF_FILE=""
    for GENUS in Schiffornis Laniocera Pachyramphus Tityra; do
        POTENTIAL_REF="${MITO_DIR}/${RAG_GENE}_${GENUS}.fasta"
        if [ -f "${POTENTIAL_REF}" ] && [ -s "${POTENTIAL_REF}" ]; then
            REF_FILE="${POTENTIAL_REF}"
            echo "Using reference from ${GENUS}: ${REF_FILE}"
            break
        fi
    done
    
    if [ -z "${REF_FILE}" ]; then
        echo "Warning: No reference sequence found for ${RAG_GENE}, skipping..."
        continue
    fi
    
    # Create BLAST database from assembly
    echo "Creating BLAST database..."
    makeblastdb -in "${ASSEMBLY_UNCOMPRESSED}" \
                -dbtype nucl \
                -out "${OUTPUT_DIR}/Tityra_assembly_db" \
                -parse_seqids \
                > /dev/null 2>&1
    
    # Run BLAST search
    echo "Running BLAST search..."
    BLAST_OUTPUT="${OUTPUT_DIR}/${RAG_GENE}_blast_hits.txt"
    blastn -query "${REF_FILE}" \
           -db "${OUTPUT_DIR}/Tityra_assembly_db" \
           -out "${BLAST_OUTPUT}" \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
           -evalue 1e-10 \
           -max_target_seqs 5
    
    # Check if we found hits
    if [ ! -s "${BLAST_OUTPUT}" ]; then
        echo "No BLAST hits found for ${RAG_GENE}"
        continue
    fi
    
    echo "BLAST hits found:"
    cat "${BLAST_OUTPUT}"
    
    # Extract the best hit (highest bitscore)
    BEST_HIT=$(sort -k12,12nr "${BLAST_OUTPUT}" | head -1)
    
    if [ -z "${BEST_HIT}" ]; then
        echo "No valid hits found for ${RAG_GENE}"
        continue
    fi
    
    # Parse the best hit
    CONTIG=$(echo "${BEST_HIT}" | awk '{print $2}')
    START=$(echo "${BEST_HIT}" | awk '{print $9}')
    END=$(echo "${BEST_HIT}" | awk '{print $10}')
    PIDENT=$(echo "${BEST_HIT}" | awk '{print $3}')
    QCOV=$(echo "${BEST_HIT}" | awk '{print $13}')
    
    echo ""
    echo "Best hit for ${RAG_GENE}:"
    echo "  Contig: ${CONTIG}"
    echo "  Position: ${START}-${END}"
    echo "  Identity: ${PIDENT}%"
    echo "  Query coverage: ${QCOV}%"
    
    # Extract sequence using samtools faidx
    echo "Indexing assembly..."
    samtools faidx "${ASSEMBLY_UNCOMPRESSED}"
    
    # Determine strand and extract
    if [ ${START} -lt ${END} ]; then
        REGION="${CONTIG}:${START}-${END}"
        STRAND="forward"
    else
        REGION="${CONTIG}:${END}-${START}"
        STRAND="reverse"
    fi
    
    echo "Extracting sequence (${STRAND} strand)..."
    samtools faidx "${ASSEMBLY_UNCOMPRESSED}" "${REGION}" > "${OUTPUT_DIR}/${RAG_GENE}_extracted.fasta"
    
    # Reverse complement if needed
    if [ "${STRAND}" = "reverse" ]; then
        echo "Reverse complementing sequence..."
        seqkit seq -r -p "${OUTPUT_DIR}/${RAG_GENE}_extracted.fasta" > "${OUTPUT_DIR}/${RAG_GENE}_extracted_rc.fasta"
        mv "${OUTPUT_DIR}/${RAG_GENE}_extracted_rc.fasta" "${OUTPUT_DIR}/${RAG_GENE}_extracted.fasta"
    fi
    
    # Rename header to Tityra_leucura
    sed -i "s/>.*/>Tityra_leucura ${RAG_GENE} extracted from assembly/" "${OUTPUT_DIR}/${RAG_GENE}_extracted.fasta"
    
    # Create final output with proper naming
    RAG_NORMALIZED=$(echo "${RAG_GENE}" | sed 's/-//g')  # RAG-1 -> RAG1
    cp "${OUTPUT_DIR}/${RAG_GENE}_extracted.fasta" "${OUTPUT_DIR}/${RAG_NORMALIZED}_Tityra_leucura.fasta"
    
    echo "Extracted sequence saved to: ${OUTPUT_DIR}/${RAG_NORMALIZED}_Tityra_leucura.fasta"
    
    # Show the sequence length
    SEQ_LENGTH=$(grep -v ">" "${OUTPUT_DIR}/${RAG_NORMALIZED}_Tityra_leucura.fasta" | tr -d '\n' | wc -c)
    echo "Sequence length: ${SEQ_LENGTH} bp"
    
done

echo ""
echo "==================================================================="
echo "RAG gene extraction completed"
echo "Completed: $(date)"
echo "==================================================================="
echo ""
echo "Extracted sequences:"
ls -lh "${OUTPUT_DIR}"/RAG*_Tityra_leucura.fasta 2>/dev/null || echo "No sequences extracted"
