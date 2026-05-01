#!/bin/bash

# Script to combine mitochondrial genes by individual specimen using voucher IDs
# This script creates concatenated gene sequences per individual

set -e

# Define paths
WD="/media/inter/mkapun/projects/Tityra"
DATA_DIR="${WD}/data/Mito"
COMBINED_DIR="${WD}/results/combined_mito_genes"
OUTPUT_DIR="${WD}/results/individual_concatenated"
MITOFINDER_TLEUC="${WD}/results/mitofinder/Tityra_leucura/Tityra_leucura/Tityra_leucura_MitoFinder_megahit_mitfi_Final_Results/Tityra_leucura_mtDNA_contig_1_genes_NT.fasta"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "=== Combining genes by individual specimen ==="
echo "Starting at $(date)"

# Define gene list (including synonyms)
GENES=("ATP6" "ATP8" "COI" "CO2" "CO3" "CYTB" "ND1" "ND2" "ND3" "ND4" "ND4L" "ND5")

# Extract voucher IDs from TSV files
echo "Step 1: Extracting voucher IDs and specimen information..."

# Create temporary file to store specimen information
SPECIMEN_INFO="${OUTPUT_DIR}/specimen_info.tmp"
> "${SPECIMEN_INFO}"

# Process all TSV files to extract voucher IDs and link them to genes
for TSV_FILE in ${DATA_DIR}/*.tsv; do
    # Get gene name and species from filename (format: GENE_Species.tsv)
    BASENAME=$(basename "${TSV_FILE}" .tsv)
    GENE=$(echo "${BASENAME}" | cut -d'_' -f1)
    SPECIES_FILE=$(echo "${BASENAME}" | cut -d'_' -f2-)
    
    # Extract voucher IDs from GenBankTitle (field 4)
    # Look for patterns like "voucher_XXX" or specimen identifiers
    awk -F'\t' -v gene="${GENE}" -v species_file="${SPECIES_FILE}" 'NR>1 {
        # Get accession and title
        accession = $2
        title = $4
        species = $3
        
        # Extract voucher ID from title
        # Look for "voucher" followed by identifier
        if (match(title, /voucher[_: ]+([A-Za-z0-9:_-]+)/, arr)) {
            voucher = arr[1]
        } else if (match(title, /[A-Z]+[:]?[A-Z0-9]+[:]?[0-9]+/, arr)) {
            # Try to extract museum/collection codes
            voucher = arr[0]
        } else {
            # Use accession as fallback
            voucher = accession
        }
        
        # Remove common gene name suffixes from voucher ID
        # Remove patterns like "_cytochrome", "_NADH", "_COX", etc.
        gsub(/_[Cc]ytochrome.*$/, "", voucher)
        gsub(/_NADH.*$/, "", voucher)
        gsub(/_[Cc]ytochrome.*$/, "", voucher)
        gsub(/_[Cc][Oo][XxIi].*$/, "", voucher)
        gsub(/_[Aa][Tt][Pp].*$/, "", voucher)
        gsub(/_ND[0-9].*$/, "", voucher)
        gsub(/_gene.*$/, "", voucher)
        gsub(/_mitochondrial.*$/, "", voucher)
        gsub(/_partial.*$/, "", voucher)
        gsub(/_complete.*$/, "", voucher)
        
        # Clean up trailing underscores and spaces
        gsub(/_+$/, "", voucher)
        gsub(/[ \t]+/, "_", voucher)
        
        # Print: voucher_id, gene, accession, species, species_file
        print voucher "\t" gene "\t" accession "\t" species "\t" species_file
    }' "${TSV_FILE}" >> "${SPECIMEN_INFO}"
done

echo "  Found $(wc -l < ${SPECIMEN_INFO}) gene-specimen combinations"

# Get unique voucher IDs
VOUCHER_LIST="${OUTPUT_DIR}/voucher_list.txt"
cut -f1 "${SPECIMEN_INFO}" | sort -u > "${VOUCHER_LIST}"
echo "  Found $(wc -l < ${VOUCHER_LIST}) unique specimens"

# Create gene presence matrix
echo "Step 2: Building gene presence matrix..."
GENE_MATRIX="${OUTPUT_DIR}/gene_matrix.txt"
echo -ne "Voucher_ID" > "${GENE_MATRIX}"
for GENE in "${GENES[@]}"; do
    echo -ne "\t${GENE}" >> "${GENE_MATRIX}"
done
echo "" >> "${GENE_MATRIX}"

# For each voucher, check which genes are available
while IFS= read -r VOUCHER; do
    echo -ne "${VOUCHER}" >> "${GENE_MATRIX}"
    for GENE in "${GENES[@]}"; do
        # Check if this voucher has this gene
        COUNT=$(grep -c "^${VOUCHER}[[:space:]]${GENE}[[:space:]]" "${SPECIMEN_INFO}" || true)
        echo -ne "\t${COUNT}" >> "${GENE_MATRIX}"
    done
    echo "" >> "${GENE_MATRIX}"
done < "${VOUCHER_LIST}"

echo "  Gene presence matrix created: ${GENE_MATRIX}"

# Analyze missing data
echo "Step 3: Analyzing missing data patterns..."
STATS_FILE="${OUTPUT_DIR}/missing_data_stats.txt"
echo "=== Missing Data Statistics ===" > "${STATS_FILE}"
echo "" >> "${STATS_FILE}"

# Calculate statistics
awk -F'\t' 'NR>1 {
    total_genes = NF - 1
    present = 0
    for (i=2; i<=NF; i++) {
        if ($i > 0) present++
    }
    missing = total_genes - present
    print present "\t" missing "\t" $1
}' "${GENE_MATRIX}" | sort -rn -k1,1 > "${OUTPUT_DIR}/specimen_completeness.txt"

echo "Specimens sorted by gene completeness:" >> "${STATS_FILE}"
echo "(Present_Genes  Missing_Genes  Voucher_ID)" >> "${STATS_FILE}"
head -20 "${OUTPUT_DIR}/specimen_completeness.txt" >> "${STATS_FILE}"

# Determine optimal threshold
echo "" >> "${STATS_FILE}"
echo "=== Threshold Analysis ===" >> "${STATS_FILE}"
for MIN_GENES in 1 2 3 4 5 6 8 10 12; do
    COUNT=$(awk -v min="${MIN_GENES}" '$1 >= min' "${OUTPUT_DIR}/specimen_completeness.txt" | wc -l)
    echo "Specimens with >= ${MIN_GENES} genes: ${COUNT}" >> "${STATS_FILE}"
done

cat "${STATS_FILE}"

# Select optimal threshold (balance between completeness and sample size)
# Default: minimum 2 genes per specimen
MIN_GENES_THRESHOLD=2

echo ""
echo "Step 4: Extracting sequences for specimens with >= ${MIN_GENES_THRESHOLD} genes..."

# Get vouchers that meet the threshold
awk -v min="${MIN_GENES_THRESHOLD}" '$1 >= min {print $3}' "${OUTPUT_DIR}/specimen_completeness.txt" > "${OUTPUT_DIR}/candidates.txt"

# For each voucher, add species information
> "${OUTPUT_DIR}/candidates_with_species.txt"
while IFS= read -r VOUCHER; do
    SPECIES=$(grep "^${VOUCHER}[[:space:]]" "${SPECIMEN_INFO}" | head -1 | cut -f4 | sed 's/ /_/g')
    if [ -n "${SPECIES}" ]; then
        echo "${SPECIES}	${VOUCHER}" >> "${OUTPUT_DIR}/candidates_with_species.txt"
    fi
done < "${OUTPUT_DIR}/candidates.txt"

# Select maximum 2 random specimens per species
> "${OUTPUT_DIR}/selected_vouchers.txt"
for SPECIES in $(cut -f1 "${OUTPUT_DIR}/candidates_with_species.txt" | sort -u); do
    grep "^${SPECIES}[[:space:]]" "${OUTPUT_DIR}/candidates_with_species.txt" | 
        shuf | head -2 | cut -f2 >> "${OUTPUT_DIR}/selected_vouchers.txt"
done

# Add special voucher for MitoFinder Tityra leucura
echo "Tityra_leucura_MitoFinder" >> "${OUTPUT_DIR}/selected_vouchers.txt"

SELECTED_COUNT=$(wc -l < "${OUTPUT_DIR}/selected_vouchers.txt")
echo "  Selected ${SELECTED_COUNT} specimens (max 2 per species, >= ${MIN_GENES_THRESHOLD} genes)"

# For each selected voucher, extract sequences directly from TSV files
CONCAT_FILE="${OUTPUT_DIR}/concatenated_sequences.fasta"
> "${CONCAT_FILE}"

echo "  Extracting sequences..."
COUNT=0

while IFS= read -r VOUCHER; do
    COUNT=$((COUNT + 1))
    if [ $COUNT -le 3 ] || [ $((COUNT % 10)) -eq 0 ]; then
        echo "    Processing $COUNT: ${VOUCHER}"
    fi
    
    # Special handling for Tityra leucura from MitoFinder
    if [ "${VOUCHER}" == "Tityra_leucura_MitoFinder" ]; then
        SPECIES="Tityra_leucura"
        CONCAT_SEQ=""
        GENE_ORDER=""
        
        # Map gene names to MitoFinder nomenclature
        for GENE in "${GENES[@]}"; do
            MITOFINDER_GENE="${GENE}"
            # Convert gene names to MitoFinder format
            case "${GENE}" in
                "COI") MITOFINDER_GENE="COX1" ;;
                "CO2") MITOFINDER_GENE="COX2" ;;
                "CO3") MITOFINDER_GENE="COX3" ;;
            esac
            
            # Extract sequence from MitoFinder FASTA
            if [ -f "${MITOFINDER_TLEUC}" ]; then
                SEQ=$(awk -v gene="${MITOFINDER_GENE}" '
                    /^>/ {
                        if (index($0, "@" gene) > 0) {
                            found=1
                            next
                        } else {
                            found=0
                        }
                    }
                    found {
                        gsub(/[[:space:]]/, "", $0)
                        seq = seq $0
                    }
                    END {
                        print seq
                    }
                ' "${MITOFINDER_TLEUC}")
                
                if [ -n "${SEQ}" ]; then
                    CONCAT_SEQ="${CONCAT_SEQ}${SEQ}"
                    GENE_ORDER="${GENE_ORDER}${GENE}+"
                fi
            fi
        done
        
        # Write concatenated sequence
        if [ -n "${CONCAT_SEQ}" ]; then
            GENE_ORDER=${GENE_ORDER%+}
            echo ">${SPECIES}_${VOUCHER} genes:${GENE_ORDER}" >> "${CONCAT_FILE}"
            echo "${CONCAT_SEQ}" >> "${CONCAT_FILE}"
        fi
        continue
    fi
    
    # Get species name for this voucher
    SPECIES=$(grep "^${VOUCHER}[[:space:]]" "${SPECIMEN_INFO}" | head -1 | cut -f4 | sed 's/ /_/g')
    
    # Extract sequences for each gene
    CONCAT_SEQ=""
    GENE_ORDER=""
    
    for GENE in "${GENES[@]}"; do
        # Get accession and species info for this voucher and gene
        INFO=$(grep "^${VOUCHER}[[:space:]]${GENE}[[:space:]]" "${SPECIMEN_INFO}" | head -1)
        ACCESSION=$(echo "${INFO}" | cut -f3)
        SPECIES_FILE=$(echo "${INFO}" | cut -f5)
        
        if [ -n "${ACCESSION}" ] && [ -n "${SPECIES_FILE}" ]; then
            # Find the TSV file for this gene and species filename
            TSV_FILE="${DATA_DIR}/${GENE}_${SPECIES_FILE}.tsv"
            
            if [ -f "${TSV_FILE}" ]; then
                # Extract sequence from TSV file (column 13 contains the sequence)
                SEQ=$(awk -F'\t' -v acc="${ACCESSION}" '
                    $2 == acc {
                        # Remove all whitespace and newlines from sequence
                        gsub(/[[:space:]]/, "", $13)
                        print $13
                        exit
                    }
                ' "${TSV_FILE}")
                
                if [ -n "${SEQ}" ]; then
                    CONCAT_SEQ="${CONCAT_SEQ}${SEQ}"
                    GENE_ORDER="${GENE_ORDER}${GENE}+"
                    if [ $COUNT -le 2 ]; then
                        echo "      Found ${GENE} sequence (length: ${#SEQ})"
                    fi
                else
                    if [ $COUNT -le 2 ]; then
                        echo "      No sequence found for ${GENE} (acc: ${ACCESSION})"
                    fi
                fi
            else
                if [ $COUNT -le 2 ]; then
                    echo "      TSV file not found: ${TSV_FILE}"
                fi
            fi
        else
            if [ $COUNT -le 2 ] && [ -n "${GENE}" ]; then
                echo "      No info for ${GENE} (acc: '${ACCESSION}', file: '${SPECIES_FILE}')"
            fi
        fi
    done
    
    # Write concatenated sequence if we have any data
    if [ -n "${CONCAT_SEQ}" ]; then
        GENE_ORDER=${GENE_ORDER%+}  # Remove trailing +
        echo ">${SPECIES}_${VOUCHER} genes:${GENE_ORDER}" >> "${CONCAT_FILE}"
        echo "${CONCAT_SEQ}" >> "${CONCAT_FILE}"
    fi
    
done < "${OUTPUT_DIR}/selected_vouchers.txt"

echo ""
echo "Step 5: Creating gene partition file for phylogenetic analysis..."
PARTITION_FILE="${OUTPUT_DIR}/gene_partitions.txt"
> "${PARTITION_FILE}"

# Analyze the concatenated sequences to determine gene boundaries
echo "# Gene partitions for concatenated alignment" > "${PARTITION_FILE}"
echo "# Format: gene_name = start-end" >> "${PARTITION_FILE}"
echo "#" >> "${PARTITION_FILE}"

# This will be generated after alignment
echo "# Note: Gene boundaries need to be determined after alignment" >> "${PARTITION_FILE}"
echo "# Check individual sequences for gene order" >> "${PARTITION_FILE}"

echo ""
echo "=== Summary ==="
echo "Total specimens processed: $(wc -l < ${VOUCHER_LIST})"
echo "Specimens selected (>= ${MIN_GENES_THRESHOLD} genes): ${SELECTED_COUNT}"
echo "Concatenated sequences written to: ${CONCAT_FILE}"
echo ""
echo "Output files:"
echo "  - Gene presence matrix: ${GENE_MATRIX}"
echo "  - Missing data statistics: ${STATS_FILE}"
echo "  - Concatenated sequences: ${CONCAT_FILE}"
echo "  - Partition file: ${PARTITION_FILE}"
echo ""
echo "Completed at $(date)"

# Cleanup temporary files (keep specimen_info for debugging)
rm -f "${OUTPUT_DIR}/specimen_completeness.txt" "${OUTPUT_DIR}/selected_vouchers.txt"

echo ""
echo "Next steps:"
echo "1. Align concatenated sequences with MAFFT"
echo "2. Use gene partitions for partitioned phylogenetic analysis"
echo "3. Consider different minimum gene thresholds if needed"
