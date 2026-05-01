#!/bin/bash

###############################################################################
# MitoFinder Pipeline for Mitochondrial Genome Extraction from Assembly
# Description: Extract mitochondrial genomes from assembled genome
#              for Tityra using reference mitochondrial genomes
###############################################################################

set -euo pipefail

# Set working directory
WD=/media/inter/mkapun/projects/Tityra

# Output directory
RESULTS="${WD}/results/mitofinder_assembly"
mkdir -p "${RESULTS}"

# Create log file
LOG="${WD}/logs/mitofinder_assembly.log"
mkdir -p "${WD}/logs"

# Redirect output to log file
#exec > >(tee -a "${LOG}") 2>&1

echo "==================================================================="
echo "MitoFinder Pipeline for Mitochondrial Genome Extraction from Assembly"
echo "Started: $(date)"
echo "==================================================================="

# Mitochondrial genome accessions for reference species
declare -A MITO_GENOMES=(
    ["Tyrannus_savana"]="NC_051025.1"
    ["Pitta_sordida"]="NC_051463.1"
    ["Piprites_chloris"]="MN356298.1"
    ["Oxyruncus_cristatus"]="NC_053052.1"
    ["Onychorhynchus_coronatus"]="NC_053085.1"
    ["Cephalopterus_ornatus"]="NC_051008.1"
    ["Grallaria_varia"]="NC_051009.1"
    ["Neopipo_cinnamomea"]="NC_053075.1"
    ["Tachuris_rubrigastra"]="MN356165.1"
    ["Pachyramphus_minor"]="NC_051035.1"
)

# Download reference mitochondrial genomes
echo "Downloading reference mitochondrial genomes..."
sh ${WD}/shell/download_refseq_genomes.sh

# MitoFinder singularity container
MITOFINDER="/opt/bioinformatics/containers/mitofinder_v1.4.1.sif"

# RefSeq GenBank directory
REFSEQ_GB_DIR="${WD}/data/refseq_genbank"

# Assembled genome
ASSEMBLY_FILE="${WD}/results/denovo/output/Tityra_ILL.fa.gz"

# Combine all reference mitogenomes into a single GenBank file
echo ""
echo "==================================================================="
echo "Combining reference mitochondrial genomes..."
echo "==================================================================="

COMBINED_REF="${REFSEQ_GB_DIR}/all_references_combined.gb"
rm -f "${COMBINED_REF}"

for species in "${!MITO_GENOMES[@]}"; do
    REF_GB="${REFSEQ_GB_DIR}/${species}.gb"
    if [ -f "${REF_GB}" ]; then
        echo "Adding ${species} to combined reference..."
        cat "${REF_GB}" >> "${COMBINED_REF}"
    else
        echo "Warning: ${species}.gb not found, skipping..."
    fi
done

echo "Combined reference created: ${COMBINED_REF}"

###############################################################################
# Step 1: Extract mitochondrial genome from Tityra assembled genome
###############################################################################

echo ""
echo "==================================================================="
echo "Step 1: Extracting mitochondrial genome from Tityra assembly"
echo "==================================================================="

# Create output directory
TITYRA_OUT="${RESULTS}/Tityra_leucura"
mkdir -p "${TITYRA_OUT}"

# Check if assembly exists
if [ ! -f "${ASSEMBLY_FILE}" ]; then
    echo "Error: Assembly file not found: ${ASSEMBLY_FILE}"
    exit 1
fi

# Check if combined reference exists
if [ ! -f "${COMBINED_REF}" ]; then
    echo "Error: Combined reference not found: ${COMBINED_REF}"
    exit 1
fi

echo "Using assembly: ${ASSEMBLY_FILE}"
echo "Using combined reference with ${#MITO_GENOMES[@]} mitochondrial genomes"

# Prepare uncompressed assembly for MitoFinder
ASSEMBLY_UNZIPPED="${TITYRA_OUT}/Tityra_ILL.fa"
if [ ! -f "${ASSEMBLY_UNZIPPED}" ]; then
    echo "Uncompressing assembly..."
    gunzip -c "${ASSEMBLY_FILE}" > "${ASSEMBLY_UNZIPPED}"
fi

# Run MitoFinder
echo "Running MitoFinder on assembly..."
cd "${TITYRA_OUT}"

singularity exec --bind "${WD}:/work" "${MITOFINDER}" mitofinder \
    -a "/work/results/mitofinder_assembly/Tityra_leucura/Tityra_ILL.fa" \
    -j "Tityra_leucura" \
    -o 5 \
    -r "/work/data/refseq_genbank/all_references_combined.gb" \
    -p 10 \
    --override

echo "Completed Tityra mitochondrial genome extraction from assembly"

echo ""
echo "==================================================================="
echo "MitoFinder Pipeline Completed"
echo "Finished: $(date)"
echo "==================================================================="
