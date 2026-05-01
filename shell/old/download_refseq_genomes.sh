#!/bin/bash

###############################################################################
# Download Mitochondrial Genomes in GenBank Format
# For species in the Tityra phylogeny analysis
###############################################################################

set -euo pipefail

WD=/media/inter/mkapun/projects/Tityra
DOWNLOAD_DIR="${WD}/data/refseq_genbank"
mkdir -p "${DOWNLOAD_DIR}"

echo "==================================================================="
echo "Downloading Mitochondrial Genomes in GenBank Format"
echo "Started: $(date)"
echo "==================================================================="

# Mitochondrial genome accession numbers from RefSeq
# Format: "Species_name:NC_accession"
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

# Download mitochondrial genomes in GenBank format
for species in "${!MITO_GENOMES[@]}"; do
    accession="${MITO_GENOMES[$species]}"
    output_file="${DOWNLOAD_DIR}/${species}.gb"
    
    echo ""
    echo "-------------------------------------------------------------------"
    echo "Downloading ${species} mitochondrial genome (${accession})..."
    
    if [ -f "${output_file}" ]; then
        echo "Already downloaded: ${output_file}"
        continue
    fi
    
    # Download using efetch from NCBI
    echo "Downloading mitochondrial genome in GenBank format..."
    curl -s "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log\$=seqview&db=nuccore&report=gbwithparts&id=${accession}&withparts=on" \
        -o "${output_file}"
    
    # Check if download was successful
    if [ ! -s "${output_file}" ]; then
        echo "Error: Failed to download ${species}"
        rm -f "${output_file}"
        continue
    fi
    
    echo "Downloaded: ${output_file}"
done

echo ""
echo "==================================================================="
echo "Download complete!"
echo "Completed: $(date)"
echo "==================================================================="
