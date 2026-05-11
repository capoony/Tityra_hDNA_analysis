#!/bin/bash
###############################################################################
# Download complete mitochondrial genomes from GenBank
# Author: Generated for Tityra project
# Description: Search and download complete mitogenomes for species in tree
#              Uses accession numbers from Tyrannini (download_refseq_genomes.sh)
###############################################################################

WD="/media/inter/mkapun/projects/Tityra"

# Create output directory
mkdir -p "${WD}/data/refseq_genbank"

# Activate Python environment with Biopython
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs

echo "==================================================================="
echo "Downloading Mitochondrial Genomes from GenBank"
echo "Started: $(date)"
echo "==================================================================="

# Mitochondrial genome accession numbers from RefSeq/GenBank
# These are the Tyrannini taxa from download_refseq_genomes.sh
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
    output_gb="${WD}/data/refseq_genbank/${species}.gb"
    output_fasta="${WD}/data/refseq_genbank/${species}.fasta"
    
    echo ""
    echo "-------------------------------------------------------------------"
    echo "Downloading ${species} mitochondrial genome (${accession})..."
    
    # Download GenBank format
    if [ -f "${output_gb}" ]; then
        echo "GenBank file already exists: ${output_gb}"
    else
        echo "Downloading GenBank format..."
        curl -s "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log\$=seqview&db=nuccore&report=gbwithparts&id=${accession}&withparts=on" \
            -o "${output_gb}"
        
        if [ ! -s "${output_gb}" ]; then
            echo "Error: Failed to download GenBank for ${species}"
            rm -f "${output_gb}"
        else
            echo "Downloaded: ${output_gb}"
        fi
    fi
    
    # Download FASTA format
    if [ -f "${output_fasta}" ]; then
        echo "FASTA file already exists: ${output_fasta}"
    else
        echo "Downloading FASTA format..."
        curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${accession}&rettype=fasta&retmode=text" \
            -o "${output_fasta}"
        
        if [ ! -s "${output_fasta}" ]; then
            echo "Error: Failed to download FASTA for ${species}"
            rm -f "${output_fasta}"
        else
            echo "Downloaded: ${output_fasta}"
        fi
    fi
done

echo ""
echo "==================================================================="
echo "Download complete!"
echo "Completed: $(date)"
echo "==================================================================="

conda deactivate

## Note: Tityra genome will be added after reorientation in main.sh
## The reorient_circularize_tityra.sh script will create the corrected files
