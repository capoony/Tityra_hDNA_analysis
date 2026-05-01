#!/bin/bash

###############################################################################
# Graph-Based BUSCO Extraction Pipeline
# Reduces reference bias by using variation graphs from multiple references
# Author: Martin Kapun
# Date: April 2026
###############################################################################

set -euo pipefail

# Input parameters
WD=$1
THREADS=${2:-50}

# Export WD so Python scripts can access it via os.environ
export WD

# Always use vg v1.59.0 from local directory
VG="${WD}/results/graph_alignment/vg"

if [ ! -f "${VG}" ]; then
    echo "Downloading vg v1.59.0..."
    mkdir -p ${WD}/results/graph_alignment
    cd ${WD}/results/graph_alignment
    wget -q https://github.com/vgteam/vg/releases/download/v1.59.0/vg
    chmod +x vg
    cd ${WD}
fi

echo "Using vg: ${VG}"
${VG} version

echo "==================================================================="
echo "Graph-Based BUSCO Extraction Pipeline"
echo "Working directory: ${WD}"
echo "Threads: ${THREADS}"
echo "==================================================================="

###############################################################################
# 1. Setup and download references
###############################################################################

echo "Step 1: Setting up references for graph construction..."

mkdir -p ${WD}/results/graph_alignment/references
cd ${WD}/results/graph_alignment

# Use multiple close references to build graph - reuse from phylogeny if available
# Primary reference: Pachyramphus minor (closest relative)
if [ -f "${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz" ]; then
    echo "Reusing Pachyramphus_minor from data/ref..."
    ln -sf ${WD}/data/ref/GCA_013397135.1_ASM1339713v1_genomic.fna.gz \
        ${WD}/results/graph_alignment/references/Pachyramphus_minor.fna.gz
elif [ -f "${WD}/results/phylogeny/data/Pachyramphus_minor.fna.gz" ]; then
    echo "Reusing Pachyramphus_minor from phylogeny/data..."
    ln -sf ${WD}/results/phylogeny/data/Pachyramphus_minor.fna.gz \
        ${WD}/results/graph_alignment/references/Pachyramphus_minor.fna.gz
fi

# Additional references from same family (Cotingidae/Tyrannidae)
declare -A REFS=(
    ["Oxyruncus_cristatus"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/399/855/GCA_013399855.1_ASM1339985v1/GCA_013399855.1_ASM1339985v1_genomic.fna.gz"
    ["Tyrannus_savana"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/399/735/GCA_013399735.1_ASM1339973v1/GCA_013399735.1_ASM1339973v1_genomic.fna.gz"
    ["Piprites_chloris"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/399/295/GCA_013399295.1_ASM1339929v1/GCA_013399295.1_ASM1339929v1_genomic.fna.gz"
)

# Download or reuse additional references
for name in "${!REFS[@]}"; do
    url="${REFS[$name]}"
    
    # Check if genome exists in phylogeny data directory
    if [ -f "${WD}/results/phylogeny/data/${name}.fna.gz" ]; then
        echo "Reusing ${name} from phylogeny/data..."
        ln -sf ${WD}/results/phylogeny/data/${name}.fna.gz \
            ${WD}/results/graph_alignment/references/${name}.fna.gz
    elif [ ! -f "${WD}/results/graph_alignment/references/${name}.fna.gz" ]; then
        echo "Downloading ${name}..."
        wget -q -O "${WD}/results/graph_alignment/references/${name}.fna.gz" "${url}"
    else
        echo "${name} already present in graph_alignment/references"
    fi
done

###############################################################################
# 2. Run BUSCO on all references to identify conserved regions
###############################################################################

echo "Step 2: Checking for existing BUSCO results or running BUSCO on reference genomes..."

source /opt/anaconda3/etc/profile.d/conda.sh 2>/dev/null || source /opt/mamba/etc/profile.d/conda.sh 2>/dev/null
conda activate busco_5.4.3

mkdir -p ${WD}/results/graph_alignment/busco_refs

# Map reference file names to potential genome names in phylogeny pipeline
declare -A ref_to_genome=(
    ["Pachyramphus_minor"]="Pachyramphus_minor"
    ["Oxyruncus_cristatus"]="Oxyruncus_cristatus"
    ["Tyrannus_savana"]="Tyrannus_savana"
    ["Piprites_chloris"]="Piprites_chloris"
)

for ref_gz in ${WD}/results/graph_alignment/references/*.fna.gz; do
    ref_name=$(basename ${ref_gz} .fna.gz)
    
    # Try to find corresponding genome name
    genome_name=""
    for key in "${!ref_to_genome[@]}"; do
        if [[ "${ref_name}" == *"${key}"* ]]; then
            genome_name="${ref_to_genome[$key]}"
            break
        fi
    done
    
    # Check if BUSCO results already exist from phylogeny pipeline
    if [ -n "${genome_name}" ] && [ -d "${WD}/results/phylogeny/BUSCO/${genome_name}" ]; then
        echo "Reusing existing BUSCO results for ${genome_name}..."
        ln -sf ${WD}/results/phylogeny/BUSCO/${genome_name} \
               ${WD}/results/graph_alignment/busco_refs/${ref_name}
        continue
    fi
    
    # Otherwise, check if we already ran BUSCO in this directory
    if [ -d "${WD}/results/graph_alignment/busco_refs/${ref_name}" ]; then
        echo "BUSCO already run for ${ref_name}, skipping..."
        continue
    fi
    
    # Need to run BUSCO
    echo "Running BUSCO on ${ref_name}..."
    
    # Decompress reference
    gunzip -c ${ref_gz} > ${WD}/results/graph_alignment/references/${ref_name}.fna
    
    cd ${WD}/results/graph_alignment/busco_refs
    
    busco -i ${WD}/results/graph_alignment/references/${ref_name}.fna \
        -o ${ref_name} \
        -m genome \
        -c ${THREADS} \
        -l aves_odb10 \
        -f \
        --quiet
done

conda deactivate

###############################################################################
# 3. Extract BUSCO gene regions from references and build graphs
###############################################################################

echo "Step 3: Extracting BUSCO regions and building variation graphs..."

mkdir -p ${WD}/results/graph_alignment/busco_regions
mkdir -p ${WD}/results/graph_alignment/graphs

# Python script to extract BUSCO coordinates and sequences
python3 << 'PYTHON_EXTRACT'
import os
import sys
from Bio import SeqIO

wd = os.environ['WD']
busco_dir = f"{wd}/results/graph_alignment/busco_refs"
output_dir = f"{wd}/results/graph_alignment/busco_regions"

# Find all complete BUSCO genes across references
busco_genes = set()
busco_refs = {}

for busco_result in os.listdir(busco_dir):
    if not busco_result.endswith('_busco'):
        continue
    
    ref_name = busco_result.replace('_busco', '')
    seq_dir = f"{busco_dir}/{busco_result}/run_aves_odb10/busco_sequences/single_copy_busco_sequences"
    
    if not os.path.exists(seq_dir):
        print(f"Warning: No sequences found for {ref_name}", file=sys.stderr)
        continue
    
    busco_refs[ref_name] = seq_dir
    
    # Find all .fna (nucleotide coding sequence) files
    for fname in os.listdir(seq_dir):
        if fname.endswith('.fna'):
            gene_id = fname.replace('.fna', '')
            busco_genes.add(gene_id)

print(f"Found {len(busco_genes)} BUSCO genes across {len(busco_refs)} references", file=sys.stderr)

# For each BUSCO gene, collect sequences from all references
extracted = 0
for gene in sorted(busco_genes):
    gene_seqs = []
    
    for ref_name, seq_dir in busco_refs.items():
        fna_file = f"{seq_dir}/{gene}.fna"
        
        if not os.path.exists(fna_file):
            continue
        
        # Read the BUSCO-extracted coding sequence
        for record in SeqIO.parse(fna_file, "fasta"):
            gene_seqs.append(f">{ref_name}\n{str(record.seq)}\n")
            break  # Only take first sequence (should be only one)
    
    # Write multi-FASTA only if we have sequences from at least 2 references
    if len(gene_seqs) >= 2:
        with open(f"{output_dir}/{gene}.fasta", 'w') as out:
            out.writelines(gene_seqs)
        extracted += 1

print(f"Extracted {extracted} BUSCO genes with >=2 reference sequences", file=sys.stderr)
PYTHON_EXTRACT

echo "Extracted $(ls ${WD}/results/graph_alignment/busco_regions/*.fasta 2>/dev/null | wc -l) BUSCO gene regions"

###############################################################################
# 4. Align BUSCO sequences with MAFFT
###############################################################################

echo "Step 4: Aligning BUSCO sequences with MAFFT..."

mkdir -p ${WD}/results/graph_alignment/aligned

# Initialize and activate MAFFT conda environment
source /opt/anaconda3/etc/profile.d/conda.sh || source /opt/mamba/etc/profile.d/conda.sh
conda deactivate
conda activate mafft-7.487

for busco_fa in ${WD}/results/graph_alignment/busco_regions/*.fasta; do
    gene=$(basename ${busco_fa} .fasta)
    aligned_fa=${WD}/results/graph_alignment/aligned/${gene}_aligned.fasta
    
    if [ -f "${aligned_fa}" ]; then
        echo "${gene} already aligned, skipping..."
        continue
    fi
    
    echo "Aligning ${gene}..."
    mafft --auto --thread 4 ${busco_fa} > ${aligned_fa} 2>/dev/null
done

conda deactivate

###############################################################################
# 5. Build variation graphs using vg
###############################################################################

echo "Step 5: Building variation graphs with vg..."

# VG already set at top of script, activate conda environment
source /opt/anaconda3/etc/profile.d/conda.sh || source /opt/mamba/etc/profile.d/conda.sh
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs

# Build graphs from aligned sequences
for aligned_fa in ${WD}/results/graph_alignment/aligned/*_aligned.fasta; do
    gene=$(basename ${aligned_fa} _aligned.fasta)
    
    if [ -f "${WD}/results/graph_alignment/graphs/${gene}.xg" ]; then
        echo "Graph for ${gene} already exists, skipping..."
        continue
    fi
    
    echo "Building graph for ${gene}..."
    
    # Construct variation graph from aligned MSA
    ${VG} construct \
        -M ${aligned_fa} \
        -m 32 \
        > ${WD}/results/graph_alignment/graphs/${gene}.vg 2>/dev/null
    
    # Index for mapping
    ${VG} index \
        -x ${WD}/results/graph_alignment/graphs/${gene}.xg \
        ${WD}/results/graph_alignment/graphs/${gene}.vg 2>/dev/null
    
    ${VG} prune \
        ${WD}/results/graph_alignment/graphs/${gene}.vg \
        > ${WD}/results/graph_alignment/graphs/${gene}.pruned.vg 2>/dev/null
    
    ${VG} index \
        -g ${WD}/results/graph_alignment/graphs/${gene}.gcsa \
        -k 16 \
        ${WD}/results/graph_alignment/graphs/${gene}.pruned.vg 2>/dev/null
    
    # Clean up intermediate files
    rm -f ${WD}/results/graph_alignment/graphs/${gene}.vg
    rm -f ${WD}/results/graph_alignment/graphs/${gene}.pruned.vg
done

###############################################################################
# 6. Map Tityra reads to graphs and extract consensus
###############################################################################

echo "Step 6: Mapping Tityra reads to variation graphs..."

mkdir -p ${WD}/results/graph_alignment/mapped
mkdir -p ${WD}/results/graph_alignment/consensus

READS_PE1="${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz"
READS_PE2="${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz"
READS_MERGED="${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz"

for graph_xg in ${WD}/results/graph_alignment/graphs/*.xg; do
    gene=$(basename ${graph_xg} .xg)
    
    if [ -f "${WD}/results/graph_alignment/consensus/${gene}_Tityra.fasta" ]; then
        echo "Consensus for ${gene} already exists, skipping..."
        continue
    fi
    
    echo "Mapping reads to ${gene}..."
    
    # Map paired-end reads
    ${VG} map \
        -x ${graph_xg} \
        -g ${WD}/results/graph_alignment/graphs/${gene}.gcsa \
        -f ${READS_PE1} -f ${READS_PE2} \
        -t ${THREADS} \
        > ${WD}/results/graph_alignment/mapped/${gene}.gam 2>/dev/null
    
    # Map merged reads
    ${VG} map \
        -x ${graph_xg} \
        -g ${WD}/results/graph_alignment/graphs/${gene}.gcsa \
        -f ${READS_MERGED} \
        -t ${THREADS} \
        >> ${WD}/results/graph_alignment/mapped/${gene}.gam 2>/dev/null
    
    # Pack alignments for variant calling
    ${VG} pack \
        -x ${graph_xg} \
        -g ${WD}/results/graph_alignment/mapped/${gene}.gam \
        -o ${WD}/results/graph_alignment/mapped/${gene}.pack \
        2>/dev/null
    
    # Check mapping stats (for debugging)
    echo "  Read count: $(${VG} stats -a ${WD}/results/graph_alignment/mapped/${gene}.gam | grep 'Total alignments' || echo 'N/A')"
    
    # Call variants with relaxed filters for low-coverage ancient DNA
    # -d 1: minimum depth of 1 (default is 4)
    # -s 1: minimum support for variant of 1 read (default is 1)
    ${VG} call \
        ${graph_xg} \
        -k ${WD}/results/graph_alignment/mapped/${gene}.pack \
        -d 1 \
        -s 1 \
        > ${WD}/results/graph_alignment/mapped/${gene}_raw.vcf \
        2>/dev/null
    
    # Extract reference sequence from BUSCO file for consensus generation
    # Use Pachyramphus_minor as reference (primary reference in graph)
    samtools faidx ${WD}/results/graph_alignment/busco_regions/${gene}.fasta \
        ref1_Pachyramphus_minor \
        > ${WD}/results/graph_alignment/mapped/${gene}_ref.fasta 2>/dev/null
    
    # Filter variants: only keep those with depth > 5
    bcftools view \
        -i 'FORMAT/DP>5' \
        ${WD}/results/graph_alignment/mapped/${gene}_raw.vcf \
        -o ${WD}/results/graph_alignment/mapped/${gene}_filtered.vcf \
        2>/dev/null
    
    # Compress and index filtered VCF for bcftools
    bgzip -c ${WD}/results/graph_alignment/mapped/${gene}_filtered.vcf > ${WD}/results/graph_alignment/mapped/${gene}.vcf.gz 2>/dev/null
    tabix -p vcf ${WD}/results/graph_alignment/mapped/${gene}.vcf.gz 2>/dev/null
    
    # Index reference FASTA
    samtools faidx ${WD}/results/graph_alignment/mapped/${gene}_ref.fasta 2>/dev/null
    
    # Apply variants to get Tityra consensus
    bcftools consensus \
        -f ${WD}/results/graph_alignment/mapped/${gene}_ref.fasta \
        ${WD}/results/graph_alignment/mapped/${gene}.vcf.gz \
        -o ${WD}/results/graph_alignment/consensus/${gene}_consensus_tmp.fasta \
        2>/dev/null
    
    # Rename sequence header
    sed "s/>ref1_Pachyramphus_minor/>${gene}|Tityra_leucura_graph/" \
        ${WD}/results/graph_alignment/consensus/${gene}_consensus_tmp.fasta \
        > ${WD}/results/graph_alignment/consensus/${gene}_Tityra.fasta
    
    # Report variant count
    variant_count=$(bcftools view -H ${WD}/results/graph_alignment/mapped/${gene}.vcf.gz | wc -l)
    echo "  Variants (DP>5): ${variant_count}"
    
    # Clean up intermediate files
    rm -f ${WD}/results/graph_alignment/mapped/${gene}.gam
    rm -f ${WD}/results/graph_alignment/mapped/${gene}_ref.fasta
    rm -f ${WD}/results/graph_alignment/mapped/${gene}_ref.fasta.fai
    rm -f ${WD}/results/graph_alignment/mapped/${gene}_raw.vcf
    rm -f ${WD}/results/graph_alignment/mapped/${gene}_filtered.vcf
    rm -f ${WD}/results/graph_alignment/consensus/${gene}_consensus_tmp.fasta
done

###############################################################################
# 7. Combine all BUSCO genes for phylogeny
###############################################################################

echo "Step 7: Creating combined BUSCO gene set..."

mkdir -p ${WD}/results/phylogeny_graph/busco_aa

# Convert nucleotide sequences to protein using BioPython
python3 << 'PYTHON_TRANSLATE'
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq

wd = os.environ['WD']
consensus_dir = f"{wd}/results/graph_alignment/consensus"
output_dir = f"{wd}/results/phylogeny_graph/busco_aa"

translated = 0
for consensus_file in os.listdir(consensus_dir):
    if not consensus_file.endswith('_Tityra.fasta'):
        continue
    
    gene = consensus_file.replace('_Tityra.fasta', '')
    
    # Read consensus nucleotide sequence
    consensus_path = f"{consensus_dir}/{consensus_file}"
    
    try:
        for record in SeqIO.parse(consensus_path, "fasta"):
            # Translate in frame 1 (BUSCO genes should be proper CDS)
            protein_seq = record.seq.translate()
            
            # Write protein sequence
            with open(f"{output_dir}/{gene}.faa", 'w') as out:
                out.write(f">{gene}|Tityra_leucura_graph\n")
                out.write(f"{str(protein_seq)}\n")
            
            translated += 1
            break  # Only translate first sequence
    except Exception as e:
        print(f"Warning: Could not translate {gene}: {e}", file=sys.stderr)
        continue

print(f"Translated {translated} nucleotide consensus sequences to protein", file=sys.stderr)
PYTHON_TRANSLATE

echo "==================================================================="
echo "Graph-based BUSCO extraction complete!"
echo "Generated $(ls ${WD}/results/graph_alignment/consensus/*_Tityra.fasta 2>/dev/null | wc -l) consensus sequences"
echo "Protein sequences in: ${WD}/results/phylogeny_graph/busco_aa/"
echo "==================================================================="

conda deactivate
