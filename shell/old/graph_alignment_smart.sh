#!/bin/bash

###############################################################################
# Smart Two-Stage Graph-Based Variant Calling
# Stage 1: Regular BWA mapping to identify 500 best genes (by coverage)
# Stage 2: Concatenated graph mapping only for those 500 genes
###############################################################################

set -euo pipefail

WD=${1:-/media/inter/mkapun/projects/Tityra}
THREADS=${2:-50}
TOP_N=500  # Number of best genes to use for graph mapping
MIN_DEPTH=5

echo "==================================================================="
echo "Smart Two-Stage Graph-Based Variant Calling"
echo "Working directory: ${WD}"
echo "Threads: ${THREADS}"
echo "Stage 1: BWA mapping to all genes → identify top ${TOP_N} genes"
echo "Stage 2: Graph mapping to top ${TOP_N} genes only"
echo "==================================================================="

# Setup paths
VG="${WD}/results/graph_alignment/vg"
CONDA_ENV=/media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs

BWA=${CONDA_ENV}/bin/bwa
SAMTOOLS=${CONDA_ENV}/bin/samtools
BCFTOOLS=${CONDA_ENV}/bin/bcftools
BGZIP=${CONDA_ENV}/bin/bgzip
TABIX=${CONDA_ENV}/bin/tabix

PE1=${WD}/data/trimmed/Tityra_leucura_1_trimmed.fastq.gz
PE2=${WD}/data/trimmed/Tityra_leucura_2_trimmed.fastq.gz
MERGED=${WD}/data/trimmed/Tityra_leucura_merged.fastq.gz

STAGE1_DIR=${WD}/results/graph_alignment/stage1_bwa
STAGE2_DIR=${WD}/results/graph_alignment/stage2_graph

mkdir -p ${STAGE1_DIR}
mkdir -p ${STAGE2_DIR}
mkdir -p ${WD}/results/graph_alignment/consensus
mkdir -p ${WD}/results/graph_alignment/busco_regions
mkdir -p ${WD}/results/graph_alignment/aligned
mkdir -p ${WD}/results/graph_alignment/busco_refs

###############################################################################
# STAGE 0: Extract and align BUSCO genes (if needed)
###############################################################################

# Check if we have aligned genes, if not extract and align
ALIGNED_COUNT=$(ls ${WD}/results/graph_alignment/aligned/*.fasta 2>/dev/null | wc -l)

if [ ${ALIGNED_COUNT} -eq 0 ]; then
    echo ""
    echo "==================================================================="
    echo "STAGE 0: Extracting and Aligning BUSCO Genes"
    echo "==================================================================="
    
    # Reuse BUSCO results from phylogeny pipeline
    echo "Linking BUSCO results..."
    for genome in Pachyramphus_minor Oxyruncus_cristatus Tyrannus_savana; do
        if [ -d "${WD}/results/phylogeny/BUSCO/${genome}" ]; then
            ln -sf ${WD}/results/phylogeny/BUSCO/${genome} \
                   ${WD}/results/graph_alignment/busco_refs/${genome}
        fi
    done
    
    # Extract BUSCO gene sequences
    echo "Extracting BUSCO genes..."
    
    python3 << 'PYTHON_EXTRACT'
import os
import sys
from Bio import SeqIO
from collections import defaultdict

wd = os.environ['WD']
busco_dir = f"{wd}/results/graph_alignment/busco_refs"
output_dir = f"{wd}/results/graph_alignment/busco_regions"
phylo_data = f"{wd}/results/phylogeny/data"

# Find complete BUSCO genes
busco_genes = defaultdict(list)

for ref_name in os.listdir(busco_dir):
    ref_path = os.path.join(busco_dir, ref_name)
    full_table = os.path.join(ref_path, "run_aves_odb10/full_table.tsv")
    
    if not os.path.exists(full_table):
        continue
    
    with open(full_table) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue
            gene_id, status = fields[0], fields[1]
            if status == "Complete":
                busco_genes[gene_id].append(ref_name)

# Select genes in at least 2 references
selected_genes = {gene: refs for gene, refs in busco_genes.items() if len(refs) >= 2}
print(f"Found {len(selected_genes)} genes in at least 2 references")

# Extract sequences
for gene_id, ref_list in selected_genes.items():
    sequences = []
    
    for ref_name in ref_list:
        ref_path = os.path.join(busco_dir, ref_name)
        gene_fasta = os.path.join(ref_path, f"run_aves_odb10/busco_sequences/single_copy_busco_sequences/{gene_id}.fna")
        
        if os.path.exists(gene_fasta):
            for record in SeqIO.parse(gene_fasta, "fasta"):
                record.id = f"ref{len(sequences)+1}_{ref_name}"
                record.description = ""
                sequences.append(record)
    
    if len(sequences) >= 2:
        with open(f"{output_dir}/{gene_id}.fasta", "w") as f:
            SeqIO.write(sequences, f, "fasta")

print(f"Extracted {len(os.listdir(output_dir))} gene files")
PYTHON_EXTRACT
    
    # Align with MAFFT
    echo "Aligning genes with MAFFT..."
    
    source /opt/anaconda3/etc/profile.d/conda.sh 2>/dev/null || source /opt/mamba/etc/profile.d/conda.sh 2>/dev/null
    conda activate mafft-7.487
    
    for gene_file in ${WD}/results/graph_alignment/busco_regions/*.fasta; do
        gene=$(basename ${gene_file} .fasta)
        output_file=${WD}/results/graph_alignment/aligned/${gene}.fasta
        
        if [ ! -f "${output_file}" ]; then
            mafft --auto --quiet ${gene_file} > ${output_file} 2>/dev/null
        fi
    done
    
    conda deactivate
    
    ALIGNED_COUNT=$(ls ${WD}/results/graph_alignment/aligned/*.fasta 2>/dev/null | wc -l)
    echo "Aligned ${ALIGNED_COUNT} genes"
fi

if [ ${ALIGNED_COUNT} -eq 0 ]; then
    echo "ERROR: No aligned genes found"
    exit 1
fi

echo "Found ${ALIGNED_COUNT} aligned genes"

###############################################################################
# STAGE 1: BWA Mapping to identify best genes
###############################################################################

echo ""
echo "==================================================================="
echo "STAGE 1: BWA Mapping to All Genes"
echo "==================================================================="

# Extract Pachyramphus sequences from all aligned genes
echo "Extracting Pachyramphus reference sequences..."

if [ ! -f "${STAGE1_DIR}/pachyramphus_all.fasta" ]; then
    for aligned in ${WD}/results/graph_alignment/aligned/*.fasta; do
        gene=$(basename ${aligned} .fasta)
        # Extract only Pachyramphus sequence
        grep -A 1000000 ">ref1_Pachyramphus_minor" ${aligned} | \
            sed '/^>ref[2-9]/,$d' | \
            sed "s/>ref1_Pachyramphus_minor.*/>ref1_Pachyramphus_minor#${gene}/" \
            >> ${STAGE1_DIR}/pachyramphus_all.fasta
    done
    echo "Extracted $(grep -c '^>' ${STAGE1_DIR}/pachyramphus_all.fasta) genes"
else
    echo "Reference already exists"
fi

# Index reference
echo "Indexing BWA reference..."
if [ ! -f "${STAGE1_DIR}/pachyramphus_all.fasta.bwt" ]; then
    ${BWA} index ${STAGE1_DIR}/pachyramphus_all.fasta
fi

${SAMTOOLS} faidx ${STAGE1_DIR}/pachyramphus_all.fasta

# Map reads
echo "Mapping paired-end reads with BWA..."
if [ ! -f "${STAGE1_DIR}/mapped_pe.bam" ]; then
    ${BWA} mem -t ${THREADS} ${STAGE1_DIR}/pachyramphus_all.fasta \
        ${PE1} ${PE2} | \
        ${SAMTOOLS} view -@ ${THREADS} -bS - | \
        ${SAMTOOLS} sort -@ ${THREADS} -o ${STAGE1_DIR}/mapped_pe.bam -
    ${SAMTOOLS} index ${STAGE1_DIR}/mapped_pe.bam
fi

echo "Mapping merged reads with BWA..."
if [ ! -f "${STAGE1_DIR}/mapped_merged.bam" ]; then
    ${BWA} mem -t ${THREADS} ${STAGE1_DIR}/pachyramphus_all.fasta \
        ${MERGED} | \
        ${SAMTOOLS} view -@ ${THREADS} -bS - | \
        ${SAMTOOLS} sort -@ ${THREADS} -o ${STAGE1_DIR}/mapped_merged.bam -
    ${SAMTOOLS} index ${STAGE1_DIR}/mapped_merged.bam
fi

# Merge BAMs
echo "Merging BAM files..."
if [ ! -f "${STAGE1_DIR}/all_reads.bam" ]; then
    ${SAMTOOLS} merge -@ ${THREADS} ${STAGE1_DIR}/all_reads.bam \
        ${STAGE1_DIR}/mapped_pe.bam ${STAGE1_DIR}/mapped_merged.bam
    ${SAMTOOLS} index ${STAGE1_DIR}/all_reads.bam
fi

# Calculate coverage per gene
echo "Calculating coverage per gene..."

${SAMTOOLS} depth ${STAGE1_DIR}/all_reads.bam | \
    awk '{gene=$1; sub(/#.*/, "", gene); cov[gene]+=$3; len[gene]++} 
         END {for(g in cov) print g"\t"cov[g]/len[g]"\t"len[g]}' | \
    sort -k2,2nr > ${STAGE1_DIR}/gene_coverage.txt

echo "Coverage calculated for $(wc -l < ${STAGE1_DIR}/gene_coverage.txt) genes"

# Select top N genes
echo "Selecting top ${TOP_N} genes by coverage..."

head -n ${TOP_N} ${STAGE1_DIR}/gene_coverage.txt > ${STAGE1_DIR}/top_genes.txt

echo ""
echo "Top 10 genes by coverage:"
head -n 10 ${STAGE1_DIR}/top_genes.txt | \
    awk '{printf "  %s: %.2fx coverage (%d bp)\n", $1, $2, $3}'

###############################################################################
# STAGE 2: Concatenated Graph Mapping for Top Genes
###############################################################################

echo ""
echo "==================================================================="
echo "STAGE 2: Graph Mapping for Top ${TOP_N} Genes"
echo "==================================================================="

# Extract aligned sequences for top genes only
echo "Concatenating top ${TOP_N} genes..."

rm -f ${STAGE2_DIR}/top_genes.fasta

while read gene avg_cov length; do
    gene_file="${WD}/results/graph_alignment/aligned/${gene}.fasta"
    if [ -f "${gene_file}" ]; then
        cat ${gene_file} >> ${STAGE2_DIR}/top_genes.fasta
    fi
done < ${STAGE1_DIR}/top_genes.txt

num_seqs=$(grep -c '^>' ${STAGE2_DIR}/top_genes.fasta || true)
echo "Concatenated ${num_seqs} sequences"

# Build variation graph
echo "Building variation graph..."
${VG} construct -M ${STAGE2_DIR}/top_genes.fasta -m 32 -t ${THREADS} \
    > ${STAGE2_DIR}/top.vg

echo "Indexing XG..."
${VG} index -x ${STAGE2_DIR}/top.xg ${STAGE2_DIR}/top.vg -t ${THREADS}
rm ${STAGE2_DIR}/top.vg

echo "Pruning graph..."
${VG} prune ${STAGE2_DIR}/top.xg -k 16 -t ${THREADS} > ${STAGE2_DIR}/top.pruned.vg

echo "Building GCSA index..."
${VG} index -g ${STAGE2_DIR}/top.gcsa -k 16 ${STAGE2_DIR}/top.pruned.vg -t ${THREADS}
rm ${STAGE2_DIR}/top.pruned.vg

# Map reads ONCE to concatenated graph
echo "Mapping paired-end reads to graph..."
${VG} map -x ${STAGE2_DIR}/top.xg -g ${STAGE2_DIR}/top.gcsa \
    -f ${PE1} -f ${PE2} -t ${THREADS} > ${STAGE2_DIR}/top.gam

echo "Appending merged reads..."
${VG} map -x ${STAGE2_DIR}/top.xg -g ${STAGE2_DIR}/top.gcsa \
    -f ${MERGED} -t ${THREADS} >> ${STAGE2_DIR}/top.gam

rm -f ${STAGE2_DIR}/top.gcsa*

# Call variants
echo "Packing read coverage..."
${VG} pack -x ${STAGE2_DIR}/top.xg -g ${STAGE2_DIR}/top.gam \
    -o ${STAGE2_DIR}/top.pack -t ${THREADS}

rm ${STAGE2_DIR}/top.gam

echo "Calling variants..."
${VG} call ${STAGE2_DIR}/top.xg -k ${STAGE2_DIR}/top.pack \
    -d 1 -s 1 -t ${THREADS} > ${STAGE2_DIR}/top_raw.vcf

rm ${STAGE2_DIR}/top.pack

echo "Filtering variants (DP>${MIN_DEPTH})..."
${BCFTOOLS} view -i "FORMAT/DP>${MIN_DEPTH}" ${STAGE2_DIR}/top_raw.vcf \
    -o ${STAGE2_DIR}/top_filtered.vcf

rm ${STAGE2_DIR}/top_raw.vcf

echo "Compressing and indexing VCF..."
${BGZIP} -c ${STAGE2_DIR}/top_filtered.vcf > ${STAGE2_DIR}/top.vcf.gz
${TABIX} -p vcf ${STAGE2_DIR}/top.vcf.gz

rm ${STAGE2_DIR}/top_filtered.vcf

###############################################################################
# Generate consensus for each gene
###############################################################################

echo ""
echo "==================================================================="
echo "Generating Consensus Sequences"
echo "==================================================================="

PATHS=$(${BCFTOOLS} query -l ${STAGE2_DIR}/top.vcf.gz | grep "ref1_Pachyramphus_minor" || true)

processed=0
skipped=0

for path in ${PATHS}; do
    # Extract gene name
    GENE=$(echo ${path} | sed 's/^ref1_Pachyramphus_minor#0#//' | sed 's/:.*$//')
    
    echo "Processing ${GENE}..."
    
    # Extract reference
    ${SAMTOOLS} faidx ${STAGE2_DIR}/top.xg "${path}" > ${STAGE2_DIR}/${GENE}_ref.fasta 2>/dev/null || {
        echo "  Skipping ${GENE}: failed to extract reference"
        ((skipped++))
        continue
    }
    
    # Extract variants for this gene
    ${BCFTOOLS} view -r "${path}" ${STAGE2_DIR}/top.vcf.gz > ${STAGE2_DIR}/${GENE}.vcf
    
    # Check if any variants
    variant_count=$(grep -v '^#' ${STAGE2_DIR}/${GENE}.vcf | wc -l)
    if [ ${variant_count} -eq 0 ]; then
        echo "  Skipping ${GENE}: no variants passing filters"
        rm -f ${STAGE2_DIR}/${GENE}*
        ((skipped++))
        continue
    fi
    
    # Check average depth
    AVG_DEPTH=$(${BCFTOOLS} query -f '%DP\n' ${STAGE2_DIR}/${GENE}.vcf | \
        awk '{sum+=$1; n++} END {if(n>0) print sum/n; else print 0}')
    
    if (( $(echo "${AVG_DEPTH} < ${MIN_DEPTH}" | bc -l) )); then
        echo "  Skipping ${GENE}: average depth ${AVG_DEPTH}x < ${MIN_DEPTH}x"
        rm -f ${STAGE2_DIR}/${GENE}*
        ((skipped++))
        continue
    fi
    
    # Compress and index
    ${BGZIP} -c ${STAGE2_DIR}/${GENE}.vcf > ${STAGE2_DIR}/${GENE}.vcf.gz
    ${TABIX} -p vcf ${STAGE2_DIR}/${GENE}.vcf.gz
    
    # Generate consensus
    ${BCFTOOLS} consensus -f ${STAGE2_DIR}/${GENE}_ref.fasta \
        ${STAGE2_DIR}/${GENE}.vcf.gz > ${STAGE2_DIR}/${GENE}_consensus.fasta
    
    # Rename header
    sed "s/^>.*/>Tityra_leucura_graph/" ${STAGE2_DIR}/${GENE}_consensus.fasta \
        > ${WD}/results/graph_alignment/consensus/${GENE}_Tityra.fasta
    
    echo "  Success: ${variant_count} variants, depth ${AVG_DEPTH}x"
    ((processed++))
    
    # Clean up
    rm -f ${STAGE2_DIR}/${GENE}*
done

echo ""
echo "Processed: ${processed} genes"
echo "Skipped: ${skipped} genes"

###############################################################################
# Translate to proteins
###############################################################################

echo ""
echo "==================================================================="
echo "Translating to Proteins"
echo "==================================================================="

python3 << 'PYTHON_SCRIPT'
import os
from Bio import SeqIO
from Bio.Seq import Seq

wd = os.environ['WD']
consensus_dir = f"{wd}/results/graph_alignment/consensus"

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
        protein_seq = record.seq.translate()
        
        if '*' in str(protein_seq)[:-1]:
            print(f"Skipping {gene_name}: internal stop codon")
            skipped += 1
            continue
        
        with open(protein_file, 'w') as f:
            f.write(f">{record.id}\n")
            f.write(str(protein_seq) + '\n')
        
        processed += 1
        
    except Exception as e:
        print(f"Error processing {gene_name}: {e}")
        skipped += 1

print(f"Translation complete: {processed} genes, {skipped} skipped")
PYTHON_SCRIPT

###############################################################################
# Final summary
###############################################################################

echo ""
echo "==================================================================="
echo "Pipeline Complete!"
echo "==================================================================="

CONSENSUS_COUNT=$(ls ${WD}/results/graph_alignment/consensus/*_Tityra.fasta 2>/dev/null | wc -l)
PROTEIN_COUNT=$(ls ${WD}/results/graph_alignment/consensus/*_protein.fasta 2>/dev/null | wc -l)

echo "Stage 1: Mapped to $(wc -l < ${STAGE1_DIR}/gene_coverage.txt) genes"
echo "Stage 2: Graph mapping for top ${TOP_N} genes"
echo "Final consensus sequences: ${CONSENSUS_COUNT}"
echo "Final protein sequences: ${PROTEIN_COUNT}"
echo ""
echo "Results in: ${WD}/results/graph_alignment/consensus/"
echo ""
echo "Top genes list: ${STAGE1_DIR}/top_genes.txt"
