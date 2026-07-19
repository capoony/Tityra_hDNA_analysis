#!/bin/bash
set -euo pipefail
###############################################################################
# Tityra UCE Mapping & Phylogeny Pipeline
# Run each section manually in order.
###############################################################################

###############################################################################
# Variables
###############################################################################

readonly WD="/media/inter/mkapun/projects/Tityra"
readonly DATA="${WD}/data"
readonly OUT="${WD}/results/tityra_uce"
readonly UCE_REF="${DATA}/Tityra_concat_UCEs.fasta"
readonly UCE_DIR="${DATA}/openwings_tityrinae"

source /opt/anaconda3/etc/profile.d/conda.sh

###############################################################################
# 1. Download Luke's Tityrinae UCE alignment archive
###############################################################################

mkdir -p "${UCE_DIR}"
cd "${UCE_DIR}"

wget -O Tityrinae-mafft-nexus-clean-54taxa-95per.zip \
  "https://dryad-assetstore-merritt-west.s3.us-west-2.amazonaws.com/ark%3A/13030/m5v45x0k%7C1%7Cproducer/Tityrinae-mafft-nexus-clean-54taxa-95per.zip?response-content-disposition=attachment%3B%20filename%3DTityrinae-mafft-nexus-clean-54taxa-95per.zip&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIA2KERHV5E3OITXZXC%2F20260709%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Date=20260709T133013Z&X-Amz-Expires=86400&X-Amz-SignedHeaders=host&X-Amz-Signature=5e97713c9ba4803967341b740d40e1b7726d489d4701419e98120168a725397c"

unzip -o Tityrinae-mafft-nexus-clean-54taxa-95per.zip

###############################################################################
# 2. Build a Tityra-specific UCE multi-fasta reference
###############################################################################

python "${WD}/scripts/extract_longest_tityra_uce.py"

###############################################################################
# 3. Map reads to the Tityra UCE reference with bwa-mem2
###############################################################################

mkdir -p "${OUT}/mapping"
conda activate bwa-mem2

bwa-mem2 index "${UCE_REF}"

bwa-mem2 mem -t 4 \
  "${UCE_REF}" \
  "${DATA}/trimmed2/Tityra_leucura_contaminants_merged_unmapped.fastq.gz" \
  2> "${OUT}/mapping/bwa_mem2.log" | \
  samtools view -bS -F 4 \
  -o "${OUT}/mapping/Tityra_leucura_uce_mapped.bam" -

samtools sort -o "${OUT}/mapping/Tityra_leucura_uce_mapped_sorted.bam" \
  "${OUT}/mapping/Tityra_leucura_uce_mapped.bam"
samtools index "${OUT}/mapping/Tityra_leucura_uce_mapped_sorted.bam"

samtools flagstat \
  "${OUT}/mapping/Tityra_leucura_uce_mapped.bam" \
  > "${OUT}/mapping/Tityra_leucura_uce_flagstat.txt"

samtools coverage \
  "${OUT}/mapping/Tityra_leucura_uce_mapped_sorted.bam" \
  > "${OUT}/mapping/Tityra_leucura_uce_coverage.txt"

Rscript "${WD}/scripts/analyze_coverage.R"

samtools consensus -a -f fasta \
  -o "${OUT}/mapping/Tityra_leucura_uce_consensus.fasta" \
  "${OUT}/mapping/Tityra_leucura_uce_mapped_sorted.bam"

conda deactivate

###############################################################################
# 4. Build per-locus FASTA inputs for UCE phylogeny
###############################################################################

python "${WD}/scripts/create_uce_phylo_pipeline.py"

###############################################################################
# 5. Align each locus with MAFFT in parallel
###############################################################################

mkdir -p "${OUT}/phylogeny/alignments"
conda activate mafft-7.487

export WD
align_locus() {
  locus=$1
  in="${WD}/results/tityra_uce/phylogeny/gene_sequences/${locus}.fasta"
  out="${WD}/results/tityra_uce/phylogeny/alignments/${locus}_aln.fasta"
  if [ ! -f "$out" ]; then
    mafft --auto --adjustdirection --quiet "$in" 2>/dev/null | \
      sed 's/^>_R_/>/g' > "$out"
    echo "  Aligned $locus ($(grep -c '^>' "$out") seqs)"
  fi
}
export -f align_locus

cat "${WD}/results/tityra_uce/phylogeny/loci_list.txt" | \
  /usr/bin/parallel -j 24 --will-cite align_locus {}

NUM_ALN=$(find "${WD}/results/tityra_uce/phylogeny/alignments/" -name "*_aln.fasta" | wc -l)
echo "Alignments produced: ${NUM_ALN}"
conda deactivate

###############################################################################
# 6. Build the supermatrix and partition file
###############################################################################

python "${WD}/scripts/build_supermatrix_and_clean.py"

###############################################################################
# 7. Run IQ-TREE partitioned analysis on the full supermatrix
###############################################################################

mkdir -p "${OUT}/phylogeny/trees" \
         "${OUT}/phylogeny/plots"

conda activate iq-tree-3.0.1

OUTGROUP="acachl"
SUPERMATRIX="${WD}/results/tityra_uce/phylogeny/supermatrix/supermatrix.fasta"
PARTITIONS="${WD}/results/tityra_uce/phylogeny/supermatrix/partitions.txt"
PREFIX="${WD}/results/tityra_uce/phylogeny/trees/uce_partitioned"

iqtree \
  -s "${SUPERMATRIX}" \
  -p "${PARTITIONS}" \
  -m MFP+MERGE \
  -bb 1000 \
  -alrt 1000 \
  -T 25 \
  --prefix "${PREFIX}" \
  --redo

Rscript "${WD}/scripts/plot_uce_trees.R" \
  "${PREFIX}" \
  "${WD}/results/tityra_uce/phylogeny/plots" \
  "${OUTGROUP}"

conda deactivate

###############################################################################
# 8. Alternative analysis: trim N-columns from Tityra_leucura and rerun IQ-TREE
###############################################################################

python "${WD}/scripts/trim_N_columns_and_build.py"

conda activate iq-tree-3.0.1

SUPERMATRIX_TRIM="${WD}/results/tityra_uce/phylogeny/trimmed/supermatrix/supermatrix_trimmed.fasta"
PARTITIONS_TRIM="${WD}/results/tityra_uce/phylogeny/trimmed/supermatrix/partitions_trimmed.txt"
PREFIX_TRIM="${WD}/results/tityra_uce/phylogeny/trees/uce_partitioned_trimmed"

iqtree \
  -s "${SUPERMATRIX_TRIM}" \
  -p "${PARTITIONS_TRIM}" \
  -m MFP+MERGE \
  -bb 1000 \
  -alrt 1000 \
  -T 25 \
  --prefix "${PREFIX_TRIM}" \
  --redo

Rscript "${WD}/scripts/plot_uce_trees.R" \
  "${PREFIX_TRIM}" \
  "${WD}/results/tityra_uce/phylogeny/plots" \
  "${OUTGROUP}"

conda deactivate

###############################################################################
# 9. Summary
###############################################################################

echo ""
echo "==================================================================="
echo "UCE phylogeny pipeline complete!"
echo ""
echo "=== Analysis 1 (full alignments) ==="
echo "  Tree:   ${PREFIX}.treefile"
echo "  Log:    ${PREFIX}.log"
echo "  Matrix: ${SUPERMATRIX}"
echo "  Length: 822,495 nt"
echo ""
echo "=== Analysis 2 (N-columns trimmed from Tityra_leucura) ==="
echo "  Tree:   ${PREFIX_TRIM}.treefile"
echo "  Log:    ${PREFIX_TRIM}.log"
echo "  Matrix: ${SUPERMATRIX_TRIM}"
echo "  Length: 358,055 nt"
echo ""
echo "=== Plots ==="
echo "  ${WD}/results/tityra_uce/phylogeny/plots/"
echo "==================================================================="
