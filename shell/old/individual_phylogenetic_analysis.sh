#!/bin/bash
###############################################################################
# Individual-Based Phylogenetic Analysis of Concatenated Mitochondrial Genes
# - Aligns concatenated sequences with MAFFT
# - Builds phylogenetic tree with IQ-TREE (using partitioned model)
# - Plots tree with R/ggplot2
###############################################################################

# Set working directory
WD=${WD:-/media/inter/mkapun/projects/Tityra}

# Input/Output directories
INPUT_DIR="${WD}/results/individual_concatenated"
OUTPUT_DIR="${WD}/results/individual_phylogenetic"
ALIGNMENT_DIR="${OUTPUT_DIR}/alignment"
TREE_DIR="${OUTPUT_DIR}/tree"
PLOT_DIR="${OUTPUT_DIR}/plots"

# Create output directories
mkdir -p ${ALIGNMENT_DIR}
mkdir -p ${TREE_DIR}
mkdir -p ${PLOT_DIR}

echo "==================================================================="
echo "Individual-Based Phylogenetic Analysis Pipeline"
echo "Started: $(date)"
echo "==================================================================="
echo ""
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

###############################################################################
# 1. Multiple Sequence Alignment with MAFFT
###############################################################################
echo "==================================================================="
echo "Step 1: Multiple Sequence Alignment with MAFFT"
echo "==================================================================="

# Activate MAFFT environment
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate mafft-7.487

INPUT_FILE="${INPUT_DIR}/concatenated_sequences.fasta"
OUTPUT_FILE="${ALIGNMENT_DIR}/concatenated_aligned.fasta"

if [ -f "${INPUT_FILE}" ]; then
    echo ""
    echo "Aligning concatenated sequences..."
    echo "Input: ${INPUT_FILE}"
    
    # Count sequences
    NUM_SEQS=$(grep -c "^>" "${INPUT_FILE}")
    echo "Number of sequences: ${NUM_SEQS}"
    
    # Run MAFFT with auto mode
    # Using --auto for optimal algorithm selection based on dataset size
    # --adjustdirection: automatically reverse complement sequences if needed
    # --thread 8 for faster processing
    mafft --auto \
          --adjustdirection \
          --thread 8 \
          "${INPUT_FILE}" > "${OUTPUT_FILE}"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Alignment completed"
        echo "  Output: ${OUTPUT_FILE}"
        
        # Get alignment statistics
        SEQ_LENGTH=$(grep -v "^>" "${OUTPUT_FILE}" | head -1 | wc -c)
        echo "  Alignment length: ${SEQ_LENGTH} bp"
    else
        echo "  ✗ Error: Alignment failed"
        exit 1
    fi
else
    echo "  ✗ Error: Input file not found: ${INPUT_FILE}"
    exit 1
fi

conda deactivate
echo ""

###############################################################################
# 2. Remove Gappy Regions with Gblocks
###############################################################################
echo "==================================================================="
echo "Step 2: Removing Gappy Regions with Gblocks"
echo "==================================================================="

# Load Gblocks module
module load Alignment/Gblocks-0.91b

ALIGNMENT_FILE="${ALIGNMENT_DIR}/concatenated_aligned.fasta"
GBLOCKS_OUTPUT="${ALIGNMENT_DIR}/concatenated_aligned_gblocks.fasta"

if [ -f "${ALIGNMENT_FILE}" ]; then
    echo ""
    echo "Processing concatenated alignment with Gblocks..."
    
    # Gblocks has issues with long sequence names, so create short-named version
    ALIGNMENT_SHORT="${ALIGNMENT_DIR}/concatenated_aligned_short.fasta"
    NAME_MAP="${ALIGNMENT_DIR}/name_map.tsv"
    
    # Create short-named version and save name mapping
    echo "  Creating short-named version for Gblocks..."
    awk '/^>/ {print ">seq_" ++i; next} {print}' "${ALIGNMENT_FILE}" > "${ALIGNMENT_SHORT}"
    awk 'BEGIN{i=0} /^>/{i++; names[i]=$0; next} END{for(j=1;j<=i;j++) print j"\t>seq_"j"\t"names[j]}' "${ALIGNMENT_FILE}" > "${NAME_MAP}"
    
    # Run Gblocks on short-named version
    # -t=d: DNA sequences
    # -b1=0.5: minimum sequences for conserved position (50%)
    # -b2=0.85: minimum sequences for flanking position (85%)
    # -b3=8: maximum contiguous nonconserved positions
    # -b4=10: minimum block length
    # -b5=h: allow gap positions within final blocks (half)
    Gblocks "${ALIGNMENT_SHORT}" -t=d -b1=0.5 -b2=0.85 -b3=8 -b4=10 -b5=h
    
    # Gblocks creates output with -gb extension
    GB_OUTPUT="${ALIGNMENT_SHORT}-gb"
    
    if [ -f "${GB_OUTPUT}" ]; then
        # Restore original sequence names
        echo "  Restoring original sequence names..."
        awk 'NR==FNR {map[$2]=$3; next} /^>/ {print map[$0]; next} {print}' "${NAME_MAP}" "${GB_OUTPUT}" > "${GBLOCKS_OUTPUT}"
        
        echo "  ✓ Gblocks processing completed"
        echo "  Output: ${GBLOCKS_OUTPUT}"
        
        # Get Gblocks statistics
        ORIG_LENGTH=$(grep -v "^>" "${ALIGNMENT_FILE}" | tr -d '\n' | wc -c)
        GBLOCKS_LENGTH=$(grep -v "^>" "${GBLOCKS_OUTPUT}" | tr -d '\n' | wc -c)
        echo "  Original alignment length: ${ORIG_LENGTH} bp"
        echo "  Gblocks filtered length: ${GBLOCKS_LENGTH} bp"
        echo "  Retained: $(echo "scale=2; ${GBLOCKS_LENGTH}/${ORIG_LENGTH}*100" | bc)%"
        
        # Clean up temporary files
        rm -f "${ALIGNMENT_SHORT}" "${GB_OUTPUT}"
    else
        echo "  ⚠ Warning: Gblocks output not found, using original alignment"
        cp "${ALIGNMENT_FILE}" "${GBLOCKS_OUTPUT}"
    fi
else
    echo "  ✗ Error: Alignment file not found: ${ALIGNMENT_FILE}"
    exit 1
fi

echo ""

###############################################################################
# 3. Prepare Partition File for IQ-TREE
###############################################################################
echo "==================================================================="
echo "Step 3: Preparing Partition File"
echo "==================================================================="

# Copy partition file from individual_concatenated results
PARTITION_INPUT="${INPUT_DIR}/gene_partitions.txt"
PARTITION_OUTPUT="${ALIGNMENT_DIR}/partitions.nex"

if [ -f "${PARTITION_INPUT}" ]; then
    cp "${PARTITION_INPUT}" "${PARTITION_OUTPUT}"
    echo "  ✓ Partition file copied"
    echo "  Partitions: $(wc -l < ${PARTITION_OUTPUT})"
else
    echo "  ⚠ Warning: Partition file not found, will use unpartitioned model"
    PARTITION_OUTPUT=""
fi

echo ""

###############################################################################
# 4. Phylogenetic Tree Reconstruction with IQ-TREE
###############################################################################
echo "==================================================================="
echo "Step 4: Phylogenetic Tree Reconstruction with IQ-TREE"
echo "==================================================================="

# Activate IQ-TREE environment
conda activate iq-tree-2.3.3

echo ""
echo "Building concatenated tree..."

# Change to tree directory
cd ${TREE_DIR}

# Use Gblocks-filtered alignment
FINAL_ALIGNMENT="${ALIGNMENT_DIR}/concatenated_aligned_gblocks.fasta"

# Determine IQ-TREE command based on partition file availability
if [ -n "${PARTITION_OUTPUT}" ] && [ -f "${PARTITION_OUTPUT}" ]; then
    echo "Using partitioned model (one model per gene)"
    
    # Run IQ-TREE with partitioned analysis:
    # -s: alignment file
    # -p: partition file
    # -m MFP: find best-fit model for each partition
    # -bb 1000: ultrafast bootstrap (fast approximation)
    # -nt 8: use 8 threads (AUTO can fail with small datasets)
    # -pre: output prefix
    iqtree -s "${FINAL_ALIGNMENT}" \
           -p "${PARTITION_OUTPUT}" \
           -m MFP \
           -bb 1000 \
           -nt 8 \
           -pre individual_concatenated
else
    echo "Using unpartitioned model"
    
    # Run IQ-TREE without partitions:
    iqtree -s "${FINAL_ALIGNMENT}" \
           -m MFP \
           -bb 1000 \
           -nt 8 \
           -pre individual_concatenated
fi

if [ $? -eq 0 ]; then
    echo "  ✓ Tree reconstruction completed"
    echo "  Tree file: ${TREE_DIR}/individual_concatenated.treefile"
else
    echo "  ✗ Error: Tree reconstruction failed"
    conda deactivate
    exit 1
fi

conda deactivate
echo ""

###############################################################################
# 5. Tree Visualization with R
###############################################################################
echo "==================================================================="
echo "Step 5: Tree Visualization with R"
echo "==================================================================="

# Activate R environment with required packages
conda activate r-ggtree-4.3

# Create R script for tree plotting
cat > ${PLOT_DIR}/plot_individual_tree.R << 'EOF'
#!/usr/bin/env Rscript
###############################################################################
# Plot Individual-Based Phylogenetic Tree
###############################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
  library(treeio)
  library(dplyr)
})

# Set paths
tree_file <- "TREE_FILE_PLACEHOLDER"
output_prefix <- "OUTPUT_PREFIX_PLACEHOLDER"

# Read tree
cat("Reading tree from:", tree_file, "\n")
tree <- read.iqtree(tree_file)

# Extract species names from tip labels
# Format: Genus_species_voucher
tip_data <- data.frame(
  label = tree@phylo$tip.label,
  stringsAsFactors = FALSE
) %>%
  mutate(
    # Extract species name (first two parts before voucher)
    species = gsub("^([^_]+_[^_]+)_.*", "\\1", label),
    species = gsub("_", " ", species),
    # Extract genus for coloring
    genus = gsub("^([^_]+)_.*", "\\1", label)
  )

# Count unique species
num_species <- length(unique(tip_data$species))
num_specimens <- nrow(tip_data)

cat("Number of specimens:", num_specimens, "\n")
cat("Number of species:", num_species, "\n")

###############################################################################
# Plot 1: Rectangular tree with bootstrap values
###############################################################################
cat("\nGenerating Plot 1: Rectangular tree with bootstrap support...\n")

p1 <- ggtree(tree, layout="rectangular") +
  geom_tiplab(aes(label=label), size=3, align=FALSE) +
  geom_nodelab(aes(label=label, subset=!is.na(as.numeric(label)) & as.numeric(label) >= 70),
               size=2.5, color="red", nudge_x=-0.002, nudge_y=0.3) +
  theme_tree2() +
  xlim(0, max(tree@data$x) * 1.4) +
  ggtitle("Individual-Based Mitochondrial Phylogeny",
          subtitle=paste0(num_specimens, " specimens, ", num_species, " species | Bootstrap support ≥70 shown"))

ggsave(paste0(output_prefix, "_rectangular_bootstrap.pdf"), 
       p1, width=12, height=num_specimens*0.3 + 2, limitsize=FALSE)

###############################################################################
# Plot 2: Circular tree
###############################################################################
cat("Generating Plot 2: Circular tree...\n")

p2 <- ggtree(tree, layout="circular") +
  geom_tiplab(aes(label=label), size=2.5, offset=0.001) +
  ggtitle("Individual-Based Mitochondrial Phylogeny (Circular)",
          subtitle=paste0(num_specimens, " specimens, ", num_species, " species"))

ggsave(paste0(output_prefix, "_circular.pdf"), 
       p2, width=14, height=14, limitsize=FALSE)

###############################################################################
# Plot 3: Rectangular tree colored by genus
###############################################################################
cat("Generating Plot 3: Tree colored by genus...\n")

# Merge tip data with tree
tree_data <- tree@phylo
tip_info <- tip_data %>% select(label, genus, species)

p3 <- ggtree(tree, layout="rectangular") %<+% tip_info +
  geom_tiplab(aes(label=label, color=genus), size=3, align=FALSE) +
  geom_nodelab(aes(label=label, subset=!is.na(as.numeric(label)) & as.numeric(label) >= 70),
               size=2.5, color="red", nudge_x=-0.002, nudge_y=0.3) +
  theme_tree2() +
  xlim(0, max(tree@data$x) * 1.4) +
  theme(legend.position="right") +
  ggtitle("Individual-Based Mitochondrial Phylogeny (Colored by Genus)",
          subtitle=paste0(num_specimens, " specimens, ", num_species, " species"))

ggsave(paste0(output_prefix, "_rectangular_genus.pdf"), 
       p3, width=14, height=num_specimens*0.3 + 2, limitsize=FALSE)

###############################################################################
# Plot 4: Fan/radial tree
###############################################################################
cat("Generating Plot 4: Fan tree...\n")

p4 <- ggtree(tree, layout="fan", open.angle=30) +
  geom_tiplab(aes(label=label), size=2.5, offset=0.001) +
  ggtitle("Individual-Based Mitochondrial Phylogeny (Fan)",
          subtitle=paste0(num_specimens, " specimens, ", num_species, " species"))

ggsave(paste0(output_prefix, "_fan.pdf"), 
       p4, width=14, height=14, limitsize=FALSE)

cat("\nAll plots saved successfully!\n")
cat("Output directory:", dirname(output_prefix), "\n")
EOF

# Replace placeholders in R script
sed -i "s|TREE_FILE_PLACEHOLDER|${TREE_DIR}/individual_concatenated.treefile|g" ${PLOT_DIR}/plot_individual_tree.R
sed -i "s|OUTPUT_PREFIX_PLACEHOLDER|${PLOT_DIR}/individual_tree|g" ${PLOT_DIR}/plot_individual_tree.R

# Run R script
echo ""
echo "Running R visualization script..."
Rscript ${PLOT_DIR}/plot_individual_tree.R

if [ $? -eq 0 ]; then
    echo "  ✓ Tree visualization completed"
    echo "  Plots saved to: ${PLOT_DIR}/"
else
    echo "  ✗ Warning: Tree visualization failed"
fi

conda deactivate

###############################################################################
# Summary
###############################################################################
echo ""
echo "==================================================================="
echo "Pipeline completed: $(date)"
echo "==================================================================="
echo ""
echo "Output files:"
echo "  - Alignment: ${ALIGNMENT_DIR}/concatenated_aligned.fasta"
echo "  - Gblocks filtered: ${ALIGNMENT_DIR}/concatenated_aligned_gblocks.fasta"
if [ -n "${PARTITION_OUTPUT}" ]; then
    echo "  - Partitions: ${PARTITION_OUTPUT}"
fi
echo "  - Tree file: ${TREE_DIR}/individual_concatenated.treefile"
echo "  - Log file: ${TREE_DIR}/individual_concatenated.log"
echo "  - Plots: ${PLOT_DIR}/"
echo ""
echo "PDF plots generated:"
for pdf in ${PLOT_DIR}/*.pdf; do
    if [ -f "$pdf" ]; then
        echo "  - $(basename $pdf)"
    fi
done
echo ""
