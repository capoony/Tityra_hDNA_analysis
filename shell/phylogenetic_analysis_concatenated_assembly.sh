#!/bin/bash

###############################################################################
# Phylogenetic Analysis of Concatenated Mitochondrial Genes (Assembly-based)
# Description: Perform phylogenetic inference on concatenated gene alignments
#              using partitioned analysis in IQ-TREE
###############################################################################

set -euo pipefail

# Set working directory
WD=/media/inter/mkapun/projects/Tityra

# Input/Output directories
INPUT_DIR="${WD}/results/concatenated_genes_assembly"
OUTPUT_DIR="${WD}/results/phylogenetic_analysis_concatenated_assembly"
TREE_DIR="${OUTPUT_DIR}/trees"
PLOT_DIR="${OUTPUT_DIR}/plots"

# Create output directories
mkdir -p "${TREE_DIR}"
mkdir -p "${PLOT_DIR}"

# Log file
LOG="${WD}/logs/phylogenetic_analysis_concatenated_assembly.log"
mkdir -p "${WD}/logs"

exec > >(tee -a "${LOG}") 2>&1

echo "==================================================================="
echo "Phylogenetic Analysis of Concatenated Genes (Assembly-based)"
echo "Started: $(date)"
echo "==================================================================="
echo ""
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Input files
CONCATENATED_FILE="${INPUT_DIR}/concatenated_all_genes.fasta"
PARTITION_FILE="${INPUT_DIR}/partitions.txt"

# Check input files
if [ ! -f "${CONCATENATED_FILE}" ]; then
    echo "Error: Concatenated sequences file not found: ${CONCATENATED_FILE}"
    exit 1
fi

if [ ! -f "${PARTITION_FILE}" ]; then
    echo "Error: Partition file not found: ${PARTITION_FILE}"
    exit 1
fi

echo "Using concatenated alignment: ${CONCATENATED_FILE}"
echo "Using partitions: ${PARTITION_FILE}"
echo ""

###############################################################################
# 1. Run IQ-TREE with Partitioned Analysis
###############################################################################
echo "==================================================================="
echo "Step 1: Phylogenetic Inference with IQ-TREE (Partitioned)"
echo "==================================================================="

# Activate IQ-TREE environment
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate iq-tree-2.3.3

cd "${TREE_DIR}"

# Check if tree already exists
if [ -f "concatenated_all_genes.treefile" ]; then
    echo "Tree file already exists, skipping IQ-TREE analysis..."
    echo "  Using existing tree: concatenated_all_genes.treefile"
    echo "  (Delete this file if you want to rerun the analysis)"
    echo ""
else
    echo "Running IQ-TREE with partitioned analysis..."
    echo "This may take a while..."
    echo ""

    # Run IQ-TREE with:
    # -s: alignment file
    # -p: partition file
    # -m MFP: model finder plus (automatic model selection per partition)
    # -bb 1000: ultrafast bootstrap with 1000 replicates
    # -alrt 1000: SH-aLRT test with 1000 replicates
    # -nt AUTO: automatic thread detection
    # -pre: prefix for output files

    iqtree2 \
        -s "${CONCATENATED_FILE}" \
        -p "${PARTITION_FILE}" \
        -m MFP \
        -bb 1000 \
        -alrt 1000 \
        -nt AUTO \
        -pre concatenated_all_genes

    if [ $? -eq 0 ]; then
        echo "  ✓ IQ-TREE analysis completed successfully"
    else
        echo "  ✗ Error: IQ-TREE analysis failed"
        exit 1
    fi
fi

echo ""
echo "Tree files generated:"
ls -lh concatenated_all_genes.treefile
echo ""

# Return to original directory
cd "${WD}"

###############################################################################
# 2. Plot Phylogenetic Trees with R
###############################################################################
echo "==================================================================="
echo "Step 2: Plotting Phylogenetic Trees"
echo "==================================================================="

# Deactivate iq-tree conda environment
conda deactivate

# Export directories for R script
export TREE_DIR
export PLOT_DIR
export TREE_FILE="${TREE_DIR}/concatenated_all_genes.treefile"

# Run R script for tree plotting
Rscript - << 'R_SCRIPT'

# Load required libraries
suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
})

# Get tree file from environment
tree_file <- Sys.getenv("TREE_FILE")
plot_dir <- Sys.getenv("PLOT_DIR")

if (!file.exists(tree_file)) {
  stop("Tree file not found: ", tree_file)
}

cat("Reading tree file:", tree_file, "\n")

# Read tree
tree <- read.tree(tree_file)

cat("Tree successfully loaded\n")
cat("Number of tips:", length(tree$tip.label), "\n")
cat("Tips:", paste(tree$tip.label, collapse=", "), "\n\n")

# Root tree with Onychorhynchus as outgroup
onychorhynchus_tips <- grep("^Onychorhynchus", tree$tip.label, value = TRUE)

if (length(onychorhynchus_tips) > 0) {
  tree <- root(tree, outgroup = onychorhynchus_tips, resolve.root = TRUE, edgelabel = TRUE)
  tree <- ladderize(tree)
  cat("Rooted with", length(onychorhynchus_tips), "Onychorhynchus tip(s) as outgroup\n")
  cat("Tree ladderized for better visualization\n\n")
} else {
  cat("Warning: No Onychorhynchus tips found, tree not rooted\n\n")
}

# Extract genus from species names for coloring
tip_data <- data.frame(
  label = tree$tip.label,
  genus = sapply(strsplit(tree$tip.label, "_"), `[`, 1),
  is_tityra = grepl("^Tityra", tree$tip.label),
  is_outgroup = grepl("^Onychorhynchus", tree$tip.label),
  stringsAsFactors = FALSE
)

# Create rectangular tree plot
cat("Creating rectangular tree plot...\n")
p_rect <- ggtree(tree, layout = "roundrect", ladderize = TRUE, right = TRUE) %<+% tip_data +
  geom_tiplab(aes(color = genus, 
                  fontface = ifelse(is_tityra, "bold", "plain"),
                  size = ifelse(is_tityra, 4, 3)), 
              align = FALSE) +
  geom_nodepoint(color = "darkblue", size = 2, shape = 21, fill = "lightblue") +
  geom_nodelab(aes(label = label),
               size = 2.5,
               hjust = 1.2,
               vjust = -0.3,
               color = "blue") +
  geom_tippoint(aes(subset = is_tityra), 
                color = "red", size = 3, shape = 17) +
  geom_tippoint(aes(subset = is_outgroup), 
                color = "darkgreen", size = 3, shape = 15) +
  theme_tree2() +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold")) +
  ggtitle("Concatenated Mitochondrial Genes Phylogeny - Rooted with Onychorhynchus") +
  scale_color_discrete(name = "Genus") +
  scale_size_identity()

# Save rectangular plot
output_file <- file.path(plot_dir, "concatenated_tree_rectangular.pdf")
ggsave(output_file, p_rect, width = 14, height = 12)
cat("Saved rectangular tree:", output_file, "\n\n")

# Create circular tree plot
cat("Creating circular tree plot...\n")
p_circ <- ggtree(tree, layout = "circular") %<+% tip_data +
  geom_tiplab(aes(color = genus,
                  fontface = ifelse(is_tityra, "bold", "plain"),
                  size = ifelse(is_tityra, 3.5, 2.5)), 
              align = FALSE) +
  geom_nodepoint(color = "darkblue", size = 2, shape = 21, fill = "lightblue") +
  geom_nodelab(aes(label = label),
               size = 2,
               hjust = 0.5,
               vjust = -0.5,
               color = "blue") +
  geom_tippoint(aes(subset = is_tityra), 
                color = "red", size = 3, shape = 17) +
  geom_tippoint(aes(subset = is_outgroup), 
                color = "darkgreen", size = 3, shape = 15) +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold")) +
  ggtitle("Concatenated Mitochondrial Genes Phylogeny - Circular (Rooted with Onychorhynchus)") +
  scale_color_discrete(name = "Genus") +
  scale_size_identity()

# Save circular plot
output_file <- file.path(plot_dir, "concatenated_tree_circular.pdf")
ggsave(output_file, p_circ, width = 14, height = 14)
cat("Saved circular tree:", output_file, "\n\n")

# Create a plot with bootstrap support highlighted
cat("Creating bootstrap support plot...\n")

# Parse bootstrap values
tree_with_support <- tree
bootstrap_values <- as.numeric(tree_with_support$node.label)

# Get number of tips
ntips <- length(tree$tip.label)

# Create data frame with support categories
node_data <- data.frame(
  node = (ntips + 1):(ntips + tree$Nnode),
  support = bootstrap_values,
  support_category = cut(bootstrap_values, 
                         breaks = c(0, 50, 70, 90, 100),
                         labels = c("< 50", "50-70", "70-90", "90-100"),
                         include.lowest = TRUE)
)

p_support <- ggtree(tree, layout = "roundrect", ladderize = TRUE, right = TRUE) %<+% tip_data +
  geom_tiplab(aes(color = genus,
                  fontface = ifelse(is_tityra, "bold", "plain"),
                  size = ifelse(is_tityra, 4, 3)),
              align = FALSE) +
  geom_nodepoint(color = "darkblue", size = 3, shape = 21, fill = "lightblue") +
  geom_nodelab(aes(label = label),
               size = 2.5,
               hjust = 1.2,
               vjust = -0.3,
               color = "blue") +
  geom_tippoint(aes(subset = is_tityra),
                color = "red", size = 3, shape = 17) +
  geom_tippoint(aes(subset = is_outgroup), 
                color = "darkgreen", size = 3, shape = 15) +
  theme_tree2() +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold")) +
  ggtitle("Concatenated Genes Phylogeny with Bootstrap Support (Rooted with Onychorhynchus)") +
  scale_color_discrete(name = "Genus") +
  scale_size_identity()

# Save support plot
output_file <- file.path(plot_dir, "concatenated_tree_support.pdf")
ggsave(output_file, p_support, width = 14, height = 12)
cat("Saved bootstrap support tree:", output_file, "\n\n")

cat("Tree plotting completed!\n")

R_SCRIPT

echo ""
echo "==================================================================="
echo "Analysis completed successfully!"
echo "Completed: $(date)"
echo "==================================================================="
echo ""
echo "Output files:"
echo "  Tree: ${TREE_DIR}/concatenated_all_genes.treefile"
echo "  Plots: ${PLOT_DIR}/concatenated_tree_*.pdf"
echo ""
