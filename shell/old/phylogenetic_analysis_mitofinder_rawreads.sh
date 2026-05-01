#!/bin/bash
###############################################################################
# Phylogenetic Analysis of Combined Mitochondrial Genes
# - Performs multiple sequence alignment with MAFFT
# - Builds phylogenetic trees with IQ-TREE
# - Plots trees with R/ggplot2
###############################################################################

# Set working directory
WD=${WD:-/media/inter/mkapun/projects/Tityra}

# Input/Output directories
INPUT_DIR="${WD}/results/combined_mito_genes"
OUTPUT_DIR="${WD}/results/phylogenetic_analysis"
ALIGNMENT_DIR="${OUTPUT_DIR}/alignments"
TREE_DIR="${OUTPUT_DIR}/trees"
PLOT_DIR="${OUTPUT_DIR}/plots"

# Create output directories
mkdir -p ${ALIGNMENT_DIR}
mkdir -p ${TREE_DIR}
mkdir -p ${PLOT_DIR}

echo "==================================================================="
echo "Phylogenetic Analysis Pipeline"
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

# List of genes to process
GENES=(ATP6 ATP8 CO2 CO3 COI CYTB ND1 ND2 ND3 ND4 ND4L ND5)

for GENE in "${GENES[@]}"; do
    INPUT_FILE="${INPUT_DIR}/${GENE}_combined.fasta"
    OUTPUT_FILE="${ALIGNMENT_DIR}/${GENE}_aligned.fasta"
    
    if [ -f "${INPUT_FILE}" ]; then
        echo ""
        echo "-------------------------------------------------------------------"
        echo "Aligning ${GENE}..."
        echo "-------------------------------------------------------------------"
        
        # Run MAFFT with auto mode and reverse complement detection
        # --adjustdirection: automatically reverse complement sequences if needed
        mafft --auto \
              --adjustdirection \
              --thread 4 \
              --quiet \
              "${INPUT_FILE}" > "${OUTPUT_FILE}.tmp"
        
        if [ $? -eq 0 ]; then
            # Remove _R_ prefix added by MAFFT for reverse-complemented sequences
            sed 's/>_R_/>/g' "${OUTPUT_FILE}.tmp" > "${OUTPUT_FILE}"
            rm "${OUTPUT_FILE}.tmp"
            echo "  ✓ ${GENE} alignment completed"
            echo "  Output: ${OUTPUT_FILE}"
        else
            echo "  ✗ Error: ${GENE} alignment failed"
        fi
    else
        echo "  ⚠ Warning: Input file not found: ${INPUT_FILE}"
    fi
done

conda deactivate
echo ""
echo "Multiple sequence alignment completed"
echo ""

###############################################################################
# 2. Remove Gappy Regions with Gblocks
###############################################################################
echo "==================================================================="
echo "Step 2: Removing Gappy Regions with Gblocks"
echo "==================================================================="

# Load Gblocks module
module load Alignment/Gblocks-0.91b

for GENE in "${GENES[@]}"; do
    ALIGNMENT_FILE="${ALIGNMENT_DIR}/${GENE}_aligned.fasta"
    GBLOCKS_OUTPUT="${ALIGNMENT_DIR}/${GENE}_aligned_gblocks.fasta"
    
    if [ -f "${ALIGNMENT_FILE}" ]; then
        echo ""
        echo "-------------------------------------------------------------------"
        echo "Processing ${GENE} with Gblocks..."
        echo "-------------------------------------------------------------------"
        
        # Run Gblocks
        # -t=d: DNA sequences
        # -b1=0.5: minimum sequences for conserved position (50%)
        # -b2=0.85: minimum sequences for flanking position (85%)
        # -b3=8: maximum contiguous nonconserved positions
        # -b4=10: minimum block length
        # -b5=h: allow gap positions within final blocks (half)
        Gblocks "${ALIGNMENT_FILE}" -t=d -b1=0.5 -b2=0.85 -b3=8 -b4=10 -b5=h
        
        # Gblocks creates output with -gb extension
        GB_OUTPUT="${ALIGNMENT_FILE}-gb"
        
        if [ -f "${GB_OUTPUT}" ]; then
            # Convert Gblocks output back to FASTA format
            mv "${GB_OUTPUT}" "${GBLOCKS_OUTPUT}"
            echo "  ✓ ${GENE} Gblocks processing completed"
            echo "  Output: ${GBLOCKS_OUTPUT}"
        else
            echo "  ⚠ Warning: Gblocks output not found, using original alignment"
            cp "${ALIGNMENT_FILE}" "${GBLOCKS_OUTPUT}"
        fi
    else
        echo "  ⚠ Warning: Alignment file not found: ${ALIGNMENT_FILE}"
    fi
done

echo ""
echo "Gblocks processing completed"
echo ""

###############################################################################
# 3. Phylogenetic Tree Reconstruction with IQ-TREE
###############################################################################
echo "==================================================================="
echo "Step 3: Phylogenetic Tree Reconstruction with IQ-TREE"
echo "==================================================================="

# Activate IQ-TREE environment (using version 2.3.3)
conda activate iq-tree-2.3.3

for GENE in "${GENES[@]}"; do
    # Use Gblocks-processed alignment
    ALIGNMENT_FILE="${ALIGNMENT_DIR}/${GENE}_aligned_gblocks.fasta"
    
    if [ -f "${ALIGNMENT_FILE}" ]; then
        echo ""
        echo "-------------------------------------------------------------------"
        echo "Building tree for ${GENE}..."
        echo "-------------------------------------------------------------------"
        
        # Change to tree directory to keep IQ-TREE output organized
        cd ${TREE_DIR}
        
        # Run IQ-TREE with:
        # -s: alignment file
        # -m MFP: automatically select best-fit model
        # -bb 1000: ultrafast bootstrap with 1000 replicates
        # -nt 4: use 4 threads
        # -pre: output prefix
        iqtree -s "${ALIGNMENT_FILE}" \
               -m MFP \
               -bb 1000 \
               -nt 4 \
               -pre "${GENE}" \
               -quiet
        
        if [ $? -eq 0 ]; then
            echo "  ✓ ${GENE} tree reconstruction completed"
            echo "  Output tree: ${TREE_DIR}/${GENE}.treefile"
        else
            echo "  ✗ Error: ${GENE} tree reconstruction failed"
        fi
        
        cd ${WD}
    else
        echo "  ⚠ Warning: Alignment file not found: ${ALIGNMENT_FILE}"
    fi
done

conda deactivate
echo ""
echo "Phylogenetic tree reconstruction completed"
echo ""

###############################################################################
# 4. Plot Trees with R/ggplot2
###############################################################################
echo "==================================================================="
echo "Step 4: Plotting Phylogenetic Trees with R/ggplot2"
echo "==================================================================="

# Create R script for tree plotting
R_SCRIPT="${OUTPUT_DIR}/plot_trees.R"

cat > ${R_SCRIPT} << 'RSCRIPT_EOF'
#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
})

# Get directories from environment
tree_dir <- Sys.getenv("TREE_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")

# List of genes
genes <- c("ATP6", "ATP8", "CO2", "CO3", "COI", "CYTB", 
           "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5")

cat("Plotting phylogenetic trees...\n\n")

for (gene in genes) {
  tree_file <- file.path(tree_dir, paste0(gene, ".treefile"))
  
  if (file.exists(tree_file)) {
    cat("Plotting", gene, "tree...\n")
    
    # Read tree
    tree <- read.tree(tree_file)
    
    # Root tree with Onychorhynchus as outgroup
    onychorhynchus_tips <- grep("^Onychorhynchus", tree$tip.label, value = TRUE)
    
    if (length(onychorhynchus_tips) > 0) {
      # Root with Onychorhynchus species
      tree <- root(tree, outgroup = onychorhynchus_tips, resolve.root = TRUE)
      cat("  Rooted with", length(onychorhynchus_tips), "Onychorhynchus tip(s)\n")
    } else {
      cat("  Warning: No Onychorhynchus tips found, tree not rooted\n")
    }
    
    # Extract genus from species names for coloring
    tip_data <- data.frame(
      label = tree$tip.label,
      genus = sapply(strsplit(tree$tip.label, "_"), `[`, 1),
      is_target = grepl("^Tityra_leucura", tree$tip.label),
      stringsAsFactors = FALSE
    )
    
    # Create plot with Tityra_leucura highlighted
    p <- ggtree(tree, layout = "rectangular") %<+% tip_data +
      geom_tiplab(aes(color = genus, fontface = ifelse(is_target, "bold", "plain"),
                      size = ifelse(is_target, 4, 3)), align = FALSE) +
      geom_nodelab(aes(label = label),
                   size = 2.5,
                   hjust = 1.2,
                   vjust = -0.3,
                   color = "blue") +
      geom_tippoint(aes(subset = is_target), color = "red", size = 3, shape = 17) +
      theme_tree2() +
      theme(legend.position = "right",
            plot.title = element_text(size = 14, face = "bold")) +
      ggtitle(paste("Phylogenetic Tree:", gene)) +
      scale_color_discrete(name = "Genus") +
      scale_size_identity()
    
    # Save plot
    output_file <- file.path(plot_dir, paste0(gene, "_tree.pdf"))
    ggsave(output_file, p, width = 12, height = 10)
    cat("  Saved:", output_file, "\n")
    
    # Also create a circular layout
    p_circular <- ggtree(tree, layout = "circular") %<+% tip_data +
      geom_tiplab(aes(color = genus, fontface = ifelse(is_target, "bold", "plain"),
                      size = ifelse(is_target, 3.5, 2.5)), align = FALSE) +
      geom_nodelab(aes(label = label),
                   size = 2,
                   hjust = 0.5,
                   vjust = -0.5,
                   color = "blue") +
      geom_tippoint(aes(subset = is_target), color = "red", size = 3, shape = 17) +
      theme(legend.position = "right") +
      ggtitle(paste("Phylogenetic Tree (Circular):", gene)) +
      scale_color_discrete(name = "Genus") +
      scale_size_identity()
    
    output_file_circular <- file.path(plot_dir, paste0(gene, "_tree_circular.pdf"))
    ggsave(output_file_circular, p_circular, width = 12, height = 12)
    cat("  Saved:", output_file_circular, "\n\n")
    
  } else {
    cat("  Warning: Tree file not found:", tree_file, "\n\n")
  }
}

cat("Tree plotting completed!\n")
RSCRIPT_EOF

# Make R script executable
chmod +x ${R_SCRIPT}

# Run R script
echo "Running R script to plot trees..."
export TREE_DIR
export PLOT_DIR

Rscript ${R_SCRIPT}

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Tree plotting completed"
    echo "  Output directory: ${PLOT_DIR}"
else
    echo ""
    echo "✗ Error: Tree plotting failed"
fi

###############################################################################
# 5. Summary
###############################################################################
echo ""
echo "==================================================================="
echo "Phylogenetic Analysis Pipeline Completed: $(date)"
echo "==================================================================="
echo ""
echo "Summary:"
echo "  Alignments: ${ALIGNMENT_DIR}"
echo "  Trees: ${TREE_DIR}"
echo "  Plots: ${PLOT_DIR}"
echo ""
echo "Generated files:"
echo "  - Individual gene alignments (12 genes)"
echo "  - Gblocks-processed alignments (12 genes)"
echo "  - Individual gene trees (12 genes)"
echo "  - Tree plots (rectangular and circular layouts)"
echo ""
echo "==================================================================="
