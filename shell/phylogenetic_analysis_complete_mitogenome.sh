#!/bin/bash
###############################################################################
# Phylogenetic Analysis of Complete Mitochondrial Genomes
# Author: Generated for Tityra project
# Description: Build phylogenetic tree from complete mitogenome alignment
###############################################################################

WD="/media/inter/mkapun/projects/Tityra"

# Create output directories
mkdir -p "${WD}/results/complete_mitogenome_phylogeny/alignment"
mkdir -p "${WD}/results/complete_mitogenome_phylogeny/trees"

# Environment setup - create phylogenetics environment if it doesn't exist
CONDA_ENV="${WD}/scripts/programs/phylogenetics"
if [ ! -d "${CONDA_ENV}" ]; then
    echo "Creating phylogenetics conda environment..."
    mamba create -y -p "${CONDA_ENV}" -c bioconda -c conda-forge \
        iqtree=2.3.3 mafft raxml-ng fasttree python=3.10
    echo "Phylogenetics environment created!"
fi

# Activate environment with IQ-TREE
# Use eval to properly activate conda in bash script
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV}"

# Check which alignment files exist and run analysis on both if available
alignment_annotated="${WD}/results/complete_mitogenome_alignment/annotation_aware_aligned.fasta"
alignment_whole="${WD}/results/complete_mitogenome_alignment/whole_genome_aligned.fasta"

ALIGN_FILES=()
ALIGN_NAMES=()

if [ -f "${alignment_annotated}" ]; then
    ALIGN_FILES+=("${alignment_annotated}")
    ALIGN_NAMES+=("annotation_aware")
    echo "Found: annotation-aware alignment"
fi

if [ -f "${alignment_whole}" ]; then
    ALIGN_FILES+=("${alignment_whole}")
    ALIGN_NAMES+=("whole_genome")
    echo "Found: whole-genome alignment"
fi

if [ ${#ALIGN_FILES[@]} -eq 0 ]; then
    echo "Error: No alignment files found"
    echo "Expected: ${alignment_annotated} or ${alignment_whole}"
    exit 1
fi

echo ""
echo "==========================================="
echo "Phylogenetic Analysis of Complete Mitogenomes"
echo "Will analyze ${#ALIGN_FILES[@]} alignment(s)"
echo "==========================================="
echo ""

for i in "${!ALIGN_FILES[@]}"; do
    alignment_file="${ALIGN_FILES[$i]}"
    alignment_name="${ALIGN_NAMES[$i]}"
    
    echo ""
    echo "-------------------------------------------"
    echo "Processing: ${alignment_name}"
    echo "-------------------------------------------"
    
    # Create subdirectory for this alignment type
    mkdir -p "${WD}/results/complete_mitogenome_phylogeny/alignment/${alignment_name}"
    mkdir -p "${WD}/results/complete_mitogenome_phylogeny/trees/${alignment_name}"
    
    # Copy alignment
    cp "${alignment_file}" "${WD}/results/complete_mitogenome_phylogeny/alignment/${alignment_name}/"
    
    cd "${WD}/results/complete_mitogenome_phylogeny/trees/${alignment_name}"
    
    alignment_input="${WD}/results/complete_mitogenome_phylogeny/alignment/${alignment_name}/$(basename ${alignment_file})"

    # Check if tree already exists
    tree_file="complete_mitogenome_${alignment_name}.treefile"
    
    if [ -f "${tree_file}" ]; then
        echo ""
        echo "Tree file already exists for ${alignment_name}"
        echo "Skipping IQ-TREE analysis to preserve existing results"
        echo "Delete ${tree_file} if you want to rerun IQ-TREE"
        echo ""
    else
        echo ""
        echo "Step 1: Running IQ-TREE with automatic model selection..."
        echo "-------------------------------------------"

        # Check if this is annotation-aware alignment and partitions file exists
        partitions_file="${WD}/results/complete_mitogenome_alignment/partitions.txt"
        
        if [[ "${alignment_name}" == "annotation_aware" ]] && [[ -f "${partitions_file}" ]]; then
            echo "Using partitioned model (separate models per gene)"
            echo "Partitions file: ${partitions_file}"
            
            # Copy partitions file to working directory
            cp "${partitions_file}" ./
            
            # Remove old checkpoint files if they exist to avoid conflicts
            rm -f complete_mitogenome_${alignment_name}.model.gz
            rm -f complete_mitogenome_${alignment_name}.ckp.gz
            
            # Run IQ-TREE with partitions:
            # -s: input alignment
            # -spp: partition file with separate models per partition
            # -m: automatic model selection (ModelFinder) for each partition
            # -bb: ultrafast bootstrap with 1000 replicates
            # -nt: number of threads
            # -pre: prefix for output files
            
        iqtree2 -s "${alignment_input}" \
            -spp partitions.txt \
            -m MFP \
            -bb 1000 \
            -nt AUTO \
            -pre complete_mitogenome_${alignment_name}
        else
            if [[ "${alignment_name}" == "annotation_aware" ]]; then
                echo "Warning: Partitions file not found at ${partitions_file}"
                echo "Running with single model"
            else
                echo "Using single model (whole-genome alignment)"
            fi
            
            # Remove old checkpoint files if they exist to avoid conflicts
            rm -f complete_mitogenome_${alignment_name}.model.gz
            rm -f complete_mitogenome_${alignment_name}.ckp.gz
            
            # Run IQ-TREE with:
            # -s: input alignment
            # -m: automatic model selection (ModelFinder)
            # -bb: ultrafast bootstrap with 1000 replicates
            # -nt: number of threads
            # -pre: prefix for output files
            
            iqtree2 -s "${alignment_input}" \
                -m MFP \
                -bb 1000 \
                -nt AUTO \
                -pre complete_mitogenome_${alignment_name}
        fi

        echo ""
        echo "Step 2: Generating tree summary..."
        echo "-------------------------------------------"

        # Check if tree was generated
        if [ -f "complete_mitogenome_${alignment_name}.treefile" ]; then
            echo "Tree successfully generated!"
            
            # Copy important files with descriptive names
            cp complete_mitogenome_${alignment_name}.treefile complete_mitogenome_${alignment_name}_ML.tree
            cp complete_mitogenome_${alignment_name}.contree complete_mitogenome_${alignment_name}_consensus.tree 2>/dev/null || true
            
            # Print tree summary
            echo ""
            echo "=== Tree Summary ==="
            cat complete_mitogenome_${alignment_name}.iqtree | grep -A 20 "MAXIMUM LIKELIHOOD TREE" || true
            
            echo ""
            echo "=== Best Model ==="
            cat complete_mitogenome_${alignment_name}.iqtree | grep -A 5 "Best-fit model" || true
        else
            echo "Error: Tree file not generated for ${alignment_name}!"
        fi
    fi
    
    echo ""
    echo "Step 3: Plotting phylogenetic tree with ggtree..."
    echo "-------------------------------------------"
    
    if [ -f "complete_mitogenome_${alignment_name}.treefile" ]; then
        # Create plots directory
        mkdir -p "${WD}/results/complete_mitogenome_phylogeny/plots/${alignment_name}"
        
        # Export variables for R
        export TREE_FILE="${WD}/results/complete_mitogenome_phylogeny/trees/${alignment_name}/complete_mitogenome_${alignment_name}.treefile"
        export PLOT_DIR="${WD}/results/complete_mitogenome_phylogeny/plots/${alignment_name}"
        export ALIGNMENT_TYPE="${alignment_name}"
        
        Rscript - << 'R_SCRIPT'

# Load required libraries
suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
})

# Get parameters from environment
tree_file <- Sys.getenv("TREE_FILE")
plot_dir <- Sys.getenv("PLOT_DIR")
alignment_type <- Sys.getenv("ALIGNMENT_TYPE")

if (!file.exists(tree_file)) {
  stop("Tree file not found: ", tree_file)
}

cat("Reading tree file:", tree_file, "\n")

# Read tree
tree <- read.tree(tree_file)

cat("Tree successfully loaded\n")
cat("Number of tips:", length(tree$tip.label), "\n")
cat("Tips:", paste(tree$tip.label, collapse=", "), "\n\n")

# Replace underscores with spaces for nicer labels
tree$tip.label <- gsub("_", " ", tree$tip.label)

# Root tree with Pitta sordida as outgroup (most distant)
pitta_tips <- grep("^Pitta sordida$", tree$tip.label, value = TRUE)

if (length(pitta_tips) > 0) {
  tree <- root(tree, outgroup = pitta_tips, resolve.root = TRUE, edgelabel = TRUE)
  tree <- ladderize(tree)
  cat("Rooted with Pitta sordida as outgroup\n")
  cat("Tree ladderized for better visualization\n\n")
} else {
  cat("Warning: Pitta sordida not found, tree not rooted\n\n")
}

# Extract genus from species names for coloring
tip_data <- data.frame(
  label = tree$tip.label,
  genus = sapply(strsplit(tree$tip.label, " "), `[`, 1),
  is_tityra = grepl("^Tityra", tree$tip.label),
  is_outgroup = grepl("^Pitta", tree$tip.label),
  stringsAsFactors = FALSE
)

# Create title based on alignment type
if (alignment_type == "annotation_aware") {
  title_suffix <- "Annotation-Aware (Protein-Coding Genes)"
} else {
  title_suffix <- "Whole Genome"
}

# Create rectangular tree plot
cat("Creating rectangular tree plot...\n")
p_rect <- ggtree(tree, layout = "roundrect", ladderize = TRUE, right = TRUE) %<+% tip_data +
  geom_tiplab(aes(color = genus, 
                  fontface = ifelse(is_tityra, "bold.italic", "italic"),
                  size = ifelse(is_tityra, 4, 3.5)), 
              align = FALSE, offset = 0.0001) +
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
  ggtitle(paste0("Complete Mitogenome Phylogeny - ", title_suffix, " (Rooted with Pitta sordida)")) +
  scale_color_discrete(name = "Genus") +
  scale_size_identity()

# Save rectangular plot
output_file <- file.path(plot_dir, paste0("complete_mitogenome_", alignment_type, "_rectangular.pdf"))
ggsave(output_file, p_rect, width = 14, height = 10)
cat("Saved rectangular tree:", output_file, "\n\n")

# Create circular tree plot
cat("Creating circular tree plot...\n")
p_circ <- ggtree(tree, layout = "circular") %<+% tip_data +
  geom_tiplab(aes(color = genus,
                  fontface = ifelse(is_tityra, "bold.italic", "italic"),
                  size = ifelse(is_tityra, 3.5, 3)), 
              align = FALSE, offset = 0.0005) +
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
  ggtitle(paste0("Complete Mitogenome Phylogeny - ", title_suffix, " (Circular, Rooted with Pitta sordida)")) +
  scale_color_discrete(name = "Genus") +
  scale_size_identity()

# Save circular plot
output_file <- file.path(plot_dir, paste0("complete_mitogenome_", alignment_type, "_circular.pdf"))
ggsave(output_file, p_circ, width = 12, height = 12)
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
                  fontface = ifelse(is_tityra, "bold.italic", "italic"),
                  size = ifelse(is_tityra, 4, 3.5)),
              align = FALSE, offset = 0.0001) +
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
  ggtitle(paste0("Complete Mitogenome Phylogeny - ", title_suffix, " with Bootstrap Support")) +
  scale_color_discrete(name = "Genus") +
  scale_size_identity()

# Save support plot
output_file <- file.path(plot_dir, paste0("complete_mitogenome_", alignment_type, "_support.pdf"))
ggsave(output_file, p_support, width = 14, height = 10)
cat("Saved bootstrap support tree:", output_file, "\n\n")

cat("Tree plotting completed!\n")

R_SCRIPT

        echo "Tree plots saved to: ${PLOT_DIR}"
    else
        echo "Skipping plots - tree file not found"
    fi
    
    echo ""
    echo "Analysis for ${alignment_name} complete!"

done

echo ""
echo "==========================================="
echo "Phylogenetic analysis complete!"
echo "==========================================="
echo "Results in: ${WD}/results/complete_mitogenome_phylogeny/trees/"
echo ""
echo "Output directories:"
for alignment_name in "${ALIGN_NAMES[@]}"; do
    echo "  - ${alignment_name}/"
done
echo ""

conda deactivate
