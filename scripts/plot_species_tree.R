#!/usr/bin/env Rscript

# Plot ASTRAL species tree from BUSCO phylogeny analysis

library(ape)

# Read the tree
tree_file <- "/media/inter/mkapun/projects/Tityra/results/BUSCO_aves_phylogeny/astral/species_tree.tre"
tree <- read.tree(tree_file)

# Output file
output_pdf <- "/media/inter/mkapun/projects/Tityra/results/BUSCO_aves_phylogeny/astral/species_tree.pdf"

# Create PDF
pdf(output_pdf, width = 10, height = 8)

# Plot tree with support values
plot(tree,
     main = "ASTRAL Species Tree (303 BUSCO genes, 14 bird species)",
     cex = 1.2,
     font = 3,
     edge.width = 2,
     label.offset = 0.01
)

# Add support values to nodes
nodelabels(tree$node.label,
     frame = "circle",
     cex = 0.8,
     bg = "lightblue"
)

# Add axis
axisPhylo(cex.axis = 0.8)

dev.off()

cat("Tree plot saved to:", output_pdf, "\n")

# Also create a simple text summary
cat("\n=== Tree Summary ===\n")
cat("Number of tips:", Ntip(tree), "\n")
cat("Number of nodes:", Nnode(tree), "\n")
cat("Tree is rooted:", is.rooted(tree), "\n")
cat("\nSpecies:\n")
print(tree$tip.label)
