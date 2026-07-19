## load R parameters from Commandline
args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
plot_dir <- args[2]
outgroup <- args[3]
fileprefix <- args[4]
# Optional 5th argument: tab-separated file mapping old tip label -> new label
label_map_file <- ifelse(length(args) >= 5 && file.exists(args[5]), args[5], NULL)
# Manacus_manacus_KJ810477.1

# Load required libraries
suppressPackageStartupMessages({
    library(ape)
    library(ggtree)
    library(ggplot2)
    library(dplyr)
})

# # Get parameters from environment
# tree_file <- Sys.getenv("TREE_FILE")
# plot_dir <- Sys.getenv("PLOT_DIR")
# alignment_type <- Sys.getenv("ALIGNMENT_TYPE")

if (!file.exists(tree_file)) {
    stop("Tree file not found: ", tree_file)
}

cat("Reading tree file:", tree_file, "\n")

# Read tree
tree <- read.tree(tree_file)

cat("Tree successfully loaded\n")
cat("Number of tips:", length(tree$tip.label), "\n")
cat("Tips:", paste(tree$tip.label, collapse = ", "), "\n\n")

# Replace underscores with spaces for nicer labels
tree$tip.label <- gsub("_", " ", tree$tip.label)

# Root tree with specified outgroup FIRST (before relabelling, so the
# outgroup pattern still matches the original, space-converted tip labels).
outgroup_tips <- grep(outgroup, tree$tip.label, value = TRUE)

if (length(outgroup_tips) > 0) {
    tree <- root(tree, outgroup = outgroup_tips, resolve.root = TRUE, edgelabel = TRUE)
    tree <- ladderize(tree)
    cat("Rooted with", outgroup, "as outgroup\n")
    cat("Tree ladderized for better visualization\n\n")
} else {
    cat("Warning:", outgroup, "not found, tree not rooted\n\n")
}

# Apply tip-label relabelling from mapping file (old -> new). The mapping
# file uses the ORIGINAL (underscored) tip labels on the left, which we
# convert to space form to match the current tree tip labels.
if (!is.null(label_map_file)) {
    cat("Applying tip-label mapping from:", label_map_file, "\n")
    map <- read.delim(label_map_file, header = FALSE,
                      col.names = c("old", "new"),
                      colClasses = "character", stringsAsFactors = FALSE)
    # Convert keys (stored with underscores) to space form to match tree
    map$old_space <- gsub("_", " ", map$old)
    lab_lookup <- setNames(map$new, map$old_space)
    idx <- match(tree$tip.label, names(lab_lookup))
    replaced <- sum(!is.na(idx))
    tree$tip.label[!is.na(idx)] <- lab_lookup[na.omit(idx)]
    cat("Relabelled", replaced, "tips\n\n")
}

# Extract genus from species names for coloring
tip_data <- data.frame(
    label = tree$tip.label,
    genus = sapply(strsplit(tree$tip.label, " "), `[`, 1),
    is_tityra = grepl("^Tityra", tree$tip.label),
    is_outgroup = grepl(paste0("^", outgroup), tree$tip.label),
    stringsAsFactors = FALSE
)

# Create rectangular tree plot
cat("Creating rectangular tree plot...\n")
p_rect <- ggtree(tree, layout = "roundrect", ladderize = TRUE, right = TRUE) %<+% tip_data +
    geom_tiplab(
        aes(
            color = genus,
            fontface = ifelse(is_tityra, "bold.italic", "italic"),
            size = ifelse(is_tityra, 4, 3.5)
        ),
        align = FALSE, offset = 0.0001
    ) +
    geom_nodepoint(color = "darkblue", size = 2, shape = 21, fill = "lightblue") +
    geom_nodelab(aes(label = label),
        size = 2.5,
        hjust = 1.2,
        vjust = -0.3,
        color = "blue"
    ) +
    geom_tippoint(aes(subset = is_tityra),
        color = "red", size = 3, shape = 17
    ) +
    geom_tippoint(aes(subset = is_outgroup),
        color = "darkgreen", size = 3, shape = 15
    ) +
    theme_tree2() +
    theme(
        legend.position = "right",
        plot.title = element_text(size = 14, face = "bold")
    ) +
    ggtitle(paste0("Phylogeny (Rooted with ", outgroup, ")")) +
    scale_color_discrete(name = "Genus") +
    scale_size_identity()

# Save rectangular plot
output_file <- file.path(plot_dir, paste0(fileprefix, "_rectangular.pdf"))
ggsave(output_file, p_rect, width = 14, height = 10)
cat("Saved rectangular tree:", output_file, "\n\n")

# Create circular tree plot
cat("Creating circular tree plot...\n")
p_circ <- ggtree(tree, layout = "circular") %<+% tip_data +
    geom_tiplab(
        aes(
            color = genus,
            fontface = ifelse(is_tityra, "bold.italic", "italic"),
            size = ifelse(is_tityra, 3.5, 3)
        ),
        align = FALSE, offset = 0.0005
    ) +
    geom_nodepoint(color = "darkblue", size = 2, shape = 21, fill = "lightblue") +
    geom_nodelab(aes(label = label),
        size = 2,
        hjust = 0.5,
        vjust = -0.5,
        color = "blue"
    ) +
    geom_tippoint(aes(subset = is_tityra),
        color = "red", size = 3, shape = 17
    ) +
    geom_tippoint(aes(subset = is_outgroup),
        color = "darkgreen", size = 3, shape = 15
    ) +
    theme(
        legend.position = "right",
        plot.title = element_text(size = 14, face = "bold")
    ) +
    ggtitle(paste0("Phylogeny (Circular, Rooted with ", outgroup, ")")) +
    scale_color_discrete(name = "Genus") +
    scale_size_identity()

# Save circular plot
output_file <- file.path(plot_dir, paste0(fileprefix, "_circular.pdf"))
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
        include.lowest = TRUE
    )
)

p_support <- ggtree(tree, layout = "roundrect", ladderize = TRUE, right = TRUE) %<+% tip_data +
    geom_tiplab(
        aes(
            color = genus,
            fontface = ifelse(is_tityra, "bold.italic", "italic"),
            size = ifelse(is_tityra, 4, 3.5)
        ),
        align = FALSE, offset = 0.0001
    ) +
    geom_nodepoint(color = "darkblue", size = 3, shape = 21, fill = "lightblue") +
    geom_nodelab(aes(label = label),
        size = 2.5,
        hjust = 1.2,
        vjust = -0.3,
        color = "blue"
    ) +
    geom_tippoint(aes(subset = is_tityra),
        color = "red", size = 3, shape = 17
    ) +
    geom_tippoint(aes(subset = is_outgroup),
        color = "darkgreen", size = 3, shape = 15
    ) +
    theme_tree2() +
    theme(
        legend.position = "right",
        plot.title = element_text(size = 14, face = "bold")
    ) +
    ggtitle(paste0("Phylogeny with Bootstrap Support")) +
    scale_color_discrete(name = "Genus") +
    scale_size_identity()

# Save support plot
output_file <- file.path(plot_dir, paste0(fileprefix, "_support.pdf"))
ggsave(output_file, p_support, width = 14, height = 10)
cat("Saved bootstrap support tree:", output_file, "\n\n")

cat("Tree plotting completed!\n")
