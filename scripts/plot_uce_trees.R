#!/usr/bin/env Rscript
# Plot UCE phylogeny trees (both full and trimmed analyses)
# Usage: Rscript plot_uce_trees.R <prefix> <plot_dir> <outgroup_pattern>

suppressPackageStartupMessages({
    library(ape)
    library(ggtree)
    library(ggplot2)
    library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript plot_uce_trees.R <prefix> <plot_dir> [outgroup_pattern]")
}

prefix          <- args[1]
plot_dir        <- args[2]
outgroup_pat    <- ifelse(length(args) >= 3, args[3], "Tr_musculus")

tree_file <- paste0(prefix, ".treefile")
if (!file.exists(tree_file)) {
    stop("Tree file not found: ", tree_file)
}

tree <- read.tree(tree_file)

# Root with outgroup
outgroup <- grep(outgroup_pat, tree$tip.label, value = TRUE)
if (length(outgroup) > 0) {
    tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
}
tree <- ladderize(tree)

# Load metadata and prepare subspecies information
# Try to find Meta.csv relative to tree file location or project root
tree_dir <- dirname(normalizePath(tree_file))
possible_meta_paths <- c(
    file.path(tree_dir, "Meta.csv"),
    file.path(dirname(tree_dir), "..", "data", "openwings_tityrinae", "Meta.csv"),
    file.path(getwd(), "data", "openwings_tityrinae", "Meta.csv")
)

meta_file <- ""
for (path in possible_meta_paths) {
    if (file.exists(path)) {
        meta_file <- path
        break
    }
}

meta_data <- data.frame()
if (file.exists(meta_file)) {
    cat("Loading metadata from:", meta_file, "\n")
    meta_data <- read.csv(meta_file, stringsAsFactors = FALSE)
    cat("Loaded", nrow(meta_data), "metadata records\n")
    # Ensure character columns for comparison
    meta_data$Genus <- as.character(meta_data$Genus)
    meta_data$Species <- as.character(meta_data$Species)
    meta_data$Subspecies <- as.character(meta_data$Subspecies)
} else {
    cat("WARNING: Meta.csv not found at:", meta_file, "\n")
}

# Function to create display name with subspecies
create_display_name <- function(tip_label, meta) {
    # Try to match directly with Coded_as column in Meta.csv
    if (nrow(meta) > 0 && "Coded_as" %in% names(meta)) {
        match_row <- meta[meta$Coded_as == tip_label, ]
        if (nrow(match_row) > 0) {
            genus <- match_row$Genus[1]
            species <- match_row$Species[1]
            subspecies <- match_row$Subspecies[1]
            
            # Abbreviate genus to first letter + period (e.g., T. for Tityra)
            #genus <- paste0(substr(genus, 1, 1), ".")
            
            # Build display name with subspecies if available
            if (!is.na(subspecies) && subspecies != "") {
                return(paste(genus, species, subspecies, sep = " "))
            } else {
                return(paste(genus, species, sep = " "))
            }
        }
    }
    # Fallback: abbreviate genus in any case
    parts <- strsplit(tip_label, "_")[[1]]
    if (length(parts) >= 2) {
        #genus <- paste0(substr(parts[1], 1, 1), ".")
        return(paste(parts[1], parts[2], sep = " "))
    }
    # Ultimate fallback: just replace underscores with spaces
    return(gsub("_", " ", tip_label))
}

# Prepare tip labels
tip_data <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE)
cat("Sample tip labels:", head(tip_data$label, 5), "\n")
# Count all labels with their display names
tip_data$display <- sapply(tip_data$label, create_display_name, meta = meta_data)
cat("Sample display labels:", head(tip_data$display, 5), "\n")
# Find unique alternatives
display_counts <- table(tip_data$display)
cat("Number of unique display labels:", length(display_counts), "\n")
cat("Most common display labels:\n")
print(sort(display_counts, decreasing = TRUE)[1:10])
tip_data$genus     <- sub("_.*", "", tip_data$label)
tip_data$is_target <- grepl("^Tityra_leucura", tip_data$label)

# Color palette
genera <- sort(unique(tip_data$genus))
pal    <- setNames(colorRampPalette(brewer.pal(9, "Set1"))(length(genera)), genera)

make_plot <- function(tr, td, layout) {
    ggtree(tr, layout = layout, linewidth = 0.6) %<+% td +
        geom_tiplab(aes(label = display, colour = genus), size = 3, hjust = -0.05) +
        geom_tippoint(aes(colour = genus, shape = is_target, size = is_target)) +
        geom_nodelab(aes(label = label), size = 2.2, colour = "blue",
                     hjust = 1.2, vjust = -0.3) +
        scale_colour_manual(values = pal, name = "Genus") +
        scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), guide = "none") +
        scale_size_manual(values  = c("FALSE" = 1.5, "TRUE" = 4),  guide = "none") +
        geom_treescale(fontsize = 3) +
        theme_tree2() +
        theme(legend.position = "right",
              plot.margin = unit(c(0.5, 8, 0.5, 0.5), "cm")) +
        coord_cartesian(clip = "off")
}

n <- length(tree$tip.label)
basename <- gsub("\\.treefile$", "", basename(tree_file))

p_rect <- make_plot(tree, tip_data, "roundrect")
ggsave(file.path(plot_dir, paste0(basename, "_roundrect.pdf")),
       p_rect, width = 20, height = 0.4 * n + 4, limitsize = FALSE)
ggsave(file.path(plot_dir, paste0(basename, "_roundrect.png")),
       p_rect, width = 20, height = 0.4 * n + 4, dpi = 300, limitsize = FALSE)

cat("Plots written to:", plot_dir, "\n")
cat("  ", basename, "_roundrect.pdf/png\n")