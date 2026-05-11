#!/bin/bash
###############################################################################
# Visualize complete mitochondrial genome alignments
# Author: Generated for Tityra project
# Description: Create visualization of mitogenome alignment
###############################################################################

WD="/media/inter/mkapun/projects/Tityra"

# Create output directory
mkdir -p "${WD}/results/complete_mitogenome_alignment/plots"

# Activate environment
conda activate /media/inter/mkapun/projects/MuseomicsWorkshop2025/scripts/programs

# Check which alignment file exists
alignment_file_annotated="${WD}/results/complete_mitogenome_alignment/annotation_aware_aligned.fasta"
alignment_file_whole="${WD}/results/complete_mitogenome_alignment/whole_genome_aligned.fasta"

# Visualize both if they exist
ALIGN_FILES=()
ALIGN_NAMES=()

if [ -f "${alignment_file_annotated}" ]; then
    ALIGN_FILES+=("${alignment_file_annotated}")
    ALIGN_NAMES+=("annotation_aware")
    echo "Found: annotation-aware alignment"
fi

if [ -f "${alignment_file_whole}" ]; then
    ALIGN_FILES+=("${alignment_file_whole}")
    ALIGN_NAMES+=("whole_genome")
    echo "Found: whole-genome alignment"
fi

if [ ${#ALIGN_FILES[@]} -eq 0 ]; then
    echo "Error: No alignment files found"
    exit 1
fi

echo "Will create visualizations for ${#ALIGN_FILES[@]} alignment(s)"
echo ""

for i in "${!ALIGN_FILES[@]}"; do
    alignment_file="${ALIGN_FILES[$i]}"
    alignment_name="${ALIGN_NAMES[$i]}"
    
    echo "Processing: ${alignment_name}"

# Create R script for visualization
cat > "${WD}/scripts/visualize_mitogenome_alignment.R" << 'EOF'
#!/usr/bin/env Rscript

###############################################################################
# Visualize Complete Mitochondrial Genome Alignment
###############################################################################

library(ggplot2)
library(Biostrings)
library(reshape2)

# Set working directory
setwd("/media/inter/mkapun/projects/Tityra")

# Get alignment file path and name from command args
args <- commandArgs(trailingOnly = TRUE)
alignment_file <- args[1]
alignment_name <- args[2]

cat(sprintf("\nProcessing alignment: %s\n", alignment_name))

alignment <- readDNAMultipleAlignment(alignment_file, format = "fasta")

cat("Alignment dimensions:\n")
cat("  Sequences:", nrow(alignment), "\n")
cat("  Length:", ncol(alignment), "\n")

# Create output directory
dir.create("results/complete_mitogenome_alignment/plots", showWarnings = FALSE, recursive = TRUE)

###############################################################################
# Plot 1: Base composition across alignment
###############################################################################

cat("\nGenerating base composition plot...\n")

# Calculate base frequencies across alignment
alignment_length <- ncol(alignment)
n_sequences <- nrow(alignment)

# Sample positions for plotting (every 100th position)
sample_freq <- max(1, floor(alignment_length / 1000))
positions <- seq(1, alignment_length, by = sample_freq)

# Calculate base composition
cons_matrix <- consensusMatrix(alignment, as.prob = TRUE)

base_comp <- data.frame(
    Position = positions,
    A = cons_matrix["A", positions],
    C = cons_matrix["C", positions],
    G = cons_matrix["G", positions],
    T = cons_matrix["T", positions]
)

base_comp_long <- melt(base_comp, id.vars = "Position", 
                       variable.name = "Base", value.name = "Frequency")

p1 <- ggplot(base_comp_long, aes(x = Position, y = Frequency, color = Base)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = c("A" = "#64B200", "C" = "#FFB300", 
                                   "G" = "#DC267F", "T" = "#3B7EA1")) +
    labs(
        title = "Base Composition across Complete Mitochondrial Genome Alignment",
        subtitle = paste0(n_sequences, " species, ", alignment_length, " bp"),
        x = "Position (bp)",
        y = "Frequency"
    ) +
    theme_bw(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom"
    )

ggsave("results/complete_mitogenome_alignment/plots/base_composition.pdf",
       p1, width = 14, height = 8)
ggsave("results/complete_mitogenome_alignment/plots/base_composition.png",
       p1, width = 14, height = 8, dpi = 300)

# Update filenames with alignment type
ggsave(sprintf("results/complete_mitogenome_alignment/plots/base_composition_%s.pdf", alignment_name),
       p1, width = 14, height = 8)
ggsave(sprintf("results/complete_mitogenome_alignment/plots/base_composition_%s.png", alignment_name),
       p1, width = 14, height = 8, dpi = 300)

###############################################################################
# Plot 2: Conservation plot across alignment
###############################################################################

cat("\nCalculating conservation scores...\n")

# Calculate conservation using consensus matrix
cons_matrix <- consensusMatrix(alignment, as.prob = TRUE)

# Calculate Shannon entropy for each position (measure of conservation)
shannon_entropy <- apply(cons_matrix[1:4, ], 2, function(x) {
    x <- x[x > 0]
    if (length(x) == 0) return(2)
    -sum(x * log2(x))
})

# Convert to conservation score (lower entropy = higher conservation)
conservation <- 2 - shannon_entropy  # Max entropy for 4 bases is 2

# Create data frame for plotting
conservation_df <- data.frame(
    Position = 1:length(conservation),
    Conservation = conservation
)

# Calculate sliding window average for smoother plot
window <- 100
conservation_df$Conservation_smooth <- stats::filter(conservation_df$Conservation, 
                                                      rep(1/window, window), 
                                                      sides = 2)

p3 <- ggplot(conservation_df, aes(x = Position, y = Conservation_smooth)) +
    geom_line(color = "steelblue", linewidth = 0.5) +
    labs(
        title = "Conservation across Complete Mitochondrial Genome",
        subtitle = paste0(n_sequences, " species aligned, smoothed over ", window, " bp windows"),
        x = "Position (bp)",
        y = "Conservation Score"
    ) +
    theme_bw(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
    )

ggsave("results/complete_mitogenome_alignment/plots/conservation_plot.pdf",
       p3, width = 14, height = 6)
ggsave("results/complete_mitogenome_alignment/plots/conservation_plot.png",
       p3, width = 14, height = 6, dpi = 300)

ggsave(sprintf("results/complete_mitogenome_alignment/plots/conservation_plot_%s.pdf", alignment_name),
       p3, width = 14, height = 6)
ggsave(sprintf("results/complete_mitogenome_alignment/plots/conservation_plot_%s.png", alignment_name),
       p3, width = 14, height = 6, dpi = 300)

###############################################################################
# Plot 3: Pairwise identity matrix
###############################################################################

cat("\nCalculating pairwise identities...\n")

# Calculate pairwise identity matrix (optimized for speed)
seq_names <- rownames(alignment)
n_seqs <- length(seq_names)
identity_matrix <- matrix(0, n_seqs, n_seqs, dimnames = list(seq_names, seq_names))

# Convert alignment to character matrix once
aln_matrix <- as.matrix(alignment)

for (i in 1:n_seqs) {
    identity_matrix[i, i] <- 100
    if (i < n_seqs) {
        for (j in (i+1):n_seqs) {
            seq1 <- aln_matrix[i, ]
            seq2 <- aln_matrix[j, ]
            # Only compare non-gap positions
            valid <- seq1 != "-" & seq2 != "-"
            if (sum(valid) > 0) {
                matches <- sum(seq1[valid] == seq2[valid])
                identity <- (matches / sum(valid)) * 100
            } else {
                identity <- 0
            }
            identity_matrix[i, j] <- identity
            identity_matrix[j, i] <- identity
        }
    }
    if (i %% 3 == 0) cat("  Processed", i, "/", n_seqs, "sequences\n")
}

# Convert to data frame for plotting
library(reshape2)
identity_df <- melt(identity_matrix)
colnames(identity_df) <- c("Species1", "Species2", "Identity")

p4 <- ggplot(identity_df, aes(x = Species1, y = Species2, fill = Identity)) +
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green",
                         midpoint = 90, limits = c(80, 100)) +
    labs(
        title = "Pairwise Sequence Identity Matrix",
        subtitle = "Complete mitochondrial genomes",
        x = "", y = "",
        fill = "% Identity"
    ) +
    theme_bw(base_size = 10) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8)
    )

ggsave("results/complete_mitogenome_alignment/plots/pairwise_identity_matrix.pdf",
       p4, width = 12, height = 11)
ggsave("results/complete_mitogenome_alignment/plots/pairwise_identity_matrix.png",
       p4, width = 12, height = 11, dpi = 300)

ggsave(sprintf("results/complete_mitogenome_alignment/plots/pairwise_identity_matrix_%s.pdf", alignment_name),
       p4, width = 12, height = 11)
ggsave(sprintf("results/complete_mitogenome_alignment/plots/pairwise_identity_matrix_%s.png", alignment_name),
       p4, width = 12, height = 11, dpi = 300)

cat("\n=== Visualization Complete ===\n")
cat(sprintf("Plots saved for: %s\n", alignment_name))
EOF

# Run R script for each alignment
Rscript "${WD}/scripts/visualize_mitogenome_alignment.R" "${alignment_file}" "${alignment_name}"

done

echo ""
echo "Visualization complete!"
echo "Plots saved in: ${WD}/results/complete_mitogenome_alignment/plots/"

conda deactivate
