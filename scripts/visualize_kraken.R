#!/usr/bin/env Rscript

# Visualize Kraken2 results at genus level with phyla highlighted
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Get working directory from command line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("Usage: Rscript visualize_kraken.R <working_directory>")
}
wd <- args[1]

# Read the kraken summary data
data <- read.csv(file.path(wd, "results/kraken2/kraken_summary.csv"), header = TRUE)

# Parse taxonomic information - extract phylum and genus from species name
data <- data %>%
    mutate(
        Phylum = sapply(strsplit(as.character(Taxon), "; "), function(x) x[1]),
        Species = sapply(strsplit(as.character(Taxon), "; "), function(x) tail(x, 1)),
        # Extract genus from species (first word)
        Genus = sapply(strsplit(as.character(Species), " "), function(x) x[1])
    ) %>%
    # Convert to long format for plotting
    pivot_longer(
        cols = c(report_PE, report_merged),
        names_to = "Sample",
        values_to = "Percentage"
    ) %>%
    # Remove NA values
    filter(!is.na(Percentage)) %>%
    # Clean sample names
    mutate(Sample = gsub("report_", "", Sample))

# Aggregate by genus across all species within each genus
genus_summary <- data %>%
    group_by(Genus, Phylum, Sample) %>%
    summarise(Total_Percentage = sum(Percentage, na.rm = TRUE), .groups = "drop")

# For plot 1: Average across samples
genus_data <- genus_summary %>%
    group_by(Genus, Phylum) %>%
    summarise(Mean_Percentage = mean(Total_Percentage, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(Mean_Percentage)) %>%
    top_n(20, Mean_Percentage)

# Set color palette for phyla
n_phyla <- length(unique(genus_data$Phylum))
phylum_colors <- setNames(
    colorRampPalette(brewer.pal(min(n_phyla, 12), "Set3"))(n_phyla),
    unique(genus_data$Phylum)
)

# Create bar plot - average across samples
p1 <- ggplot(genus_data, aes(
    x = reorder(Genus, Mean_Percentage),
    y = Mean_Percentage,
    fill = Phylum
)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
    scale_fill_manual(values = phylum_colors) +
    coord_flip() +
    labs(
        title = "Top 20 Genera in Kraken2 Classification",
        subtitle = "Aggregated by Genus, colored by Phylum",
        x = "Genus",
        y = "Mean Percentage (%)"
    ) +
    theme_bw(base_size = 24) +
    theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 20, face = "italic")
    )

# Save plot
ggsave(file.path(wd, "results/kraken2/kraken_genus_barplot.pdf"), p1, width = 12, height = 8)
ggsave(file.path(wd, "results/kraken2/kraken_genus_barplot.png"), p1, width = 12, height = 8, dpi = 300)

# For plot 2: Separate panels for PE and merged with same genera and scale
# Get top genera based on combined total across both samples
top_genera_list <- genus_summary %>%
    group_by(Genus, Phylum) %>%
    summarise(Combined_Total = sum(Total_Percentage, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(Combined_Total)) %>%
    top_n(15, Combined_Total)

# Filter data to only include top genera (Phylum already in genus_summary)
sample_data <- genus_summary %>%
    filter(Genus %in% top_genera_list$Genus) %>%
    mutate(Sample = factor(Sample,
        levels = c("merged", "PE"),
        labels = c("Merged", "Paired-End")
    ))

# Create ordered factor for consistent y-axis (descending order - highest at top)
genus_order <- top_genera_list %>%
    arrange(Combined_Total) %>%
    pull(Genus)
sample_data <- sample_data %>%
    mutate(Genus = factor(Genus, levels = genus_order))

p2 <- ggplot(sample_data, aes(
    x = Total_Percentage,
    y = Genus,
    fill = Phylum
)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
    scale_fill_manual(values = phylum_colors) +
    facet_wrap(~Sample, scales = "fixed") +
    labs(
        title = "Top Genera by Sample Type",
        subtitle = "Aggregated by Genus (Merged vs Paired-End Reads)",
        x = "Total Percentage (%)",
        y = "Genus"
    ) +
    theme_bw(base_size = 22) +
    theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 18, face = "italic"),
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold")
    )

ggsave(file.path(wd, "results/kraken2/kraken_genus_comparison.pdf"), p2, width = 14, height = 8)
ggsave(file.path(wd, "results/kraken2/kraken_genus_comparison.png"), p2, width = 14, height = 8, dpi = 300)

cat(paste0("Kraken2 visualization completed. Plots saved to ", wd, "/results/kraken2/\n"))
