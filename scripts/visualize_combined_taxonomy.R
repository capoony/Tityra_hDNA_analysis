#!/usr/bin/env Rscript

# Combined visualization: Kraken2, ECMSD proportions, and read length distributions
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(cowplot)

# Get working directory from command line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript visualize_combined_taxonomy.R <working_directory>")
}
wd <- args[1]

################################################################################
# Panel A: Kraken2 Results (PE vs Merged comparison, like kraken_genus_comparison.pdf)
################################################################################

kraken_data <- read.csv(file.path(wd, "results/kraken2/kraken_summary.csv"), header = TRUE)

# Parse taxonomic information and extract phylum and genus
kraken_data <- kraken_data %>%
  mutate(
    Phylum = sapply(strsplit(as.character(Taxon), "; "), function(x) x[1]),
    Species = sapply(strsplit(as.character(Taxon), "; "), function(x) tail(x, 1)),
    Genus = sapply(strsplit(as.character(Species), " "), function(x) x[1])
  )

# Aggregate by genus and phylum for both PE and merged
genus_summary <- kraken_data %>%
  group_by(Genus, Phylum) %>%
  summarise(
    report_PE = sum(report_PE, na.rm = TRUE),
    report_merged = sum(report_merged, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Combined_Total = report_PE + report_merged)

# Get top 15 genera by combined total
top15_genera <- genus_summary %>%
  arrange(desc(Combined_Total)) %>%
  head(15) %>%
  pull(Genus)

# Reshape for faceting by sample type
sample_data <- genus_summary %>%
  filter(Genus %in% top15_genera) %>%
  pivot_longer(cols = c(report_merged, report_PE),
               names_to = "Sample",
               values_to = "Percentage") %>%
  mutate(Sample = factor(Sample, 
                        levels = c("report_merged", "report_PE"),
                        labels = c("Merged", "Paired-End")))

# Order genera by combined total (ascending for coord_flip to show descending)
genus_order <- genus_summary %>%
  filter(Genus %in% top15_genera) %>%
  arrange(Combined_Total) %>%
  pull(Genus)

sample_data <- sample_data %>%
  mutate(Genus = factor(Genus, levels = genus_order))

# Get unique phyla and assign viridis colors
phyla <- unique(sample_data$Phylum)
n_phyla <- length(phyla)
phylum_colors <- setNames(
  viridis(n_phyla, option = "plasma"),
  phyla
)

# Panel A: Faceted PE vs Merged comparison, colored by phylum
p_kraken <- ggplot(sample_data, aes(x = Percentage, y = Genus, fill = Phylum)) +
  geom_col() +
  scale_fill_manual(values = phylum_colors) +
  facet_wrap(~ Sample, scales = "fixed") +
  labs(title = "A) Kraken2: Top Genera by Sample Type",
       x = "Percentage (%)",
       y = "Genus",
       fill = "Phylum") +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0, size = 18),
    axis.text.y = element_text(size = 14, face = "italic"),
    axis.text.x = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(face = "bold", size = 16)
  )

################################################################################
# Panel B: ECMSD Proportions (Top 20 genera)
################################################################################

ecmsd_prop <- read.table(file.path(wd, "results/ECMSD/mapping/Mito_summary_genus_proportions.txt"), 
                         header = TRUE, sep = "\t")

# Get top 20 genera by proportion
ecmsd_top20 <- ecmsd_prop %>%
  arrange(desc(Proportion)) %>%
  head(20) %>%
  mutate(Percentage = Proportion * 100)

p_ecmsd_prop <- ggplot(ecmsd_top20, aes(x = reorder(genus, Percentage), 
                                         y = Percentage, 
                                         fill = genus)) +
  geom_col() +
  scale_fill_viridis_d(option = "viridis") +
  coord_flip() +
  labs(title = "B) ECMSD: Top 20 Genera by Proportion",
       x = "Genus",
       y = "Percentage (%)") +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0, size = 18),
    axis.text.y = element_text(size = 14, face = "italic"),
    axis.text.x = element_text(size = 14),
    axis.title = element_text(size = 16)
  )

################################################################################
# Panel C: ECMSD Read Length Distribution (Top 10 Genera, FACETED)
################################################################################

ecmsd_lengths <- read.table(file.path(wd, "results/ECMSD/mapping/Mito_summary_genus.txt"), 
                            header = TRUE, sep = "\t")

# Get top 10 genera by total reads
top10_genera <- ecmsd_lengths %>%
  group_by(genus) %>%
  summarise(TotalReads = sum(TotalReads, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(TotalReads)) %>%
  head(10) %>%
  pull(genus)

# Filter for top 10 genera and set factor order
ecmsd_lengths_top10 <- ecmsd_lengths %>%
  filter(genus %in% top10_genera) %>%
  mutate(genus = factor(genus, levels = top10_genera))

# Panel C: Faceted barplot by genus (like ECMSD pipeline)
p_ecmsd_lengths <- ggplot(ecmsd_lengths_top10, aes(x = Length, 
                                                    y = TotalReads, 
                                                    color = genus)) +
  geom_col(width = 1) +
  facet_wrap(~ genus, scales = "free_y", ncol = 5) +
  scale_color_viridis_d(option = "viridis") +
  labs(title = "C) ECMSD: Read Length Distribution (Top 10 Genera)",
       x = "Read Length (bp)",
       y = "Total Reads") +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0, size = 18),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14, face = "bold.italic")
  )

################################################################################
# Combine and Save: A and B in left column, C on right
################################################################################

# Create left column with A and B stacked (aligned)
left_column <- plot_grid(p_kraken, p_ecmsd_prop, ncol = 1, rel_heights = c(1, 1), align = "v", axis = "lr")

# Combine left column with C on the right
combined_plot <- plot_grid(left_column, p_ecmsd_lengths, ncol = 2, rel_widths = c(1, 1.5))

# Save combined figure
ggsave(file.path(wd, "results/taxonomy_combined.pdf"), 
       combined_plot, width = 20, height = 12, dpi = 300)
ggsave(file.path(wd, "results/taxonomy_combined.png"), 
       combined_plot, width = 20, height = 12, dpi = 300)

cat("\n=== Combined Taxonomy Visualization Complete ===\n")
cat("Combined figure saved:\n")
cat("  - results/taxonomy_combined.pdf\n")
cat("  - results/taxonomy_combined.png\n")

