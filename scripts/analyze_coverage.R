#!/usr/bin/env Rscript

# Analyze UCE coverage and sequencing depth distribution
# This script creates histograms with mean values highlighted as red lines

# Load required libraries
library(ggplot2)
library(dplyr)

# Define input and output paths
input_file <- "/media/inter/mkapun/projects/Tityra/results/tityra_uce/mapping/Tityra_leucura_uce_coverage.txt"
output_dir <- "/media/inter/mkapun/projects/Tityra/results/tityra_uce/mapping"

# Read coverage data (header line starts with #, so skip it)
coverage_data <- read.table(input_file, header = TRUE, sep = "\t")

# Extract coverage and meandepth columns
# Note: Some rows may have zero coverage (excluded from alignment), filter those out
coverage_data_filtered <- coverage_data %>%
  filter(coverage > 0, meandepth > 0)

# Calculate statistics
mean_coverage <- mean(coverage_data_filtered$coverage)
median_coverage <- median(coverage_data_filtered$coverage)
mean_depth <- mean(coverage_data_filtered$meandepth)
median_depth <- median(coverage_data_filtered$meandepth)

# Print summary statistics
cat("\n========== Coverage Analysis Summary ==========\n")
cat(sprintf("Total UCE loci analyzed: %d\n", nrow(coverage_data_filtered)))
cat(sprintf("Loci with zero coverage (excluded): %d\n", nrow(coverage_data) - nrow(coverage_data_filtered)))
cat(sprintf("\nCoverage (%%):\n"))
cat(sprintf("  Mean:   %.2f%%\n", mean_coverage))
cat(sprintf("  Median: %.2f%%\n", median_coverage))
cat(sprintf("  Min:    %.2f%%\n", min(coverage_data_filtered$coverage)))
cat(sprintf("  Max:    %.2f%%\n", max(coverage_data_filtered$coverage)))
cat(sprintf("\nMean Depth (x):\n"))
cat(sprintf("  Mean:   %.2fx\n", mean_depth))
cat(sprintf("  Median: %.2fx\n", median_depth))
cat(sprintf("  Min:    %.2fx\n", min(coverage_data_filtered$meandepth)))
cat(sprintf("  Max:    %.2fx\n", max(coverage_data_filtered$meandepth)))
cat("\n==============================================\n\n")

# Create Coverage Distribution Plot
p1 <- ggplot(coverage_data_filtered, aes(x = coverage)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = mean_coverage, color = "red", linewidth = 1.5, linetype = "solid") +
  geom_vline(xintercept = median_coverage, color = "darkred", linewidth = 1, linetype = "dashed") +
  annotate("text", x = mean_coverage, y = Inf, 
           label = sprintf("Mean: %.1f%%", mean_coverage), 
           vjust = 1.5, hjust = 1.2, color = "red", size = 4) +
  annotate("text", x = median_coverage, y = Inf, 
           label = sprintf("Median: %.1f%%", median_coverage), 
           vjust = 2.7, hjust = 1.1, color = "darkred", size = 3.5) +
  labs(title = "Distribution of UCE Coverage",
       x = "Coverage (%)",
       y = "Frequency") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))

# Save coverage plot
ggsave(filename = file.path(output_dir, "uce_coverage_distribution.png"),
       plot = p1, width = 10, height = 6, dpi = 300)

# Create Mean Depth Distribution Plot
p2 <- ggplot(coverage_data_filtered, aes(x = meandepth)) +
  geom_histogram(bins = 50, fill = "darkgreen", color = "black", alpha = 0.7) +
  geom_vline(xintercept = mean_depth, color = "red", linewidth = 1.5, linetype = "solid") +
  geom_vline(xintercept = median_depth, color = "darkred", linewidth = 1, linetype = "dashed") +
  annotate("text", x = mean_depth, y = Inf, 
           label = sprintf("Mean: %.1fx", mean_depth), 
           vjust = 1.5, hjust = 1.2, color = "red", size = 4) +
  annotate("text", x = median_depth, y = Inf, 
           label = sprintf("Median: %.1fx", median_depth), 
           vjust = 2.7, hjust = 1.1, color = "darkred", size = 3.5) +
  labs(title = "Distribution of UCE Mean Sequencing Depth",
       x = "Mean Depth (x)",
       y = "Frequency") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))

# Save depth plot
ggsave(filename = file.path(output_dir, "uce_depth_distribution.png"),
       plot = p2, width = 10, height = 6, dpi = 300)

# Create combined scatter plot showing relationship between coverage and depth
p3 <- ggplot(coverage_data_filtered, aes(x = coverage, y = meandepth)) +
  geom_point(color = "purple", alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = mean_coverage, color = "red", linewidth = 1.5, linetype = "solid") +
  geom_hline(yintercept = mean_depth, color = "red", linewidth = 1.5, linetype = "solid") +
  annotate("label", x = mean_coverage, y = max(coverage_data_filtered$meandepth), 
           label = sprintf("Mean Coverage: %.1f%%", mean_coverage), 
           color = "red", size = 3.5, hjust = 0) +
  annotate("label", x = max(coverage_data_filtered$coverage), y = mean_depth, 
           label = sprintf("Mean Depth: %.1fx", mean_depth), 
           color = "red", size = 3.5, vjust = -0.5) +
  labs(title = "Coverage vs Mean Depth for UCE Loci",
       x = "Coverage (%)",
       y = "Mean Depth (x)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))

# Save scatter plot
ggsave(filename = file.path(output_dir, "uce_coverage_vs_depth.png"),
       plot = p3, width = 10, height = 6, dpi = 300)

cat(sprintf("\nPlots saved to %s:\n", output_dir))
cat("  - uce_coverage_distribution.png\n")
cat("  - uce_depth_distribution.png\n")
cat("  - uce_coverage_vs_depth.png\n")
cat("\nAnalysis complete!\n")