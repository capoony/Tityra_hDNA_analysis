#!/usr/bin/env Rscript

###############################################################################
# GBIF Occurrence Map for Tityra inquisitor (ssp.) and Tityra leucura
# 
# - Tityra inquisitor albitorques: blue circles
# - Tityra inquisitor pelzelni:   blue triangles
# - Tityra leucura:               red diamonds
###############################################################################

library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)

# --- Parse command line arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("Usage: Rscript plot_GBIF_map.R <working_directory>")
}
wd <- args[1]

# --- Read GBIF data ---
dat <- read.table(file.path(wd, "data/GBIF_tityra.txt"),
    header = TRUE, sep = "\t", quote = "",
    comment.char = "", fill = TRUE, stringsAsFactors = FALSE)

cat("Total records:", nrow(dat), "\n")

# --- Filter records with valid coordinates ---
dat <- dat %>%
    filter(!is.na(decimalLatitude) & !is.na(decimalLongitude)) %>%
    filter(decimalLatitude != "" & decimalLongitude != "")

cat("Records with coordinates:", nrow(dat), "\n")

# --- Assign grouping variable ---
# Tityra leucura: species == "Tityra leucura"
# Tityra inquisitor albitorques: species == "Tityra inquisitor" & infraspecificEpithet == "albitorques"
# Tityra inquisitor pelzelni:   species == "Tityra inquisitor" & infraspecificEpithet == "pelzelni"
# Some records of T. inquisitor lack infraspecificEpithet; we treat them as "unknown ssp."

dat <- dat %>%
    mutate(
        Group = case_when(
            species == "Tityra leucura"                 ~ "Tityra leucura",
            species == "Tityra inquisitor" &
                infraspecificEpithet == "albitorques"   ~ "T. inquisitor albitorques",
            species == "Tityra inquisitor" &
                infraspecificEpithet == "pelzelni"      ~ "T. inquisitor pelzelni",
            species == "Tityra inquisitor"              ~ "T. inquisitor (ssp. unknown)",
            TRUE                                        ~ species
        )
    )

# --- Summarise counts ---
counts <- dat %>% count(Group)
cat("\nOccurrence counts:\n")
print(counts)

# --- Download country data for basemap ---
countries <- ne_countries(scale = "medium", returnclass = "sf")

# --- Define shapes and colors ---
# Shape mapping (pch values):
#   16 = filled circle  (albitorques)
#   17 = filled triangle (pelzelni)
#   18 = diamond         (leucura)
#   1  = open circle     (unknown ssp.)
shape_vals <- c(
    "T. inquisitor albitorques"  = 16,
    "T. inquisitor pelzelni"     = 17,
    "Tityra leucura"             = 18,
    "T. inquisitor (ssp. unknown)" = 1
)

color_vals <- c(
    "T. inquisitor albitorques"  = "#2166AC",
    "T. inquisitor pelzelni"     = "#2166AC",
    "Tityra leucura"             = "#D73027",
    "T. inquisitor (ssp. unknown)" = "grey50"
)

# --- Filter groups actually present ---
present_groups <- intersect(names(shape_vals), unique(dat$Group))

shape_vals <- shape_vals[present_groups]
color_vals <- color_vals[present_groups]

# --- Map ---
p <- ggplot() +
    # Basemap
    geom_sf(data = countries, fill = "grey90", colour = "grey60", size = 0.2) +
    # Add a thin bounding box area
    coord_sf(
        xlim = c(-110, -30),
        ylim = c(-30, 30),
        expand = FALSE
    ) +
    # Points / symbols
    geom_point(
        data = dat %>% filter(Group %in% present_groups),
        aes(x = decimalLongitude, y = decimalLatitude,
            shape = Group, colour = Group),
        size = 2.5,
        alpha = 0.8
    ) +
    scale_shape_manual(values = shape_vals) +
    scale_colour_manual(values = color_vals) +
    # Labels
    labs(
        title = "GBIF Occurrences: Tityra inquisitor and Tityra leucura",
        x = "Longitude", y = "Latitude",
        shape = NULL, colour = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "bottom",
        legend.text = element_text(size = 12, face = "italic")
    ) +
    guides(
        shape = guide_legend(nrow = 2, override.aes = list(size = 3))
    )

# --- Save ---
outdir <- file.path(wd, "results")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(outdir, "GBIF_tityra_map.pdf"), p, width = 10, height = 8, dpi = 300)
ggsave(file.path(outdir, "GBIF_tityra_map.png"), p, width = 10, height = 8, dpi = 300)

cat("\nMap saved to:\n")
cat("  - results/GBIF_tityra_map.pdf\n")
cat("  - results/GBIF_tityra_map.png\n")
cat("\nDone!\n")