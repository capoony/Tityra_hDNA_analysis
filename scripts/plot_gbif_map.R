#!/usr/bin/env Rscript

###############################################################################
# GBIF Occurrence Map for Tityra inquisitor (two subspecies) and Tityra leucura
# - Tityra inquisitor albitorques: blue circles
# - Tityra inquisitor pelzelni:    blue triangles
# - Tityra leucura:               red squares
###############################################################################

library(ggplot2)
library(dplyr)
library(maps)
library(viridis)

# Get working directory from command line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("Usage: Rscript plot_gbif_map.R <working_directory>")
}
wd <- args[1]

################################################################################
# Read and filter GBIF data
################################################################################

gbif <- read.table(
    file.path(wd, "data/GBIF_tityra.txt"),
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
)

cat("Total records:", nrow(gbif), "\n")

# Filter for records that have valid coordinates
gbif_coord <- gbif %>%
    filter(
        !is.na(decimalLatitude),
        !is.na(decimalLongitude),
        decimalLatitude != "",
        decimalLongitude != "",
        decimalLatitude >= -60,
        decimalLatitude <= 30,
        decimalLongitude >= -120,
        decimalLongitude <= -30
    ) %>%
    mutate(
        decimalLatitude  = as.numeric(decimalLatitude),
        decimalLongitude = as.numeric(decimalLongitude)
    ) %>%
    filter(
        !is.na(decimalLatitude),
        !is.na(decimalLongitude)
    )

cat("Records with valid coordinates:", nrow(gbif_coord), "\n")

################################################################################
# Assign groups for plotting
################################################################################

# Create a grouping column
gbif_coord <- gbif_coord %>%
    mutate(
        Group = case_when(
            species == "Tityra leucura" ~ "Tityra leucura",
            species == "Tityra inquisitor" & infraspecificEpithet == "albitorques" ~ "Tityra inquisitor albitorques",
            species == "Tityra inquisitor" & infraspecificEpithet == "pelzelni" ~ "Tityra inquisitor pelzelni",
            species == "Tityra inquisitor" & infraspecificEpithet == "" ~ "Tityra inquisitor (subsp. unknown)",
            TRUE ~ "Other"
        )
    )

# Print summary
cat("\nRecords per group:\n")
print(table(gbif_coord$Group))

################################################################################
# Get world map data
################################################################################

world_map <- map_data("world", 
    regions = c(
        "Mexico", "Guatemala", "Belize", "Honduras", "El Salvador",
        "Nicaragua", "Costa Rica", "Panama", "Colombia", "Venezuela",
        "Guyana", "Suriname", "French Guiana", "Ecuador", "Peru",
        "Brazil", "Bolivia", "Paraguay", "Uruguay", "Argentina",
        "Chile", "Trinidad and Tobago"
    )
)

# Also get full world map for broader context
world_map_full <- map_data("world")

################################################################################
# Define colors and shapes
################################################################################

# Colors: blue shades for T. inquisitor, red for T. leucura
group_colors <- c(
    "Tityra inquisitor albitorques"     = "#1a6b8a",  # medium blue
    "Tityra inquisitor pelzelni"        = "#4da6c4",  # light blue
    "Tityra inquisitor (subsp. unknown)" = "#999999",  # grey
    "Tityra leucura"                    = "#d73027"   # red
)

# Shapes: circles for albitorques, triangles for pelzelni, squares for leucura
group_shapes <- c(
    "Tityra inquisitor albitorques"      = 16,   # filled circle
    "Tityra inquisitor pelzelni"         = 17,   # filled triangle
    "Tityra inquisitor (subsp. unknown)" = 1,    # open circle
    "Tityra leucura"                     = 15    # filled square
)

# Alpha/transparency
group_alpha <- c(
    "Tityra inquisitor albitorques"      = 0.5,
    "Tityra inquisitor pelzelni"         = 0.5,
    "Tityra inquisitor (subsp. unknown)" = 0.3,
    "Tityra leucura"                     = 0.8
)

# Size
group_size <- c(
    "Tityra inquisitor albitorques"      = 1.5,
    "Tityra inquisitor pelzelni"         = 1.5,
    "Tityra inquisitor (subsp. unknown)" = 0.8,
    "Tityra leucura"                     = 2.5
)

################################################################################
# Create the map
################################################################################

p <- ggplot() +
    # Base map
    geom_polygon(
        data = world_map_full,
        aes(x = long, y = lat, group = group),
        fill = "grey85",
        colour = "grey50",
        linewidth = 0.2
    ) +
    # Add country borders for focus region
    geom_polygon(
        data = world_map,
        aes(x = long, y = lat, group = group),
        fill = "grey90",
        colour = "grey40",
        linewidth = 0.3
    ) +
    # Occurrence points
    geom_point(
        data = gbif_coord,
        aes(
            x = decimalLongitude,
            y = decimalLatitude,
            color = Group,
            shape = Group,
            alpha = Group,
            size = Group
        )
    ) +
    # Set manual scales
    scale_color_manual(
        values = group_colors,
        name = NULL
    ) +
    scale_shape_manual(
        values = group_shapes,
        name = NULL
    ) +
    scale_alpha_manual(
        values = group_alpha,
        name = NULL
    ) +
    scale_size_manual(
        values = group_size,
        name = NULL
    ) +
    # Coordinate system (Central/South America)
    coord_fixed(
        xlim = c(-90, -30),
        ylim = c(-30, 20),
        ratio = 1.3
    ) +
    # Labels
    labs(
        title = "GBIF Occurrence Records",
        subtitle = expression(italic("Tityra inquisitor") ~ "(two subspecies) and" ~ italic("Tityra leucura")),
        x = "Longitude",
        y = "Latitude"
    ) +
    # Theme
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 13, hjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(size = 12, face = "italic"),
        legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey80", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "aliceblue", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)
    ) +
    # Add legend guide with combined color/shape/alpha/size
    guides(
        color = guide_legend(
            override.aes = list(
                size = c(3, 3, 2, 4),
                alpha = 1
            )
        )
    )

################################################################################
# Save the map
################################################################################

output_dir <- file.path(wd, "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
    file.path(output_dir, "GBIF_Tityra_map.pdf"),
    plot = p,
    width = 10,
    height = 12,
    dpi = 300
)

ggsave(
    file.path(output_dir, "GBIF_Tityra_map.png"),
    plot = p,
    width = 10,
    height = 12,
    dpi = 300
)

cat("\n=== GBIF Occurrence Map Complete ===\n")
cat("Output saved to:\n")
cat("  - results/GBIF_Tityra_map.pdf\n")
cat("  - results/GBIF_Tityra_map.png\n")
cat("\nRecords plotted:\n")
print(table(gbif_coord$Group))