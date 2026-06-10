#!/bin/sh
# ==============================================================================
# Script:       42_alluvial.sh
# Description:  Visualise chromosomal inversion genotype combinations across the
#               seven high-confidence LGCs using an alluvial plot.
#               Each ribbon tracks one individual across LGCs; strata show
#               genotype proportions.
#               Genotype states: A = ancestral, H = heterozygous, D = derived.
#
# Input:        Per-LGC genotype lists in ${base_path} (one individual ID per
#               line): <LGC>_ancestral.txt, <LGC>_het.txt, <LGC>_derived.txt
# Output:       ${base_path}/alluvial_plot.png and .pdf
#
# ==============================================================================

. /etc/profile.d/modules.sh

module load R/4.5

Rscript - <<'EOF'

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggalluvial)
})

# ==============================================================================
# Settings
# ==============================================================================

base_path <- "local/pathr"

# Fixed display order across the plot (high-confidence inversions)
lgc_order <- c("LGC02", "LGC15", "LGC13.01", "LGC17.01",
               "LGC18.01", "LGC18.03", "LGC19.01")

# Base colour per LGC; light (A) / dark (D) tones derived automatically
lgc_base <- c(
  "LGC02"    = "#577590",
  "LGC15"    = "#43aa8b",
  "LGC13.01" = "#90be6d",
  "LGC17.01" = "#f9c74f",
  "LGC18.01" = "#f8961e",
  "LGC18.03" = "#f3722c",
  "LGC19.01" = "#f94144"
)

common_theme <- theme_bw() +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(colour = "black", fill = NA),
    legend.title      = element_text(size = 14, face = "bold"),
    legend.position   = "right",
    legend.text       = element_text(size = 12),
    legend.key.size   = unit(0.6, "cm"),
    axis.title        = element_text(size = 14),
    axis.text         = element_text(size = 12, colour = "black"),
    axis.ticks        = element_line(colour = "black"),
    plot.margin       = margin(10, 10, 10, 10, unit = "points"),
    plot.tag          = element_text(size = 16, face = "bold")
  )

# ==============================================================================
# Section 1: Load and reshape genotype data
# ==============================================================================
# Read each per-LGC genotype list (one individual ID per line) into a long
# table of (category, individual), then pivot to one row per individual with a
# logical column per "<LGC> <state>" category.

x <- list()
for (lgc in lgc_order) {
  x[[paste0(lgc, " A")]] <- readLines(file.path(base_path, paste0(lgc, "_ancestral.txt")))
  x[[paste0(lgc, " H")]] <- readLines(file.path(base_path, paste0(lgc, "_het.txt")))
  x[[paste0(lgc, " D")]] <- readLines(file.path(base_path, paste0(lgc, "_derived.txt")))
}

df_long <- data.frame(
  ID   = rep(names(x), lengths(x)),
  indv = unlist(x, use.names = FALSE),
  stringsAsFactors = FALSE
)
categories <- unique(df_long$ID)

df <- df_long %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = ID, values_from = value, values_fill = FALSE)

cat("Data loaded:", nrow(df), "individuals x",
    length(categories), "genotype categories\n")

# ==============================================================================
# Section 2: Build genotype matrix  (0 = A, 1 = H, 2 = D, NA = missing)
# ==============================================================================

prepare_genotype_matrix <- function(df, lgc_order) {
  mat <- matrix(NA_real_, nrow = nrow(df), ncol = length(lgc_order),
                dimnames = list(df$indv, lgc_order))
  
  for (i in seq_along(lgc_order)) {
    lgc   <- lgc_order[i]
    a_col <- paste0(lgc, " A")
    h_col <- paste0(lgc, " H")
    d_col <- paste0(lgc, " D")
    
    has_value <- df[[a_col]] | df[[h_col]] | df[[d_col]]
    mat[, i]  <- ifelse(!has_value, NA,
                        ifelse(df[[a_col]], 0,
                               ifelse(df[[h_col]], 1,
                                      ifelse(df[[d_col]], 2, NA))))
  }
  mat
}

geno_mat <- prepare_genotype_matrix(df, lgc_order)
cat("Genotype matrix:", nrow(geno_mat), "individuals x", ncol(geno_mat), "LGCs\n")

# ==============================================================================
# Section 3: Build colour palette  (light = A, base = H, dark = D)
# ==============================================================================

build_color_palette <- function(lgc_base, lgc_order) {
  node_colors <- c()
  for (lgc in lgc_order) {
    rgb_vals <- col2rgb(lgc_base[lgc]) / 255
    light <- rgb(pmin(rgb_vals[1] + 0.2, 1),
                 pmin(rgb_vals[2] + 0.2, 1),
                 pmin(rgb_vals[3] + 0.2, 1))
    dark  <- rgb(pmax(rgb_vals[1] - 0.2, 0),
                 pmax(rgb_vals[2] - 0.2, 0),
                 pmax(rgb_vals[3] - 0.2, 0))
    node_colors <- c(node_colors,
                     setNames(c(light, lgc_base[[lgc]], dark),
                              paste0(lgc, c("_A", "_H", "_D"))))
  }
  node_colors
}

node_colors <- build_color_palette(lgc_base, lgc_order)

# ==============================================================================
# Section 4: Figure 4 - Alluvial plot
# ==============================================================================

alluvial_df <- as.data.frame(geno_mat) %>%
  rownames_to_column("indv") %>%
  pivot_longer(-indv, names_to = "LGC", values_to = "genotype") %>%
  mutate(
    genotype = case_when(
      genotype == 0 ~ "A",
      genotype == 1 ~ "H",
      genotype == 2 ~ "D",
      TRUE          ~ NA_character_),
    genotype = factor(genotype, levels = c("A", "H", "D")),
    node_id  = paste0(LGC, "_", as.character(genotype)),
    LGC      = factor(LGC, levels = lgc_order)
  ) %>%
  filter(!is.na(genotype)) %>%
  mutate(indv = factor(indv))

alluvial <- ggplot(
  alluvial_df,
  aes(x = LGC, stratum = genotype, alluvium = indv, fill = node_id)) +
  geom_flow(alpha = 0.4, color = "white", linewidth = 0.2, aes.flow = "forward") +
  geom_stratum(width = 0.3, color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3.5, fontface = "bold", color = "white") +
  scale_fill_manual(values = node_colors) +
  labs(title = "", x = NULL, y = "Number of individuals") +
  theme_minimal(base_size = 12) +
  common_theme +
  theme(
    plot.title         = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position    = "none",
    axis.title.y       = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 12),
    panel.grid.major.x = element_blank())

alluvial

ggsave(file.path(base_path, "fig4_alluvial_plot.png"),
       plot = alluvial, width = 12, height = 8, dpi = 600, bg = "white")
ggsave(file.path(base_path, "fig4_alluvial_plot.pdf"),
       plot = alluvial, width = 12, height = 8, dpi = 600, bg = "white")

cat("Saved: fig4_alluvial_plot.png/.pdf\n")
cat("Done.\n")
