# ==============================================================================
# Script:       14b.pcangsd.R
# Description:  Visualise population structure from PCAngsd covariance matrices.
#
# Input:        Genome-wide .cov file and per-chromosome .cov files from 14a.pcangsd.sh;
#                    pop_info_160.txt sample metadata table
# ==============================================================================

library(data.table)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)

# ── Paths ──────────────────────────────────────────────────────────────────────
COV_GENOME <- "PCAngsd/Pmax_160_merged_chr.cov"
COV_DIR    <- "PCAngsd/chr/"
POP_FILE   <- "PCAngsd/pop_info_160.txt"
OUT_DIR    <- "Figures/"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Colour palette and legend order ───────────────────────────────────────────
col_palette  <- c("#0477BF", "#38E0F2", "#29A655", "#7FBF50",
                  "#F28891", "#FFFE48", "#FFA75C", "#DC143C")
legend_order <- c("SHE", "ORK", "MNC", "MOR", "FRA", "CLY", "TAR", "CDB")

# ── Common ggplot theme ────────────────────────────────────────────────────────
common_theme <- theme_bw() +
  theme(
    panel.grid.major    = element_blank(),
    panel.grid.minor    = element_blank(),
    panel.border        = element_rect(colour = "black", fill = NA),
    legend.title        = element_text(size = 22, face = "bold"),
    legend.position     = "right",
    legend.text         = element_text(size = 20),
    legend.key.size     = unit(0.8, "cm"),
    legend.spacing.y    = unit(0.5, "cm"),
    axis.title          = element_text(size = 22, face = "bold"),
    axis.text           = element_text(size = 20, colour = "black"),
    axis.ticks          = element_line(colour = "black"),
    plot.margin         = margin(20, 20, 20, 20, unit = "points"),
    plot.title          = element_text(size = 24, face = "bold"),
    plot.tag            = element_text(size = 20, face = "bold")
  )

# ── Sample metadata ────────────────────────────────────────────────────────────
pop_info <- read_delim(POP_FILE, delim = "\t",
                       escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(Location = if_else(Location == "MFR", "MOR", Location))

# ==============================================================================
# PART 1: Genome-wide PCA
# ==============================================================================

cat("\n[Genome-wide PCA]\n")

cov <- as.matrix(read.table(COV_GENOME))
e   <- eigen(cov)

pca_sum     <- sum(e$values)
var_exp     <- (e$values[1:4] / pca_sum) * 100
names(var_exp) <- paste0("PC", 1:4)
cat("Variance explained (PC1–PC4):\n"); print(round(var_exp, 2))

pca_df <- as.data.frame(e$vectors[, 1:3]) %>%
  setNames(c("PC1", "PC2", "PC3")) %>%
  bind_cols(pop_info)

# Helper: create a single PCA scatter plot
make_pca_plot <- function(df, xpc, ypc, vx, vy) {
  ggplot(df, aes(x = .data[[xpc]], y = .data[[ypc]],
                 fill  = factor(Location, levels = legend_order),
                 label = Sample)) +
    geom_point(pch = 21, size = 4, colour = "#000033", alpha = 0.9) +
    scale_fill_manual(values = col_palette) +
    labs(x    = paste0(xpc, " (", round(vx, 1), "%)"),
         y    = paste0(ypc, " (", round(vy, 1), "%)"),
         fill = "Locations") +
    common_theme
}

p1 <- make_pca_plot(pca_df, "PC1", "PC2", var_exp[1], var_exp[2]) + labs(tag = "A")
p2 <- make_pca_plot(pca_df, "PC2", "PC3", var_exp[2], var_exp[3]) + labs(tag = "B")

print(p1 + p2 + plot_layout(guides = "collect"))

# Save genome-wide PCA
genome_plot <- p1 + p2 + plot_layout(guides = "collect")
ggsave(file.path(OUT_DIR, "Pmax_160_pca.pdf"), genome_plot, width = 12, height = 6, dpi = 300)
cat("Saved:", file.path(OUT_DIR, "Pmax_160_pca.pdf"), "\n")

# ==============================================================================
# PART 2: Per-chromosome PCA panel
# ==============================================================================

cat("\n[Per-chromosome PCA]\n")

cov_files <- sort(list.files(COV_DIR, pattern = "\\.cov$", full.names = TRUE))
cat("Chromosome .cov files found:", length(cov_files), "\n")

# Minimal theme for faceted chromosome panels
chr_theme <- theme_cowplot() +
  theme(
    legend.position  = "none",
    panel.grid       = element_blank(),
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    strip.text.x     = element_text(size = 20, face = "bold"),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

# Build one plot per chromosome
chr_plots <- vector("list", length(cov_files))

for (i in seq_along(cov_files)) {

  cov_chr <- as.matrix(fread(cov_files[[i]], header = FALSE))
  e_chr   <- eigen(cov_chr)

  pca_chr <- as.data.frame(e_chr$vectors[, 1:2]) %>%
    setNames(c("PC1", "PC2")) %>%
    bind_cols(pop_info) %>%
    mutate(chromosome = paste0("chr ", i))

  chr_plots[[i]] <- ggplot(pca_chr,
      aes(x = PC1, y = PC2,
          fill  = factor(Location, levels = legend_order),
          label = Sample)) +
    geom_point(pch = 21, size = 4, colour = "#000033", alpha = 0.9) +
    scale_fill_manual(values = col_palette) +
    labs(fill = "Locations") +
    facet_wrap(~chromosome) +
    chr_theme
}

# Combine chromosome panels; add shared axis labels and legend
chr_combined <- wrap_plots(chr_plots) +
  plot_layout(ncol = 5, byrow = TRUE, guides = "collect") &
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 20),
    legend.title    = element_text(size = 20, face = "bold")
  )

# Add shared axis labels via patchwork annotation
chr_combined <- wrap_elements(chr_combined) +
  labs(tag = "PC1") +
  theme(plot.tag.position = c(0.5, 0),
        plot.tag          = element_text(size = 24, face = "bold", vjust = 1))

chr_combined

ggsave(file.path(OUT_DIR, "Pmax_160_pca_chr.png"), chr_combined, width = 20, height = 16, dpi = 300)
cat("Saved:", file.path(OUT_DIR, "Pmax_160_pca_chr.png"), "\n")

cat("\nDone. Figures saved to:", OUT_DIR, "\n")
