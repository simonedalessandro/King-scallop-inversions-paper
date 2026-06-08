#!/bin/sh
# ==============================================================================
# Script:       25.double_inversion_plot.sh
# Description:  Multi-panel figure for the LGC18.01_02 double inversion, which
#               resolves into five local-PCA clusters (rather than the usual
#               three karyotypes).
#
#               Panel layout:
#                 A — genome-wide MDS (axis 14) with inversion window highlighted
#                 B — local PCA coloured by the five k-means clusters
#                 C — per-cluster heterozygosity (dot + boxplot)
#
# Input:        From the LGC18.01_02 working directory:
#                 mds_matrix.txt        — per-window MDS values
#                 sample.ID.txt         — sample IDs (--012 row order)
#                 cluster_<1-5>.list.txt — k-means cluster membership lists
#                 LGC18.01_02_pca_coords.txt — local PCA coordinates
#                 LGC18.01_02.het       — per-individual heterozygosity
#
# Output:       LGC18.01_02.png and LGC18.01_02.pdf
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE job configuration
# ------------------------------------------------------------------------------
#$ -N LGC18_plot
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_rss=40G
#$ -pe sharedmem 6
#$ -o o_files
#$ -e e_files

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load R/4.5

# ------------------------------------------------------------------------------
# Run R
# ------------------------------------------------------------------------------

Rscript - <<'EOF'

library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggdist)
library(data.table)
library(ggrastr)

# ── Paths ──────────────────────────────────────────────────────────────────────
base_dir <- "/simone/pmax2023/out/20.lostruct/inversions/LGC18.01_02"
LGC      <- "LGC18.01_02"

message("==> Processing: ", LGC)

# ── Colour palette ─────────────────────────────────────────────────────────────
# One colour per local-PCA cluster (1–5)
col_palette <- c(
  "1" = "#1B3A5C",
  "2" = "#2E86AB",
  "3" = "#E07B39",
  "4" = "#8B1A1A",
  "5" = "#F21A00"
)
dosage_levels <- c("1", "2", "3", "4", "5")

# ── Load shared data ───────────────────────────────────────────────────────────
mds  <- read.table(file.path(base_dir, "mds_matrix.txt"),
                   header = TRUE, row.names = NULL)
indv <- read.table(file.path(base_dir, "sample.ID.txt"),
                   header = FALSE)$V1
cat("Samples loaded:", length(indv), "\n")

# ── Cluster assignments ────────────────────────────────────────────────────────
# Map each cluster list file to a colour/dosage level
cluster_files <- list(
  "1" = file.path(base_dir, "cluster_3.list.txt"),
  "2" = file.path(base_dir, "cluster_1.list.txt"),
  "3" = file.path(base_dir, "cluster_4.list.txt"),
  "4" = file.path(base_dir, "cluster_2.list.txt"),
  "5" = file.path(base_dir, "cluster_5.list.txt")
)

cluster_info <- map_dfr(names(cluster_files), function(cl) {
  ids <- read.table(cluster_files[[cl]], header = FALSE)$V1
  data.frame(id = ids, geno = cl, stringsAsFactors = FALSE)
}) %>%
  mutate(geno = factor(geno, levels = dosage_levels))

cat("Cluster sizes:\n"); print(table(cluster_info$geno))

# ==============================================================================
# PLOT A — Genome-wide MDS (axis 14)
# ==============================================================================

p1 <- mds %>%
  mutate(status = if_else(chrom == 18 & mds14 >= 0.15,
                          "outlier", "non-outlier")) %>%
  ggplot(aes(x = midpos / 1e6, y = mds14)) +
  geom_point(aes(colour = status), size = 1.2, alpha = 0.6) +
  scale_colour_manual(name   = "",
                      values = c("non-outlier" = "grey40",
                                 "outlier"     = "#D1495B")) +
  facet_grid(~ chrom, scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        plot.tag     = element_text(face = "bold")) +
  labs(x = "Mbp", y = "MDS14", tag = "A")

rm(mds); gc()

# ==============================================================================
# PLOT B — Local PCA
# ==============================================================================

pca_coords_df <- read.table(file.path(base_dir,
                                      paste0(LGC, "_pca_coords.txt")),
                            header = TRUE)
pc1_var <- 6.2
pc2_var <- 2.0

pca_df <- pca_coords_df %>%
  left_join(cluster_info, by = "id") %>%
  mutate(geno = factor(geno, levels = dosage_levels))

p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = geno)) +
  geom_point(shape = 21, alpha = 0.5, stroke = 0.1,
             size = 3, colour = "black") +
  scale_fill_manual(values = col_palette, name = "Cluster") +
  xlab(paste0("PC1 (", pc1_var, "%)")) +
  ylab(paste0("PC2 (", pc2_var, "%)")) +
  theme_classic() +
  theme(legend.position = "none",
        plot.tag        = element_text(face = "bold")) +
  coord_fixed(ratio = 0.5) +
  labs(tag = "B")

gc()

# ==============================================================================
# PLOT C — Heterozygosity
# ==============================================================================

het <- read.table(file.path(base_dir, paste0(LGC, ".het")),
                  header = TRUE) %>%
  mutate(O_HET = N_SITES - `O.HOM.`,
         p_HET = O_HET / N_SITES,
         id    = indv) %>%
  select(id, p_HET)

het_df <- cluster_info %>%
  left_join(het, by = "id") %>%
  mutate(geno = factor(geno, levels = dosage_levels))

p3 <- ggplot(het_df, aes(x = geno, y = p_HET, fill = geno)) +
  ggdist::stat_dots(
    side       = "left",
    shape      = 21,
    alpha      = 0.5,
    dotsize    = 0.6,        # smaller dots
    stackratio = 0.9,        # tighter stacking
    binwidth   = 0.003,      # bin width in data units — controls column width
    overflow   = "compress"  # keeps dots within slab boundary
  ) +
  geom_boxplot(width = 0.12, outlier.color = NA,
               lwd = 0.3, colour = "black", alpha = 0.8) +
  theme_classic() +
  scale_fill_manual(values = col_palette, name = "Cluster") +
  labs(x = "", y = "Heterozygosity", tag = "C") +
  theme(axis.text.x  = element_blank(),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        plot.tag     = element_text(face = "bold"))

gc()

# ==============================================================================
# Assemble and save
# ==============================================================================

final <- p1 / (p2 | p3) +
  plot_layout(heights = c(1, 2, 2)) +
  plot_annotation(
    title = "LGC18.01-02",
    theme = theme(plot.title = element_text(size = 18, face = "bold")))

ggsave(file.path(base_dir, paste0(LGC, ".png")),
       plot = final, dpi = 600,
       width = 24, height = 28, units = "cm")

ggsave(file.path(base_dir, paste0(LGC, ".pdf")),
       plot = final, dpi = 600,
       width = 24, height = 28, units = "cm")

message("==> Saved: ", file.path(base_dir, LGC))

EOF
