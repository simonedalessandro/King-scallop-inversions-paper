#!/bin/sh
# ==============================================================================
# Script:       24.MDS_pca_het_LD_plot.sh
# Description:  Generate multi-panel figures (MDS, local PCA, heterozygosity,
#               LD heatmap) for each putative inversion.
#
#               Panel layout per inversion:
#                 A — genome-wide MDS with highlighted inversion window
#                 B — local PCA coloured by k-means cluster
#                 C — per-cluster heterozygosity
#                 D — chromosome-wide LD heatmap (all vs AA homozygotes)
#
# Input:        From 20b.lostruct.R:
#                 mds_matrix.txt     — chromosome, midpoint, MDS values
#                 inversion_info.txt — LGC ID, chromosome, MDS axis, thresholds
#               From 21b.local_pca.R:
#                 <LGC>.cov          — local covariance matrix
#                 <LGC>_geno.txt     — cluster calls per inversion
#                 pop_info_rm8.txt   — sample ID and location
#               From 22.het.sh:
#                 <LGC>.het          — per-individual heterozygosity
#               From 23.LD.sh:
#                 <LGC>.ld           — LD for all individuals
#                 <LGC>.AA.ld        — LD for common-homozygous (AA) individuals
#
# Output:       Per-inversion figures: <LGC>.png and <LGC>.pdf
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE job configuration
# ------------------------------------------------------------------------------
#$ -N MDS_plot
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=30G
#$ -pe sharedmem 4
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
library(plotly)
library(patchwork)
library(tidyverse)
library(wesanderson)
library(gghalves)

# ── Colour palette ─────────────────────────────────────────────────────────────
# Three colours for AA / AB / BB chromosomal arrangement clusters
col_palette <- c(
  wes_palette("Zissou1")[1],
  wes_palette("Zissou1")[3],
  wes_palette("Zissou1")[5]
)

# ── Load MDS matrix and inversion metadata ─────────────────────────────────────
mds  <- read.table("/simone/pmax2023/out/21.MDS_pca_het_LD/mds_matrix.txt", header = TRUE)
info <- read.table("/simone/pmax2023/out/21.MDS_pca_het_LD/inversion_info.txt", header = TRUE)

# ── Load sample metadata ───────────────────────────────────────────────────────
pop_info <- read.table("/simone/pmax2023/out/21.MDS_pca_het_LD/local_pca/pop_info_rm8.txt", header = TRUE) %>%
  rename(id = Sample, location = Location) %>%
  select(id, location)

# ==============================================================================
# PLOT A — Genome-wide MDS with highlighted inversion window
# ==============================================================================
# Points within the inversion's MDS thresholds on the relevant chromosome are
# coloured red.

mds_plot <- function(mds_df, inv_row) {

  axis_col <- paste0("mds", inv_row$MDS)
  chr      <- inv_row$chrom
  max_thr  <- inv_row$`if.MDS.max`
  min_thr  <- inv_row$`if.MDS.min`

  mds_df %>%
    mutate(
      mds_value = .data[[axis_col]],
      status    = if_else(
        chrom == chr & mds_value <= max_thr & mds_value >= min_thr,
        "outlier", "non-outlier")
    ) %>%
    ggplot(aes(x = midpos / 1e6, y = mds_value)) +
    geom_point(aes(colour = status), size = 1.2, alpha = 0.6) +
    scale_colour_manual(
      values = c("non-outlier" = "grey40", "outlier" = "red"),
      name   = "") +
    facet_grid(~chrom, scales = "free_x", space = "free_x") +
    theme_classic() +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.tag     = element_text(face = "bold")) +
    labs(x   = "Mbp",
         y   = paste0("MDS", inv_row$MDS),
         tag = "A")
}

plots_mds        <- purrr::map(seq_len(nrow(info)), ~ mds_plot(mds, info[.x, ]))
names(plots_mds) <- info$LGC

# ==============================================================================
# PLOT B — Local PCA per inversion
# ==============================================================================
# PCA computed from the covariance matrix output by PCAngsd for the inversion
# window; points coloured by k-means cluster assignment.

local_pca <- function(LGC) {

  cov_path <- file.path(
    "/simone/pmax2023/out/21.MDS_pca_het_LD/local_pca",
    paste0(LGC, ".cov"))

  cov <- as.matrix(read.table(cov_path))
  e   <- eigen(cov)

  total_var <- sum(e$values)
  var_expl  <- (e$values[1:2] / total_var) * 100
  names(var_expl) <- c("PC1", "PC2")

  pcs           <- as.data.frame(e$vectors[, 1:2])
  colnames(pcs) <- c("PC1", "PC2")

  # Load k-means cluster assignments from 21b.local_pca.R
  geno_path    <- file.path(
    "/simone/pmax2023/out/21.MDS_pca_het_LD/local_pca",
    paste0(LGC, "_geno.txt"))
  cluster_info <- read.table(geno_path, header = TRUE) %>%
    rename(id = id_inv, geno = cluster_inv)

  pca_pop <- pcs %>%
    mutate(id = pop_info$id) %>%
    left_join(pop_info,     by = "id") %>%
    left_join(cluster_info, by = "id")

  ggplot(pca_pop, aes(x = PC1, y = PC2, fill = factor(geno))) +
    geom_point(shape = 21, alpha = 0.5, size = 4,
               colour = "black", stroke = 0.1) +
    scale_fill_manual(values = col_palette, name = "geno") +
    xlab(paste0("PC1 (", round(var_expl["PC1"], 1), "%)")) +
    ylab(paste0("PC2 (", round(var_expl["PC2"], 1), "%)")) +
    coord_fixed(0.5) +
    theme_classic() +
    theme(legend.position = "none",
          plot.tag        = element_text(face = "bold")) +
    labs(tag = "B")
}

pca_plots        <- purrr::map(info$LGC, local_pca)
names(pca_plots) <- info$LGC

# ==============================================================================
# PLOT C — Heterozygosity per genotype cluster
# ==============================================================================
# Half-point / half-boxplot showing the distribution of observed heterozygosity
# per chromosomal arrangement cluster for each inversion region.

het_plot <- function(LGC) {

  # Load cluster assignments
  geno_path    <- file.path(
    "/simone/pmax2023/out/21.MDS_pca_het_LD/local_pca",
    paste0(LGC, "_geno.txt"))
  cluster_info <- read.table(geno_path, header = TRUE) %>%
    rename(id = id_inv, geno = cluster_inv)

  # Load heterozygosity output from 22.het.sh and compute p_HET
  het_path <- file.path(
    "/simone/pmax2023/out/21.MDS_pca_het_LD/het",
    paste0(LGC, ".het"))
  het <- read.table(het_path, header = TRUE) %>%
    mutate(O_HET = N_SITES - `O.HOM.`,
           p_HET = O_HET / N_SITES)

  pop_loc       <- pop_info
  pop_loc$p_HET <- het$p_HET

  LGC_info <- pop_loc %>%
    left_join(cluster_info, by = "id")

  ggplot(LGC_info, aes(x = factor(geno), y = p_HET, fill = factor(geno))) +
    gghalves::geom_half_point(side = "l", shape = 21, alpha = 0.5,
                              size = 4, stroke = 0.1) +
    gghalves::geom_half_boxplot(side = "r", outlier.color = NA,
                                width = 0.3, lwd = 0.3,
                                colour = "black", alpha = 0.8) +
    scale_fill_manual(values = col_palette) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          plot.tag    = element_text(face = "bold")) +
    coord_fixed(15) +
    labs(y = "Heterozygosity", x = "", tag = "C")
}

het_plots        <- purrr::map(info$LGC, het_plot)
names(het_plots) <- info$LGC

# ==============================================================================
# PLOT D — LD heatmap (all individuals vs AA homozygotes)
# ==============================================================================
# Upper triangle: LD for all individuals; lower triangle: LD for AA homozygotes.
# Overlaying the two datasets allows visual comparison of LD block structure
# between the full cohort and the non-inverted arrangement.

ld_plot <- function(LGC) {

  # Load LD for all individuals
  ld_path <- file.path(
    "/simone/pmax2023/out/21.MDS_pca_het_LD/LD",
    paste0(LGC, ".ld"))
  ld <- read.table(ld_path, header = TRUE)

  df <- data.frame(BP_A = ld$BP_A, BP_B = ld$BP_B, R2 = ld$R2) %>%
    mutate(win1 = as.factor(BP_A / 1e6),
           win2 = as.factor(BP_B / 1e6)) %>%
    select(-BP_A, -BP_B)

  # Load LD for AA homozygous cluster
  ldAA_path <- file.path(
    "/simone/pmax2023/out/21.MDS_pca_het_LD/LD",
    paste0(LGC, ".AA.ld"))
  ld_AA <- read.table(ldAA_path, header = TRUE)

  df_AA <- data.frame(BP_A = ld_AA$BP_A, BP_B = ld_AA$BP_B, R2 = ld_AA$R2) %>%
    mutate(win1 = as.factor(BP_A / 1e6),
           win2 = as.factor(BP_B / 1e6)) %>%
    select(-BP_A, -BP_B)

  # Upper triangle = all individuals; lower triangle = AA homozygotes
  ggplot(df, aes(x = win1, y = win2, fill = R2)) +
    geom_tile() +
    geom_tile(data = df_AA, aes(x = win2, y = win1)) +
    scale_fill_gradientn(
      colours = c("grey95", "blue", "red"),
      values  = c(0, 0.5, 1),
      name    = "R2") +
    xlab("Mbp") + ylab("Mbp") +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks  = element_blank(),
          plot.tag    = element_text(face = "bold")) +
    coord_fixed(ratio = 1) +
    labs(tag = "D")
}

ld_plots        <- purrr::map(info$LGC, ld_plot)
names(ld_plots) <- info$LGC

# ==============================================================================
# Assemble and save multi-panel figures
# ==============================================================================
# Layout: A on top (full width); B | C | D on bottom row

save_inversion_fig <- function(LGC) {

  final <- plots_mds[[LGC]] / (pca_plots[[LGC]] | het_plots[[LGC]] | ld_plots[[LGC]]) +
    plot_annotation(
      title = LGC,
      theme = theme(plot.title = element_text(size = 18)))

  fig_dir  <- "/simone/pmax2023/out/21.MDS_pca_het_LD/figures"
  png_path <- file.path(fig_dir, paste0(LGC, ".png"))
  pdf_path <- file.path(fig_dir, paste0(LGC, ".pdf"))

  ggsave(filename = png_path, plot = final, dpi = "retina",
         width = 16, height = 8, units = "cm")
  ggsave(filename = pdf_path, plot = final, dpi = "retina",
         width = 16, height = 8, units = "cm")
}

# Generate figures for all inversions
walk(info$LGC, save_inversion_fig)

EOF
