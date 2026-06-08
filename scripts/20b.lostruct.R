# ==============================================================================
# Script:       20b.lostruct.R
# Description:  Genome-wide windowed local PCA and MDS outlier detection using lostruct.
#
#               Pipeline:
#                 1. Compute local PCA in non-overlapping 1000-SNP windows
#                 2. Compute pairwise distances between windows
#                 3. MDS (20 axes) on pairwise distances
#                 4. Detect outlier windows per MDS axis using SD-based cutoffs;
#                    candidate regions require >= 3 adjacent outlier windows
#                    (gap <= 1 window tolerated)
#                 5. Manhattan plots of all 20 MDS axes with outliers highlighted
#                 6. Correlate local PCA with global PCA to identify driver windows
#
# Input:     Pmax_160_merged_sorted_chr.bcf (indexed BCF from 20a.prepare_lostruct.sh)
#                 Pmax_160_merged_sorted_chr.012 (genotype matrix from 20a.prepare_lostruct.sh)
#
# Output: pca_matrix.txt                    — per-window PCA scores and positions
#                mds_matrix.txt                    — per-window MDS scores and positions
#                mds_outliers.txt                  — outlier windows passing adjacency filter
#                Pmax_160_mds_manhattan.png        — Manhattan plots for all 20 MDS axes
#                Pmax_160_localPC[1/2]_vs_globalPC[1/2].png — local vs global PCA correlation
#
# ==============================================================================

library(data.table)
library(lostruct)
library(ggplot2)
library(patchwork)
library(tidyverse)

options(datatable.fread.input.cmd.message = FALSE)

# ── Paths ──────────────────────────────────────────────────────────────────────
BCF_FILE  <- "/simone/pmax2023/out/20.lostruct/Pmax_160_merged_sorted_chr.bcf"
GENO_FILE <- "/simone/pmax2023/out/20.lostruct/Pmax_160_merged_sorted_chr.012"
OUT_DIR   <- "/simone/pmax2023/out/20.lostruct/"
FIG_DIR   <- "/simone/pmax2023/out/20.lostruct/figures/"

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Chromosome name → number mapping ──────────────────────────────────────────
# Maps scaffold IDs from the P. maximus assembly to two-digit chromosome numbers
chr_map <- c(
  "LR736838.1" = "01", "LR736839.1" = "02", "LR736840.1" = "03",
  "LR736841.1" = "04", "LR736842.1" = "05", "LR736843.1" = "06",
  "LR736844.1" = "07", "LR736845.1" = "08", "LR736846.1" = "09",
  "LR736847.1" = "10", "LR736848.1" = "11", "LR736849.1" = "12",
  "LR736850.1" = "13", "LR736851.1" = "14", "LR736852.1" = "15",
  "LR736853.1" = "16", "LR736854.1" = "17", "LR736855.1" = "18",
  "LR736856.1" = "19"
)

# ==============================================================================
# PART 1: Windowed local PCA
# ==============================================================================
# Divide the genome into non-overlapping 1000-SNP windows and compute the top
# 2 PCs per window using lostruct::eigen_windows.

snps <- vcf_windower(BCF_FILE, size = 1000, type = "snp", sites = vcf_positions(BCF_FILE))

pcs <- eigen_windows(snps, k = 2)
cat("PCA matrix dimensions:", dim(pcs), "\n")

# Check dimension. We have a matrix with 323 columns (3 columns of info, 160 columns with
# PC1 score for each individual, and 160 column with PC2 score for each individual).
# It has as many rows as #windows (3012 with windows of 1000 SNPs)


# Remove windows with NA values
if (anyNA(pcs)) {
  pcs <- pcs[complete.cases(pcs), ]
  warning("Removed windows with NA values from PCA results.")
}

window_pos <- region(snps)()

if (anyNA(window_pos)) {
  window_pos <- window_pos[complete.cases(window_pos), ]
  warning("Removed windows with NA values from window positions.")
}

pca_matrix <- cbind(window_pos, pcs)
n_windows  <- nrow(pca_matrix)
cat("Total windows:", n_windows, "\n")

write.table(pca_matrix, file.path(OUT_DIR, "pca_matrix.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ==============================================================================
# PART 2: Pairwise window distances and MDS
# ==============================================================================
# Compute pairwise distances between window PCA summaries, then embed in 20D via classical MDS (cmdscale).

pcdist  <- pc_dist(pcs, npc = 2)
mds_axe <- cmdscale(pcdist, k = 20)
cat("MDS matrix dimensions:", dim(mds_axe), "\n")
cat("Any NA in MDS matrix:",  anyNA(mds_axe), "\n")

# Build labelled MDS matrix with window mid-points and numeric chromosome codes
mds_matrix        <- as.data.frame(cbind(window_pos, mds_axe))
mds_matrix$midpos <- (mds_matrix$start + mds_matrix$end) / 2
colnames(mds_matrix) <- c("chrom", "start", "end", paste0("mds", 1:20), "midpos")
mds_matrix$chrom     <- chr_map[mds_matrix$chrom]

# Add sequential window index for adjacency filtering
mds_matrix <- mds_matrix %>%
  mutate(N = row_number()) %>%
  select(N, everything())

write.table(mds_matrix, file.path(OUT_DIR, "mds_matrix.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ==============================================================================
# PART 3: Outlier detection and Manhattan plots
# ==============================================================================
# Outliers are defined per MDS axis as windows with scores > cutoff * SD from
# zero (both positive and negative tails). A candidate inversion region requires
# at least 3 adjacent outlier windows on the same axis (gap of <= 1 window
# between consecutive outliers is tolerated).

# ── Outlier detection function ─────────────────────────────────────────────────
# Flags windows exceeding +/- (cutoff * SD) on each MDS axis.
# Returns a data frame of all outlier windows with their axis and direction.

outlier_detect <- function(data, mds_columns, high_cutoff = 3, low_cutoff = 3) {
  outlier_results <- data.frame()
  for (col in mds_columns) {
    if (!col %in% colnames(data)) stop(paste("Column", col, "not found."))
    sd_mds <- sd(data[[col]], na.rm = TRUE)
    pos_outliers <- data %>%
      filter(!!sym(col) >  sd_mds * high_cutoff) %>%
      mutate(outlier = "Outlier", mds_coord = paste0(col, "-pos"))
    neg_outliers <- data %>%
      filter(!!sym(col) < -sd_mds * low_cutoff) %>%
      mutate(outlier = "Outlier", mds_coord = paste0(col, "-neg"))
    outlier_results <- bind_rows(outlier_results, pos_outliers, neg_outliers)
  }
  return(outlier_results)
}

# Apply outlier detection across all 20 MDS axes
outlier_columns <- paste0("mds", 1:20)
mds_outliers    <- outlier_detect(as.data.frame(mds_matrix), outlier_columns)

# ── Adjacency filter ───────────────────────────────────────────────────────────
# Within each chromosome × MDS axis direction, identify runs of consecutive
# outlier windows. A gap of <= 1 window between outliers is tolerated.

mds_outliers <- mds_outliers %>%
  arrange(chrom, mds_coord, N) %>%
  group_by(chrom, mds_coord) %>%
  mutate(gap = N - lag(N), consecutive_group = cumsum(ifelse(is.na(gap) | gap > 2, 1, 0))) %>%
  ungroup() %>%
  arrange(N)

# Retain only runs of >= 3 adjacent outlier windows
outlier_area <- mds_outliers %>%
  group_by(chrom, mds_coord, consecutive_group) %>%
  filter(n() >= 3) %>%
  ungroup()

cat("Candidate regions identified:", n_distinct(outlier_area$consecutive_group), "\n")
cat("Chromosomes affected:", paste(sort(unique(outlier_area$chrom)), collapse = ", "), "\n")

write.table(outlier_area, file.path(OUT_DIR, "mds_outliers.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── Manhattan plots ────────────────────────────────────────────────────────────
# One panel per MDS axis; outlier windows coloured red, non-outliers black.

plot_list <- lapply(1:20, function(i) {

  the_mds           <- paste0("mds", i)
  specific_outliers <- outlier_area %>%
    filter(mds_coord %in% c(paste0(the_mds, "-pos"), paste0(the_mds, "-neg")))

  plot_data <- mds_matrix %>%
    mutate(outlier = ifelse(
      paste(chrom, start, end, midpos) %in%
        paste(specific_outliers$chrom, specific_outliers$start,
              specific_outliers$end,   specific_outliers$midpos),
      "Outlier", "Non-outlier"))

  ggplot(plot_data, aes(x = midpos, y = .data[[the_mds]])) +
    geom_point(size = 2, alpha = 0.5, shape = 21,
               aes(fill = outlier, col = outlier)) +
    scale_colour_manual(name = "",
                        values = c("Non-outlier" = "black", "Outlier" = "red")) +
    scale_fill_manual(name = "",
                      values = c("Non-outlier" = "black", "Outlier" = "red")) +
    facet_grid(cols = vars(chrom), scales = "free_x", space = "free_x") +
    labs(x = "midpos", y = the_mds) +
    theme_classic() +
    theme(axis.text.x     = element_blank(),
          axis.ticks.x    = element_blank(),
          axis.text.y     = element_blank(),
          legend.position = "none")
})

combined_plot <- wrap_plots(plot_list, ncol = 3) + plot_layout(axes = "collect")

ggsave(file.path(FIG_DIR, "Pmax_160_mds_manhattan.png"),
       combined_plot, width = 50, height = 25, units = "cm", dpi = 320)

# ==============================================================================
# PART 4: Correlation of local PCA with global PCA
# ==============================================================================
# Identifies which windows drive the global population structure signal by
# correlating per-window local PC1/PC2 scores with global PC1 and PC2.

# Load genotype matrix (first column is the row index from vcftools --012, drop it)
geno       <- fread(GENO_FILE)[, -1]
Nind       <- nrow(geno)
global_pca <- prcomp(geno)
cat("Individuals:", Nind, "\n")

# pca_matrix column layout:
#   1:3                     = chrom / start / end
#   4:6                     = lostruct summary stats (total, lam_1, lam_2)
#   7:(Nind+6)              = local PC1 scores per individual
#   (Nind+7):(2*Nind+6)     = local PC2 scores per individual
pca_mat <- as.matrix(pca_matrix)

for (global_pc in 1:2) {

  global_vec <- global_pca$x[, global_pc]

  for (local_pc in 1:2) {

    col_start <- if (local_pc == 1) 7 else Nind + 7
    label     <- paste0("localPC", local_pc, "_vs_globalPC", global_pc)

    # Compute absolute Pearson correlation between global PC scores and
    # each window's local PC scores across all individuals
    corr_vec <- vapply(seq_len(n_windows), function(i) {
      local_scores <- t(pca_mat[i, col_start:(col_start + Nind - 1)])
      abs(cor(global_vec, local_scores)[1, 1])
    }, numeric(1))

    df <- as.data.frame(pca_matrix[, 1:3]) %>%
      mutate(corr_vector = corr_vec, midpos = (start + end) / 2, chrom = chr_map[chrom])

    write.table(df, file.path(OUT_DIR, paste0(label, ".txt")), sep = "\t", row.names = FALSE, quote = FALSE)

    p <- ggplot(df, aes(x = midpos, y = corr_vector, colour = chrom)) +
      geom_point(size = 0.8, alpha = 0.8) +
      labs(x = "Mbp",
           y = paste0("|r| local PC", local_pc, " vs global PC", global_pc)) +
      facet_grid(cols = vars(chrom), scales = "free_x", space = "free_x") +
      theme_classic() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

    ggsave(file.path(FIG_DIR, paste0("Pmax_160_", label, ".png")),
           p, width = 42, height = 21, units = "cm", dpi = 320)
  }
}

cat("\nDone. Outputs written to:", OUT_DIR, "\n")
cat("\nNext step:\n")
cat("  Review Pmax_160_mds_manhattan.png to confirm candidate regions.\n")
cat("  Candidate regions (>= 3 adjacent outlier windows) are listed in mds_outliers.txt.\n")
cat("  Run 21a.local_pca.sh to extract regions and prepare input for k-means clustering.\n")
