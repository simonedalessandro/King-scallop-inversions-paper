# ==============================================================================
# Script:       21b.local_pca.R
# Description:  Local PCA and k-means clustering for each candidate inversion
#               region identified from lostruct MDS (20b.lostruct.R).
#
#               For each region:
#                 1. Load --012 genotype matrix
#                 2. Run PCA (prcomp)
#                 3. Assign individuals to chromosomal arrangement clusters
#                    via k-means on PC1 (initialised at min, midpoint, max)
#                 4. Evaluate clustering quality: between-SS / total-SS >= 95%
#                    retains region in the inversion dataset
#                 5. Save cluster assignments (AA / AB / BB lists) and PCA plot
#
# Input:        Per-region --012 genotype matrices from 21a.local_pca.sh
#               sample.ID.txt — one sample ID per line, matching --012 row order
#
# Output:       Per-region cluster assignment tables, AA/AB/BB sample lists,
#               and PCA scatter plots coloured by k-means cluster
# ==============================================================================

library(ggplot2)
library(dplyr)
library(purrr)

# ── Paths ──────────────────────────────────────────────────────────────────────
INV_DIR     <- "/simone/pmax2023/out/20.lostruct/inversions/"
SAMPLE_FILE <- "/simone/pmax2023/out/file_lists/sample.ID.txt"
OUT_DIR     <- "/simone/pmax2023/out/20.lostruct/inversions/local_pca/"
FIG_DIR     <- "/simone/pmax2023/out/20.lostruct/inversions/figures/"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Sample IDs ─────────────────────────────────────────────────────────────────
# One ID per line; order must match the row order of all --012 genotype matrices
indv <- read.table(SAMPLE_FILE, header = FALSE)[, 1]

# ── Candidate regions to process ──────────────────────────────────────────────
regions <- c(
  "LGC02",    "LGC08.01", "LGC08.02", "LGC09",
  "LGC11",    "LGC13.01", "LGC13.02", "LGC13.03",
  "LGC15",    "LGC17.01", "LGC17.02", "LGC18.01",
  "LGC18.02", "LGC18.03", "LGC18.04", "LGC18.01_02",
  "LGC19.01", "LGC19.02", "LGC01",    "LGC04",    "LGC06"
)

# ── k-means clusters to try per region ────────────────────────────────────────
# Default is 3 (homokaryotic non-inverted AA / heterokaryotic AB /
# homokaryotic inverted BB). Adjusted for regions with more complex structure
# visible on local PCA.
k_per_region <- setNames(rep(3, length(regions)), regions)
k_per_region["LGC06"]       <- 6   # 6 distinct clusters on local PCA
k_per_region["LGC18.01_02"] <- 5   # 5 distinct clusters on local PCA

# ==============================================================================
# Loop over candidate regions
# ==============================================================================

results_summary <- list()

for (code in regions) {

  geno_file <- file.path(INV_DIR, paste0(code, ".012"))

  if (!file.exists(geno_file)) {
    warning("File not found, skipping: ", geno_file)
    next
  }

  cat("\n[", code, "]\n", sep = "")

  # ── Load genotype matrix ──────────────────────────────────────────────────
  # vcftools --012 output: first column is the row index — drop it
  geno           <- read.table(geno_file)[, -1]
  rownames(geno) <- indv

  # ── PCA ───────────────────────────────────────────────────────────────────
  pca     <- prcomp(geno)
  var_exp <- pca$sdev^2 / sum(pca$sdev^2) * 100

  # ── k-means clustering on PC1 ─────────────────────────────────────────────
  k       <- k_per_region[code]
  pc1     <- pca$x[, 1]
  centres <- seq(min(pc1), max(pc1), length.out = k)

  set.seed(42)
  km <- kmeans(pc1, centers = centres)

  wss_pct <- round(km$betweenss / km$totss * 100, 1)

  # Retain regions with >= 95% between-cluster variance
  retained <- wss_pct >= 95
  cat("  Retained:", ifelse(retained, "YES", "NO — check clustering manually"), "\n")

  # ── Build cluster assignment table ────────────────────────────────────────
  cluster_df <- data.frame(id = indv, cluster = km$cluster)

  # Label clusters by PC1 centre position (lowest -> AA, highest -> BB).
  # For k <= 3 use the AA/AB/BB karyotype labels; for k > 3 (regions with more
  # complex structure) fall back to ordered numeric labels to avoid NA labels.
  centre_order <- order(km$centers)
  if (k <= 3) {
    base_labels <- if (k == 2) c("AA", "BB") else c("AA", "AB", "BB")[seq_len(k)]
  } else {
    base_labels <- as.character(seq_len(k))
  }
  cluster_labels     <- setNames(base_labels, centre_order)
  cluster_df$genotype <- cluster_labels[as.character(cluster_df$cluster)]

  write.table(cluster_df,
              file.path(OUT_DIR, paste0(code, "_cluster_assignments.txt")),
              quote = FALSE, row.names = FALSE)

  # Write per-genotype sample lists (inputs for downstream Fst / LD analyses)
  for (gt in unique(cluster_df$genotype)) {
    samples <- cluster_df$id[cluster_df$genotype == gt]
    write.table(samples,
                file.path(OUT_DIR, paste0(code, "_", gt, ".list.txt")),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }

  # ── PCA scatter plot coloured by cluster ──────────────────────────────────
  pca_df <- as.data.frame(pca$x[, 1:2]) %>%
    mutate(sample  = indv,
           cluster = factor(cluster_df$genotype))

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = cluster, label = sample)) +
    geom_point(size = 3, alpha = 0.9) +
    xlab(sprintf("PC1 (%.2f%%)", var_exp[1])) +
    ylab(sprintf("PC2 (%.2f%%)", var_exp[2])) +
    labs(title    = code,
         subtitle = sprintf("k = %d | Between-SS/Total-SS = %.1f%%", k, wss_pct),
         colour   = "Cluster") +
    theme_classic() +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 10))

  ggsave(file.path(FIG_DIR, paste0(code, "_local_pca.png")),
         p, width = 7, height = 6, dpi = 300)

  # ── Store summary ─────────────────────────────────────────────────────────
  results_summary[[code]] <- data.frame(
    region   = code,
    k        = k,
    wss_pct  = wss_pct,
    retained = retained,
    pc1_var  = round(var_exp[1], 2),
    pc2_var  = round(var_exp[2], 2)
  )
}

# ── Print and save summary table ──────────────────────────────────────────────
cat("\n==================================================\n")
cat("SUMMARY\n")
cat("==================================================\n")
summary_df <- do.call(rbind, results_summary)
print(summary_df[order(summary_df$wss_pct, decreasing = TRUE), ])

write.table(summary_df,
            file.path(OUT_DIR, "clustering_summary.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")

cat("\nDone. Outputs written to:", OUT_DIR, "\n")
cat("Figures saved to:", FIG_DIR, "\n")
cat("\nNext step:\n")
cat("  Regions with between-SS/total-SS >= 95% are retained in the inversion dataset.\n")
cat("  Cluster assignment files are inputs for downstream Fst / LD analyses.\n")
