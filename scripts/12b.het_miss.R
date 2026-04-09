# ==============================================================================
# Script:       12b.het_miss.R
# Description:  Individual-level QC for the P. maximus 168-sample dataset.
#               Computes observed heterozygosity (p_HET) per individual and
#               plots its distribution alongside per-individual missingness
# ==============================================================================

library(ggplot2)
library(dplyr)

# ── Paths ──────────────────────────────────────────────────────────────────────
HET_FILE  <- "/het&miss/Pmax_168_merged_chr.het"
MISS_FILE <- "/het&miss/Pmax_168_merged_chr.imiss"
POP_FILE  <- "/het&miss/pop_info_all.txt"
OUT_DIR   <- "/het&miss/figures/"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Sample map ─────────────────────────────────────────────────────────────────
# pop_info_all.txt has two columns: index (= PLINK FID) and sample_id
pop_info   <- read.table(POP_FILE, header = TRUE)
sample_map <- setNames(pop_info$sample_id, as.character(pop_info$index))

# ── Heterozygosity ─────────────────────────────────────────────────────────────
het <- read.table(HET_FILE, header = TRUE)

# Map FID index to sample name; compute observed heterozygosity proportion:
#   O_HET  = non-missing sites − observed homozygotes
#   p_HET  = O_HET / total non-missing sites (N.NM.)
het <- het %>%
  mutate(sample_id = sample_map[as.character(FID)],
         O_HET     = N.NM. - O.HOM.,
         p_HET     = O_HET / N.NM.) %>%
  rename(ID = IID)

# Flag individuals > 3 SD from the mean p_HET
mu        <- mean(het$p_HET, na.rm = TRUE)
sd_val    <- sd(het$p_HET,   na.rm = TRUE)
lower_thr <- mu - 3 * sd_val
upper_thr <- mu + 3 * sd_val

het <- het %>%
  mutate(is_outlier = p_HET < lower_thr | p_HET > upper_thr)

# Print outlier table and threshold summary to console
cat("\n=== Heterozygosity outliers (> 3 SD from mean) ===\n")
print(het[het$is_outlier, c("sample_id", "O.HOM.", "O_HET", "N.NM.", "p_HET")])
cat(sprintf("Mean p_HET: %.4f | SD: %.4f | Thresholds: [%.4f, %.4f]\n",
            mu, sd_val, lower_thr, upper_thr))

# Plot heterozygosity distribution
# Red rug ticks mark outlier individuals; dashed red lines show ± 3 SD thresholds
het_plot <- ggplot(het, aes(x = p_HET)) +
  geom_histogram(binwidth = 0.005, fill = "#FFCC33",
                 colour = "#996600", alpha = 0.85) +
  geom_rug(data = subset(het, is_outlier), aes(x = p_HET),
           colour = "red", linewidth = 1, inherit.aes = FALSE) +
  geom_vline(xintercept = c(lower_thr, upper_thr),
             linetype = "dashed", colour = "red", alpha = 0.5) +
  labs(title    = "Heterozygosity frequency distribution",
       subtitle = "Red ticks/lines mark individuals > 3 SD from mean",
       x        = "Observed heterozygosity (p_HET)",
       y        = "Count of individuals") +
  theme_minimal()

print(het_plot)
ggsave(file.path(OUT_DIR, "Pmax_168_het.png"), het_plot, width = 7, height = 5, dpi = 300)

# ── Missingness ────────────────────────────────────────────────────────────────
miss <- read.table(MISS_FILE, header = TRUE)

# Map FID index to sample name
miss <- miss %>%
  mutate(sample_id = sample_map[as.character(FID)])

# Flag individuals > 3 SD from the mean F_MISS
mu_miss        <- mean(miss$F_MISS, na.rm = TRUE)
sd_miss        <- sd(miss$F_MISS,   na.rm = TRUE)
lower_thr_miss <- mu_miss - 3 * sd_miss
upper_thr_miss <- mu_miss + 3 * sd_miss

miss <- miss %>%
  mutate(is_outlier = F_MISS < lower_thr_miss | F_MISS > upper_thr_miss)

# Print outlier table to console
cat("\n=== Missingness outliers (> 3 SD from mean) ===\n")
print(miss[miss$is_outlier, c("sample_id", "N_MISS", "N_GENO", "F_MISS")])

# Plot missingness distribution
# Red rug ticks mark outlier individuals; dashed red lines show ± 3 SD thresholds
miss_plot <- ggplot(miss, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.001, fill = "#FFCC33",
                 colour = "#CC9900", alpha = 0.8, boundary = 0) +
  geom_rug(data = subset(miss, is_outlier), aes(x = F_MISS),
           sides = "b", colour = "red", linewidth = 1, inherit.aes = FALSE) +
  geom_vline(xintercept = c(lower_thr_miss, upper_thr_miss),
             linetype = "dashed", colour = "red", alpha = 0.5) +
  labs(title    = "Missingness frequency distribution",
       subtitle = "Red ticks/lines mark individuals > 3 SD from mean",
       x        = "Missing genotype frequency (F_MISS)",
       y        = "Count of individuals") +
  theme_minimal()

print(miss_plot)
ggsave(file.path(OUT_DIR, "Pmax_168_miss.png"), miss_plot, width = 7, height = 5, dpi = 300)

cat("\nDone. Figures saved to:", OUT_DIR, "\n")
