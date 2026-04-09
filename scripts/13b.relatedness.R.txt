# ==============================================================================
# Script:       13b.relatedness.R
# Description:  Classify pairwise relatedness and visualise results using PLINK IBD statistics 
#                         and ngsRelate kinship coefficients.
#
# ==============================================================================

library(data.table)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(patchwork)
library(wesanderson)

# ── Paths ──────────────────────────────────────────────────────────────────────
GENOME_FILE   <- "relatedness/Pmax_168_merged_pruned.genome"
NGSREL_FILE   <- "relatedness/Pmax_168_merged_pruned.res"
POP_INFO_FILE <- "relatedness/pop_info_all.txt"
OUT_DIR       <- "relatedness/figures/"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Colour palette ─────────────────────────────────────────────────────────────
# Six relatedness categories following Manichaikul et al. (2010)
LEVEL_ORDER <- c(
  "Parent-offspring", "Full-sibling", "Second-degree",
  "Third-degree", "Unrelated", "Unclassified")

pal <- c(
  wes_palette("Cavalcanti1")[3],    # Parent-offspring
  wes_palette("Cavalcanti1")[1],    # Full-sibling
  wes_palette("Cavalcanti1")[2],    # Second-degree
  wes_palette("GrandBudapest1")[4], # Third-degree
  wes_palette("GrandBudapest2")[4], # Unrelated
  wes_palette("Chevalier1")[4]      # Unclassified
)

# ── Sample map ─────────────────────────────────────────────────────────────────
# pop_info_all.txt: index (= PLINK FID, BAM list position) | sample_id
pop_info   <- fread(POP_INFO_FILE, header = TRUE)
sample_map <- setNames(pop_info$sample_id, as.character(pop_info$index))

# ── Read PLINK and ngsRelate output ───────────────────────────────────────────
gen    <- fread(GENOME_FILE, header = TRUE)
ngsrel <- fread(NGSREL_FILE, header = TRUE)

# Verify row counts match before merging (one row per pair in both files)
stopifnot("Row count mismatch between .genome and .res" = nrow(gen) == nrow(ngsrel))

gen$R1   <- ngsrel$R1
gen$R0   <- ngsrel$R0
gen$KING <- ngsrel$KING

# ── Kinship classification (Manichaikul et al. 2010) ──────────────────────────
# Kinship coefficient: ϕ = Z1/4 + Z2/2
# Classification thresholds (expressed as powers of 1/2):
#   ϕ ≥ 1/2^(5/2) (~0.177) and < 1/2^(3/2) (~0.354): first-degree relatives
#     Z0 < 0.1             → Parent-offspring
#     0.1 < Z0 < 0.365     → Full-sibling
#   ϕ ≥ 1/2^(7/2) (~0.088): second-degree
#   ϕ ≥ 1/2^(9/2) (~0.044): third-degree
#   ϕ < 1/2^(9/2)          → Unrelated
#   Pairs not meeting any criterion → Unclassified

gen <- gen %>%
  mutate(
    kinship  = (Z1 / 4) + (Z2 / 2),
    criteria = case_when(
      kinship >= 1/2^(5/2) & kinship < 1/2^(3/2) & Z0 < 0.1                                 ~ "Parent-offspring",
      kinship >= 1/2^(5/2) & kinship < 1/2^(3/2) & Z0 > 0.1 & Z0 < 0.365                   ~ "Full-sibling",
      kinship >= 1/2^(7/2) & kinship < 1/2^(5/2) & Z0 > 0.365 & Z0 < 1 - (1/2^(3/2))      ~ "Second-degree",
      kinship >= 1/2^(9/2) & kinship < 1/2^(7/2) & Z0 > 1-(1/2^(3/2)) & Z0 < 1-(1/2^(5/2))~ "Third-degree",
      kinship < 1/2^(9/2)  & Z0 > 1 - (1/2^(5/2))                                           ~ "Unrelated",
      TRUE ~ "Unclassified"
    ),
    criteria  = factor(criteria, levels = LEVEL_ORDER),
    # Resolve sample names via FID index
    IID1_name = sample_map[as.character(FID1)],
    IID2_name = sample_map[as.character(FID2)]
  )

# ── Sanity check: warn if any FID did not resolve to a sample name ─────────────
if (any(is.na(gen$IID1_name) | is.na(gen$IID2_name))) {
  warning("Some FID values did not match pop_info_all indices — check FID1/FID2 columns.")
  print(gen[is.na(IID1_name) | is.na(IID2_name), .(FID1, IID1, FID2, IID2)] |> head(10))}

# ── Console summary ────────────────────────────────────────────────────────────
cat("\n=== PI_HAT summary ===\n");        print(summary(gen$PI_HAT))
cat("\n=== Classification counts ===\n"); print(table(gen$criteria))
cat("\nMean PI_HAT:", round(mean(gen$PI_HAT, na.rm = TRUE), 4), "\n")

# Print unclassified pairs
cat("\n=== Unclassified pairs (possible cross-contamination — flag for removal) ===\n")
unk_pairs <- gen[criteria == "Unclassified",
                 .(IID1_name, IID2_name, PI_HAT, KING, R1, R0, Z0, Z1, Z2, kinship)]
print(unk_pairs)

cat("\n=> Samples involved in unclassified pairs:\n")
unk_samples <- unique(c(unk_pairs$IID1_name, unk_pairs$IID2_name))
cat(paste(unk_samples, collapse = ", "), "\n")
cat("=> Check lab records for cross-contamination; consider removing both samples from each pair.\n")

# Print full-sibling pairs
cat("\n=== Full-sibling pairs ===\n")
print(gen[criteria == "Full-sibling",
          .(IID1_name, IID2_name, PI_HAT, KING, R1, R0, Z0, Z1, Z2, kinship)])

# ── Plot 1: PI_HAT distribution ───────────────────────────────────────────────
p_hist <- ggplot(gen, aes(x = PI_HAT)) +
  geom_histogram(binwidth = 0.03, alpha = 0.9, fill = "grey40", colour = "white") +
  xlab("Pairwise relatedness (PI_HAT)") + ylab("Count") +
  theme_few() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))

p_hist

# ── Plot 2: IBD probability triangle plots ─────────────────────────────────────
# Three pairwise plots of IBD probabilities (Z0, Z1, Z2) 

p_z0z1 <- ggplot(gen, aes(Z0, Z1, colour = criteria)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(values = setNames(pal, LEVEL_ORDER)) +
  xlab("Pr(IBD=0)") + ylab("Pr(IBD=1)") +
  theme_few() + theme(legend.title = element_blank())

p_z0z2 <- ggplot(gen, aes(Z0, Z2, colour = criteria)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(values = setNames(pal, LEVEL_ORDER)) +
  xlab("Pr(IBD=0)") + ylab("Pr(IBD=2)") + ylim(c(0, 1)) +
  theme_few() + theme(legend.title = element_blank())

p_z1z2 <- ggplot(gen, aes(Z1, Z2, colour = criteria)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(values = setNames(pal, LEVEL_ORDER)) +
  xlab("Pr(IBD=1)") + ylab("Pr(IBD=2)") + ylim(c(0, 1)) +
  theme_few() + theme(legend.title = element_blank())

p_z_combined <- p_z0z1 + p_z0z2 + p_z1z2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

p_z_combined

# ── Plot 3: KING vs log10(R1) ─────────────────────────────────────────────────
# KING-robust kinship estimator (y) vs ngsRelate R1 coefficient (x, log scale).
# Dashed line marks the first-degree KING threshold (1/2^(5/2) ≈ 0.177).

KING_THRESH <- 1/2^(5/2)

p_relate <- ggplot(gen, aes(R1, KING)) +
  geom_point(size = 4, alpha = 0.8,
             aes(colour = factor(criteria, levels = LEVEL_ORDER))) +
  scale_colour_manual(values = setNames(pal, LEVEL_ORDER)) +
  scale_x_log10() +
  geom_hline(yintercept = KING_THRESH, linetype = "dashed", alpha = 0.3) +
  xlab(expression(log[10](R1))) +
  ylab("KING") +
  labs(title = "") +
  theme_few() +
  theme(legend.title          = element_blank(),
        legend.text           = element_text(size = 16),
        axis.title            = element_text(size = 18),
        axis.text.x           = element_text(size = 16),
        axis.text.y           = element_text(size = 16),
        plot.title            = element_text(size = 14, face = "bold"),
        legend.position.inside = c(0.2, 0.9)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1),
                               keyheight    = 0.1,
                               default.unit = "inch"))

p_relate

# ── Save figures ───────────────────────────────────────────────────────────────
ggsave(file.path(OUT_DIR, "Pmax_168_PI_HAT_hist.png"),
       p_hist,       width = 6,  height = 5,  dpi = 300)
ggsave(file.path(OUT_DIR, "Pmax_168_Z_scores.png"),
       p_z_combined, width = 14, height = 5,  dpi = 300)
ggsave(file.path(OUT_DIR, "Pmax_168_R1_KING.png"),
       p_relate,     width = 10, height = 7,  dpi = 300)

cat("\nDone. Figures saved to:", OUT_DIR, "\n")
