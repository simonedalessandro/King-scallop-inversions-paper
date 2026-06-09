#!/bin/bash
# ==============================================================================
# Script:       32_map_het_DAPC_IBD.sh
# Description:  Calculate and visualise population genetic diversity statistics
#               for the inversion-free 158-sample dataset, and assemble a
#               four-panel composite figure.
#
#               Pipeline:
#                 0. Generate per-population heterozygosity with VCFtools --het
#                 A. Sampling map
#                 B. Pairwise FST + observed heterozygosity
#                 C. DAPC: LD1 vs LD2 scatter with 95 % confidence ellipses
#                 D. Isolation by Distance: Mantel test + regression plot
#
# Input:        Pecmax_no_inv_pruned_filtered.vcf  — from 31_Prune_filter.sh
#               ID_pop.txt                         — sample ID / population assignment, tab-separated, header row: ID  Location
#               locations_coordinates.xlsx         — sampling site coordinates
#               *.shp                              — Scottish marine regions
#               DistvsFst.txt                      — MANUAL input for Panel D, tab-separated, columns: Comparison  Fst  Dist  Labels
#                                                    Built by hand: the Fst column comes from a separate pairwise dartR
#                                                    gl.fst.pop run, paired with minimal marine distances (Dist, km)
#                                                    computed following Assis et al. (2013).
#
# Output:       Figure_3.png / .pdf
#               fst_pop.txt / fst_pvalue_pop.txt
#               <POP>.het                           — per-population het, one per population
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE job directives
# ------------------------------------------------------------------------------
#$ -N 32.map_het_dapc_ibd        
#$ -cwd                                             
#$ -o map_het_dapc_ibd.log    
#$ -l h_rt=04:00:00           
#$ -l h_vmem=8G               
#$ -pe sharedmem 6            

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

BASE_DIR=/simone/pmax2023/out/23.het_DAPC_IBD
MAP_DIR=${BASE_DIR}/map
HET_DIR=${BASE_DIR}/het
DAPC_DIR=${BASE_DIR}/DAPC
IBD_DIR=${BASE_DIR}/IBD
FIG_DIR=${BASE_DIR}/figures

VCF_FILE=${BASE_DIR}/Pecmax_no_inv_pruned_filtered.vcf
ID_POP=${BASE_DIR}/ID_pop.txt

COORDS_FILE=${MAP_DIR}/locations_coordinates.xlsx
COAST_FILE=${MAP_DIR}/area_management_scallop_assessment_areas/area_management_scallop_assessment_areasPolygon.shp
DISTVSFST_FILE=${IBD_DIR}/DistvsFst.txt

NCORES=${NSLOTS:-6}

mkdir -p "${MAP_DIR}" "${HET_DIR}" "${DAPC_DIR}" "${IBD_DIR}" "${FIG_DIR}"

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load vcftools/0.1.16
module load R/4.5

# ------------------------------------------------------------------------------
# Generate per-population heterozygosity (VCFtools --het)
# ------------------------------------------------------------------------------
echo "[$(date)] === Generating per-population .het files ==="

# Clear any stale per-population files from previous runs
rm -f "${HET_DIR}"/*.samples "${HET_DIR}"/*.het "${HET_DIR}"/*.het.log

# Write one <POP>.samples keep-list per population
awk -F'\t' -v outdir="${HET_DIR}" '
    NR > 1 {
        loc = $2
        if (loc == "MFR") loc = "MOR"
        print $1 > (outdir "/" loc ".samples")
    }
' "${ID_POP}"

# Run VCFtools --het for each population keep-list
for samplelist in "${HET_DIR}"/*.samples; do
    pop=$(basename "${samplelist}" .samples)
    echo "[$(date)] VCFtools --het: ${pop}"
    vcftools \
        --vcf  "${VCF_FILE}" \
        --keep "${samplelist}" \
        --het \
        --out  "${HET_DIR}/${pop}"
done

echo "[$(date)] Per-population .het files written to ${HET_DIR}"

# ------------------------------------------------------------------------------
# Run R
# ------------------------------------------------------------------------------

Rscript - \
  --args \
  "${BASE_DIR}" \
  "${NCORES}" \
  "${VCF_FILE}" \
  "${ID_POP}" \
  "${HET_DIR}" \
  "${COORDS_FILE}" \
  "${COAST_FILE}" \
  "${DISTVSFST_FILE}" \
  "${FIG_DIR}" \
<<'EOF'

# -- Load packages -------------------------------------------------------------
library(tidyverse)
library(vcfR)
library(adegenet)
library(hierfstat)
library(adespatial)
library(vegan)
library(gghalves)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggspatial)
library(readxl)
library(data.table)
library(cowplot)
library(dartR)

# ── Colour palette and legend ─────────────────────────────────────────────────
col_palette <- c("SHE" = "#0477BF", "ORK" = "#38E0F2", "MNC" = "#29A655",
                 "MOR" = "#7FBF50", "FRA" = "#F28891", "CLY" = "#E2D200",
                 "TAR" = "#FFA75C", "CDB" = "#DC143C")

legend_order <- c("SHE", "ORK", "MNC", "MOR", "FRA", "CLY", "TAR", "CDB")

legend_labels <- c("SHE" = "Shetland", "ORK" = "Orkney", "MNC" = "The Minch",
                   "MOR" = "Moray Firth", "FRA" = "Fraserburgh", "CLY" = "Clyde",
                   "TAR" = "Targets", "CDB" = "Cardigan Bay")

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
# VCF + ID_pop.txt
# ==============================================================================

vcf <- read.vcfR(VCF_FILE)

gen_pops <- read_delim(
  ID_POP,
  delim         = "\t",
  escape_double = FALSE,
  col_names     = TRUE,
  trim_ws       = TRUE) %>%
  mutate(Location = case_when(Location == "MFR" ~ "MOR", TRUE ~ Location))

# Build genlight once and assign population labels
genotypes      <- vcfR2genlight(vcf)
genotypes@pop  <- as.factor(gen_pops$Location)

# Sanity check: VCF sample order must match ID_pop.txt (assignment is positional)
if (!all(indNames(genotypes) == gen_pops$ID)) {
  warning("Sample order in VCF does not match ID_pop.txt — populations may be ",
          "mis-assigned. Ensure ID_pop.txt is in the same order as the VCF samples.")
}

# ==============================================================================
# PANEL A - Sampling map
# ==============================================================================

site_coords <- readxl::read_xlsx(COORDS_FILE) %>%
  mutate(ID = factor(ID, levels = legend_order))

sf_use_s2(FALSE)

uk_map <- ne_countries(scale = "large", returnclass = "sf") %>%
  st_make_valid() %>%
  st_crop(xmin = -11, xmax = 3, ymin = 49.5, ymax = 61.5)

# Scottish marine regions
coast <- st_read(COAST_FILE, quiet = TRUE) %>%
  st_make_valid()
if (!is.na(st_crs(coast))) coast <- st_transform(coast, 4326)

pA <- ggplot() +
  geom_sf(data = coast, fill = NA, colour = "grey30", linewidth = 0.2) +
  geom_sf(data = uk_map, fill = "grey85", colour = "grey50", linewidth = 0.3) +
  geom_point(data = site_coords,
             aes(x = Long, y = Lat, fill = ID),
             shape = 21, size = 5, colour = "black", stroke = 0.6) +
  scale_fill_manual(values = col_palette, breaks = legend_order, guide = "none") +
  annotation_scale(location = "bl", width_hint = 0.3,
                   bar_cols = c("black", "white"), text_cex = 0.6) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style  = north_arrow_fancy_orienteering(text_size = 7),
                         height = unit(0.8, "cm"), width = unit(0.8, "cm")) +
  coord_sf(xlim = c(-11, 3), ylim = c(49.5, 61.5), expand = FALSE) +
  labs(tag = "A") +
  common_theme +
  theme(axis.text        = element_blank(),
        axis.ticks       = element_blank(),
        axis.title       = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        plot.background  = element_rect(fill = "white", colour = NA),
        plot.margin      = margin(0, 0, 0, 0))

# ==============================================================================
# PANEL B - Pairwise FST + observed heterozygosity (Ho)
# ==============================================================================

# Weir & Cockerham FST with 10 000 bootstraps and 95 % CIs
fst_pop <- gl.fst.pop(genotypes, nboots = 10000, percent = 95, nclusters = NCORES)

write.table(fst_pop[["Fsts"]],
            file.path(BASE_DIR, "fst_pop.txt"),
            sep = "\t", row.names = TRUE, col.names = TRUE)

write.table(fst_pop[["Pvalues"]],
            file.path(BASE_DIR, "fst_pvalue_pop.txt"),
            sep = "\t", row.names = TRUE, col.names = TRUE)

# Per-population het files
het_files <- list.files(HET_DIR, pattern = "\\.het$", full.names = TRUE)

pop_data <- lapply(het_files, function(f) {
  pop  <- gsub(".*/(.*?)\\.het$", "\\1", f)
  dat  <- read.table(f, header = TRUE)
  dat$Population <- pop
  dat
}) %>%
  do.call(rbind, .) %>%
  mutate(
    O_HET      = N_SITES - O.HOM.,
    p_O_HET    = O_HET / N_SITES,      # Observed heterozygosity (Ho)
    E_HET      = N_SITES - E.HOM.,
    p_E_HET    = E_HET / N_SITES,      # Expected heterozygosity (He)
    Fis        = F,                    
    Population = factor(Population, levels = rev(legend_order)))

# Reverse order for left-to-right display matching the image (CDB -> SHE)
display_order       <- rev(legend_order)
het_palette_display <- col_palette[display_order]

pB <- ggplot(
  pop_data %>% mutate(Population = factor(Population, levels = display_order)),
  aes(x = Population, y = p_O_HET, fill = Population)) +
  geom_half_boxplot(
    side = "l", outlier.shape = NA, colour = "black", alpha = 0.8) +
  geom_jitter(
    position = position_nudge(x = 0.3), size = 3, alpha = 0.8,
    shape = 21, colour = "black", aes(fill = Population)) +
  scale_fill_manual(values = het_palette_display, guide = "none") +
  common_theme +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(size = 10, angle = 0, hjust = 0.5)) +
  labs(x = "", y = "Heterozygosity", tag = "B")

# ==============================================================================
# PANEL C - DAPC (LD1 vs LD2)
# ==============================================================================

dapc1   <- dapc(genotypes, genotypes@pop, n.pca = 6, n.da = 6)
percent <- dapc1$eig / sum(dapc1$eig) * 100

ld1_label <- paste0("LD1 (", round(percent[1], 1), " %)")
ld2_label <- paste0("LD2 (", round(percent[2], 1), " %)")

ind_coords <- as.data.frame(dapc1$ind.coord) %>%
  setNames(paste0("Axis", 1:ncol(.))) %>%
  mutate(
    Ind  = indNames(genotypes),
    Site = factor(genotypes@pop, levels = legend_order))

# Centroids
centroid <- ind_coords %>%
  group_by(Site) %>%
  summarise(Axis1.cen = mean(Axis1), Axis2.cen = mean(Axis2), .groups = "drop")

ind_coords <- left_join(ind_coords, centroid, by = "Site")

pC <- ggplot(ind_coords, aes(x = Axis1, y = Axis2)) +
  stat_ellipse(aes(fill = Site), geom = "polygon", alpha = 0.10, show.legend = FALSE) +
  geom_segment(
    aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site),
    linewidth = 0.4, show.legend = FALSE) +
  geom_point(aes(fill = Site), shape = 21, size = 3, colour = "black", alpha = 0.8) +
  scale_fill_manual(values = col_palette, name = "Location", breaks = legend_order) +
  scale_colour_manual(values = col_palette, name = "Location", breaks = legend_order) +
  labs(x = ld1_label, y = ld2_label, tag = "C") +
  common_theme +
  theme(legend.position = "none")

# ==============================================================================
# PANEL D - Isolation by Distance (IBD)
# ==============================================================================
# DistvsFst.txt is a MANUAL input (not computed here): the Fst column is taken
# from  pairwise dartR gl.fst.pop run and paired by hand with minimal
# marine distances (Dist, km) computed following Assis et al. (2013).
# Required columns: Comparison, Fst, Dist, Labels.

DistvsFst <- read_delim(DISTVSFST_FILE, delim = "\t",
                        escape_double = FALSE, trim_ws = TRUE)

# Split the Comparison column and compute linearised genetic distance
DistvsFst <- DistvsFst %>%
  mutate(Comparison = gsub("vs", "-vs-", Comparison)) %>%
  separate(Comparison, into = c("Pop1", "Pop2"), sep = "-vs-") %>%
  mutate(GD = Fst / (1 - Fst))

# Build symmetric distance matrices for the Mantel test
pop_levels <- sort(unique(c(DistvsFst$Pop1, DistvsFst$Pop2)))
n          <- length(pop_levels)
geo_mat    <- matrix(NA, nrow = n, ncol = n, dimnames = list(pop_levels, pop_levels))
gen_mat    <- matrix(NA, nrow = n, ncol = n, dimnames = list(pop_levels, pop_levels))

for (i in seq_len(nrow(DistvsFst))) {
  p1 <- DistvsFst$Pop1[i]
  p2 <- DistvsFst$Pop2[i]
  geo_mat[p1, p2] <- geo_mat[p2, p1] <- DistvsFst$Dist[i]
  gen_mat[p1, p2] <- gen_mat[p2, p1] <- DistvsFst$GD[i]
}

# Run Mantel test
geo_dist  <- as.dist(geo_mat)
gen_dist  <- as.dist(gen_mat)
mantel_r  <- mantel(geo_dist, gen_dist, method = "pearson", permutations = 999)

# Annotation text
r       <- mantel_r$statistic[1]
p_value <- mantel_r$signif[1]
annotation_text <- paste0(
  "r = ", round(r, 3), "\n",
  if (p_value < 0.05) "p < 0.05" else paste0("p = ", format.pval(p_value, digits = 3)))

# Data frame for plotting
df_dist <- data.frame(geodistance = DistvsFst$Dist,
gendistance = DistvsFst$GD,Fst = DistvsFst$Fst, Labels = DistvsFst$Labels)

pD <- ggplot(df_dist, aes(x = geodistance, y = gendistance)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              color = "grey", linetype = "dashed", alpha = 0.1) +
  geom_point(shape = 21, alpha = 0.8, stroke = 1, size = 3,
             colour = "black", fill = "grey20") +
  xlab("Geographic distance (km)") +
  ylab("Genetic distance (Fst / 1 - Fst)") +
  annotate("text",
           x = max(df_dist$geodistance), y = min(df_dist$gendistance),
           label = annotation_text, hjust = 1, vjust = 0, size = 4) +
  labs(tag = "D") +
  common_theme +
  theme(legend.position = "none")

# ==============================================================================
# SHARED LEGEND
# ==============================================================================

legend_df <- tibble(x = 1:8, y = 1:8, Location = factor(legend_order, levels = legend_order))

legend_plot <- ggplot(legend_df, aes(x, y, fill = Location)) +
  geom_point(shape = 21, size = 4, colour = "black") +
  scale_fill_manual(values = col_palette, breaks = legend_order,
                    labels = legend_labels,
                    name   = "Sampling\nlocation") +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme_void() +
  theme(
    legend.title    = element_text(size = 14, face = "bold"),
    legend.text     = element_text(size = 12),
    legend.key.size = unit(0.7, "cm"))

shared_legend <- cowplot::get_legend(legend_plot)

# ==============================================================================
# COMPOSE FINAL FIGURE
# ==============================================================================

right_column <- pB / pC / pD + plot_layout(heights = c(1, 1.3, 1))

composite <- (pA | right_column) + plot_layout(widths = c(1.1, 1))

final_fig <- cowplot::plot_grid(composite, shared_legend,
rel_widths = c(1, 0.08), nrow = 1)

# -- Save ----------------------------------------------------------------------
ggsave(file.path(FIG_DIR, "Figure_3.png"), plot = final_fig, width = 18, height = 12, dpi = 600, bg = "white")

ggsave(file.path(FIG_DIR, "Figure_3.pdf"), plot = final_fig, width = 18, height = 12, bg = "white")

cat("Figure saved to", FIG_DIR, "\n")

EOF
