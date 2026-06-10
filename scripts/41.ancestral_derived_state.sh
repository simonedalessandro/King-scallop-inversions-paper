#!/bin/sh
# ==============================================================================
# Script:       41_ancestral_derived_state.sh
# Description:  Assign ancestral / derived state to inversion genotypes and
#               visualise the results. Two inversion sets are handled:
#                 - high_confidence
#                 - candidate
#
#               Figures:
#                 1. anc_het_der_boxplot_inversions(.png/.pdf)
#                      Per-individual % ancestral / het / derived SNPs by karyotype.
#                 2. anc_het_der_boxplot_inversions_candidate(.png/.pdf)
#                      Same boxplot for the candidate set.
#                 3. inversion_genotype_matrix(.png/.pdf)
#                      Individual x SNP genotype tile matrix, high_confidence.
#                 4. inversion_karyotype_frequency(.png/.pdf)
#                      Per-population karyotype frequency, high_confidence.
#
# Input:        Per-LGC .traw + <LGC>_IDs_loc_clusters.txt under ${base_path}/<LGC>/
#               .traw files are produced by 40_polarise.sh.
#               traw coding (REF = ancestral): 2 = A/A, 1 = A/D, 0 = D/D.
#
# Output:       Figures written under ${base_path}.
# ==============================================================================

. /etc/profile.d/modules.sh

module load R/4.5

Rscript - <<'EOF'

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(gghalves)
})

# ==============================================================================
# Settings
# ==============================================================================

# Polarised, high-confidence inversions
high_confidence <- c("LGC02", "LGC13.01", "LGC15", "LGC17.01",
                     "LGC18.01", "LGC18.03", "LGC19.01")

# Candidate inversions that did not meet criteria
candidate <- c("LGC08.01", "LGC09", "LGC11", "LGC13.02", "LGC17.02", "LGC18.04")

base_path <- "local/path"

loc_order        <- c("CDB", "TAR", "CLY", "FRA", "MOR", "MNC", "ORK", "SHE")
karyotype_levels <- c("Ancestral homokaryotype",
                      "Heterokaryotype",
                      "Derived homokaryotype")

# Karyotype cluster: 1 = ancestral homokaryotype, 2 = heterokaryotype, 3 = derived homokaryotype
col_palette  <- c("1" = "#00798C", "2" = "#D1495B", "3" = "#EDAE49")
geno_palette <- c("Ancestral homokaryotype" = "#00798C",
                  "Heterokaryotype"         = "#D1495B",
                  "Derived homokaryotype"   = "#EDAE49")

shared_theme <- theme_classic(base_size = 11) +
  theme(
    strip.background  = element_rect(fill = "grey85", color = "grey50", linewidth = 0.4),
    strip.text        = element_text(size = 10, face = "bold", color = "grey20",
                                     margin = margin(3, 3, 3, 3)),
    axis.title.y      = element_text(size = 12, margin = margin(r = 6)),
    axis.title.x      = element_text(size = 12, margin = margin(t = 6)),
    axis.text         = element_text(size = 10, color = "grey20"),
    axis.line         = element_line(color = "grey40", linewidth = 0.3),
    axis.ticks        = element_line(color = "grey40", linewidth = 0.3),
    panel.spacing     = unit(0.4, "cm"),
    panel.border      = element_rect(color = "grey70", fill = NA, linewidth = 0.3),
    panel.background  = element_rect(fill = "white"),
    legend.title      = element_text(size = 11, face = "bold"),
    legend.text       = element_text(size = 10),
    legend.key.size   = unit(0.45, "cm"),
    legend.position   = "right",
    plot.margin       = unit(c(0.4, 0.3, 0.4, 0.3), "cm"))

# ==============================================================================
# Helpers
# ==============================================================================
# load_lgc(): read one LGC's .traw + metadata, return long format
#   (SNP, name, value + metadata + dataset). PLINK2 prefixes sample columns
#   with the FID ("0_"), which is stripped to recover the sample IDs.
load_lgc <- function(dataset_name) {
  traw_file <- file.path(base_path, dataset_name, paste0(dataset_name, ".traw"))
  meta_file <- file.path(base_path, dataset_name,
                         paste0(dataset_name, "_IDs_loc_clusters.txt"))

  meta <- fread(meta_file) %>%
    mutate(Cluster = factor(Cluster, levels = c("1", "2", "3")))

  fread(traw_file) %>%
    select(-any_of(c("CHR", "(C)M", "POS", "COUNTED", "ALT"))) %>%
    pivot_longer(cols = -SNP, names_to = "name", values_to = "value") %>%
    mutate(name = sub("^0_", "", name)) %>%
    left_join(meta, by = c("name" = "Sample")) %>%
    mutate(dataset = dataset_name)
}

# Per-individual genotype counts (one row per individual per LGC)
summarise_indiv <- function(long) {
  long %>%
    filter(!is.na(value)) %>%
    group_by(dataset, name, Location, Cluster) %>%
    summarise(
      total_genotype  = n(),
      total_ancestral = sum(value == 2),
      total_het       = sum(value == 1),
      total_derived   = sum(value == 0),
      .groups = "drop"
    )
}

# One boxplot column (half-point + half-boxplot)
create_plot <- function(data, y_var, y_lab) {
  ggplot(data, aes(x = Cluster, y = .data[[y_var]], fill = Cluster)) +
    geom_half_point(
      side = "l", shape = 21,
      alpha = 0.55, stroke = 0.15, size = 1.8
    ) +
    geom_half_boxplot(
      side = "r",
      outlier.color = NA,
      width = 0.35,
      linewidth = 0.35,
      color = "grey20",
      alpha = 0.85,
      notch = FALSE
    ) +
    scale_fill_manual(values = col_palette, name = "Inversion \nkaryotype") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
    labs(y = y_lab, x = "Inversion karyotype") +
    facet_wrap(~dataset, ncol = 1, strip.position = "top") +
    shared_theme
}

# Full 3-column boxplot figure from pre-loaded long data
make_boxplot_figure <- function(long) {
  prop_data <- summarise_indiv(long) %>%
    mutate(
      prop_ancestral = total_ancestral / total_genotype * 100,
      prop_het       = total_het       / total_genotype * 100,
      prop_derived   = total_derived   / total_genotype * 100
    )

  p1 <- create_plot(prop_data, "prop_ancestral", "% of ancestral SNPs")
  p2 <- create_plot(prop_data, "prop_het",       "% of heterozygous SNPs")
  p3 <- create_plot(prop_data, "prop_derived",   "% of derived SNPs")

  p1 + p2 + p3 +
    plot_layout(guides = "collect", widths = c(1, 1, 1, 1.1), nrow = 1) &
    theme(legend.position = "right",
          axis.title.x = element_blank())
}

# ==============================================================================
# Load both inversion sets
# ==============================================================================

hc_long <- bind_rows(lapply(high_confidence, load_lgc)) %>%
  mutate(dataset = factor(dataset, levels = high_confidence))

cand_long <- bind_rows(lapply(candidate, load_lgc)) %>%
  mutate(dataset = factor(dataset, levels = candidate))

# ==============================================================================
# FIGURE S23: anc/het/derived boxplot — high-confidence inversions
# ==============================================================================

fig_hc <- make_boxplot_figure(hc_long)

ggsave(file.path(base_path, "anc_het_der_boxplot_inversions.png"),
       plot = fig_hc, width = 16, height = 16, dpi = 300, bg = "white")
ggsave(file.path(base_path, "anc_het_der_boxplot_inversions.pdf"),
       plot = fig_hc, width = 16, height = 16, dpi = 300, bg = "white")

cat("Saved: anc_het_der_boxplot_inversions.png/.pdf\n")

# ==============================================================================
# FIGURE S24: anc/het/derived boxplot — candidate inversions
# ==============================================================================

fig_cand <- make_boxplot_figure(cand_long)

ggsave(file.path(base_path, "anc_het_der_boxplot_inversions_candidate.png"),
       plot = fig_cand, width = 16, height = 14, dpi = 300, bg = "white")
ggsave(file.path(base_path, "anc_het_der_boxplot_inversions_candidate.pdf"),
       plot = fig_cand, width = 16, height = 14, dpi = 300, bg = "white")

cat("Saved: anc_het_der_boxplot_inversions_candidate.png/.pdf\n")

# ==============================================================================
# FIGURE S25: karyotype frequency per population — high-confidence
# ==============================================================================

hc_indiv <- summarise_indiv(hc_long)

freq_data <- hc_indiv %>%
  mutate(
    Location = recode(Location, "MFR" = "MOR"),
    Location = factor(Location, levels = loc_order),
    Genotype = case_when(
      Cluster == "1" ~ "Ancestral homokaryotype",
      Cluster == "2" ~ "Heterokaryotype",
      Cluster == "3" ~ "Derived homokaryotype"),
    Genotype = factor(Genotype, levels = karyotype_levels)
  ) %>%
  group_by(dataset, Location, Genotype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(dataset, Location) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

p_bar_map <- ggplot(freq_data,
                    aes(x = Location, y = freq, fill = Genotype)) +
  geom_col(position = "stack", color = "grey30", linewidth = 0.2, width = 0.7) +
  scale_fill_manual(values = geno_palette,
                    name   = "Inversion\nkaryotype\n",
                    labels = c("Ancestral homokaryotype" = "Ancestral\nhomokaryotype",
                               "Heterokaryotype"          = "Heterokaryotype",
                               "Derived homokaryotype"    = "Derived\nhomokaryotype")) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.03))) +
  facet_wrap(~dataset, ncol = 1, strip.position = "right") +
  labs(x = "Population", y = "Inversion karyotype frequency") +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey85", color = "grey50", linewidth = 0.4),
    strip.text       = element_text(face = "bold", size = 10, color = "grey20"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.title       = element_text(face = "bold", size = 16),
    panel.spacing    = unit(0.3, "cm"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 16))

ggsave(file.path(base_path, "inversion_karyotype_frequency.png"),
       plot = p_bar_map, width = 12, height = 14, dpi = 300, bg = "white")
ggsave(file.path(base_path, "inversion_karyotype_frequency.pdf"),
       plot = p_bar_map, width = 12, height = 14, dpi = 300, bg = "white")

cat("Saved: inversion_karyotype_frequency.png/.pdf\n")
cat("\nDone.\n")

EOF
echo "Wrote $(wc -l < /mnt/user-data/outputs/41_ancestral_derived_state.sh) lines"
