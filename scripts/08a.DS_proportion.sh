#!/bin/sh
# ==============================================================================
# Script:       08a.DS_proportions.sh
# Description:  Compute per-sample proportions for downsampling high-coverage BAM files to uniform 6x using samtools view -s.
#
# Input:        Per-sample *_cov.gz files from 07.bamcov.sh
# Output:       proportion_to_downsample_10x.txt in file_lists/
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE job configuration
# ------------------------------------------------------------------------------
#$ -N DS_prop
#$ -cwd
#$ -l h_rt=02:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -o o_files
#$ -e e_files

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load R/4.5

# ------------------------------------------------------------------------------
# Run R code
# ------------------------------------------------------------------------------

Rscript - <<'EOF'

library(data.table)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)

options(scipen = 999)

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

data_path   <- "/simone/pmax2023/out/07coverage"
output_file <- "/simone/pmax2023/out/file_lists/proportion_to_downsample_10x.txt"

# ------------------------------------------------------------------------------
# Load samtools coverage output
# ------------------------------------------------------------------------------

files <- dir(data_path, pattern = "*cov.gz")

cov <- tibble(filename = files) %>%
  mutate(file_contents = map(filename, ~ fread(file.path(data_path, .)))) %>%
  unnest(cols = c(file_contents))

# ------------------------------------------------------------------------------
# Filter to P. maximus chromosomes (LR736838.1 - LR736856.1)
# ------------------------------------------------------------------------------

pmax_scaffolds <- c(
  paste0("LR73683", 8:9, ".1"),
  paste0("LR73684", 0:9, ".1"),
  paste0("LR73685", 0:6, ".1"))

cov <- cov %>%
  filter(`#rname` %in% pmax_scaffolds)

# ------------------------------------------------------------------------------
# Clean sample names
# ------------------------------------------------------------------------------

cov <- cov %>%
  mutate(filename = gsub("_mapped_sorted_RG_rmdup_cov.gz", "", filename))

# ------------------------------------------------------------------------------
# Diagnostic plots: depth and breadth distributions
# ------------------------------------------------------------------------------

hist(cov$meandepth, main = "Mean depth per scaffold", xlab = "Mean depth")
hist(cov$coverage,  main = "Breadth of coverage per scaffold", xlab = "Coverage (%)")

# ------------------------------------------------------------------------------
# Mean depth per sample across autosomes
# ------------------------------------------------------------------------------

depth_summary <- cov %>%
  group_by(filename) %>%
  summarise(
    mean_depth = mean(meandepth),
    min_depth  = min(meandepth),
    max_depth  = max(meandepth)
  )

print(depth_summary)

# ------------------------------------------------------------------------------
# Compute downsampling proportion
# Target depth: 6x
# Proportion = 6 / observed_mean_depth
# samtools view -s uses integer part as seed and decimal as fraction,
# so 0.XX is replaced with 1.XX to set seed = 1
# ------------------------------------------------------------------------------

subsample <- depth_summary %>%
  mutate(
    prop_sub = 6 / mean_depth,
    prop_sub = round(prop_sub, 2),
    prop_sub = gsub("0\\.", "1.", as.character(prop_sub))
  )

# ------------------------------------------------------------------------------
# Write output: two columns — sample_name  proportion
# ------------------------------------------------------------------------------

write.table(
  subsample[c("filename", "prop_sub")],
  file      = output_file,
  quote     = FALSE,
  col.names = FALSE,
  row.names = FALSE
)

cat("Written:", output_file, "\n")
cat("Samples to downsample:", nrow(subsample), "\n")

EOF
