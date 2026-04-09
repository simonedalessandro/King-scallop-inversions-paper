#!/bin/sh
# ==============================================================================
# Script:       12a.het_miss.sh
# Description:  Compute per-individual heterozygosity and missing data rates
#               from the full ANGSD-called SNP dataset (pre-pruning) using PLINK.
#               Output .het and .imiss files are used for individual-level QC
#               to flag potential cross-contamination and low-quality samples.
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE job configuration
# ------------------------------------------------------------------------------
#$ -N het_miss
#$ -cwd
#$ -l h_rt=02:00:00
#$ -l h_rss=8G
#$ -pe sharedmem 1
#$ -o o_files
#$ -e e_files

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load plink/1.90b7.2

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

INPUT=/simone/pmax2023/out/10.angsd/Pmax_168_merged_chr
OUTPUT=/simone/pmax2023/out/10.angsd/Pmax_168_merged_chr

# ------------------------------------------------------------------------------
# Compute heterozygosity and missingness
# ------------------------------------------------------------------------------

plink \
    --bfile  "${INPUT}" \
    --allow-extra-chr \
    --missing \
    --het \
    --out    "${OUTPUT}"

echo "[$(date)] Done. Output: ${OUTPUT}.het / ${OUTPUT}.imiss"

# ------------------------------------------------------------------------------
# Next step:
#   Download .het and .imiss to local machine and run 12b.het_miss.R for QC plots.
# ------------------------------------------------------------------------------
