#!/bin/sh
# ==============================================================================
# Script:       22.het.sh
# Description:  Compute per-individual observed heterozygosity for each putative
#               inversion region using VCFtools --het.
#               Output .het files are used in 24.MDS_pca_het_LD_plot.sh to
#               compare heterozygosity across chromosomal arrangement clusters.
#
# Input:        Per-region VCF files (<REGION>.vcf) from 21a.local_pca.sh
#               putative_inversions_list.txt — one region name per line
#
# Output:       Per-region <REGION>.het files in OUT_DIR
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE array job configuration
# ------------------------------------------------------------------------------
#$ -N het
#$ -cwd
#$ -l h_rt=6:00:00
#$ -l h_rss=12G
#$ -t 1-14
#$ -pe sharedmem 4
#$ -o o_files
#$ -e e_files

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load vcftools/0.1.16

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

INPUT_DIR=/simone/pmax2023/out/20.lostruct/inversions
OUT_DIR=/simone/pmax2023/out/21.MDS_pca_het_LD/het
REGION_LIST=/simone/pmax2023/out/file_lists/putative_inversions_list.txt

mkdir -p "${OUT_DIR}"

# ------------------------------------------------------------------------------
# Parse region name for this array task
# ------------------------------------------------------------------------------

REGION=$(sed -n "${SGE_TASK_ID}p" "${REGION_LIST}")
PREFIX=${OUT_DIR}/${REGION}

echo "=================================================="
echo "Task ID:    ${SGE_TASK_ID}"
echo "Region:     ${REGION}"
echo "Input VCF:  ${INPUT_DIR}/${REGION}.vcf"
echo "Start time: $(date)"
echo "=================================================="

# ------------------------------------------------------------------------------
# Compute per-individual heterozygosity
# ------------------------------------------------------------------------------

vcftools --vcf "${INPUT_DIR}/${REGION}.vcf" --het --out "${PREFIX}"

echo "=================================================="
echo "Finished region: ${REGION}"
echo "End time:        $(date)"
echo "=================================================="
