#!/bin/sh
# ==============================================================================
# Script:       08.DS.bam.sh
# Description:  Downsample high-coverage (10x) BAM files to uniform 6x coverage
# using SAMtools. The downsampling proportion per sample is precomputed and provided via a # coverage proportion file.
#
# Input:        Deduplicated 10x BAM files and per-sample proportion file
# Output:       Downsampled 6x BAM files
#==============================================================================

# ------------------------------------------------------------------------------
# SGE array job configuration
# ------------------------------------------------------------------------------
#$ -N DS_bam_10x
#$ -cwd
#$ -l h_rt=8:00:00
#$ -l h_vmem=6G
#$ -t 1-36
#$ -o o_files
#$ -e e_files

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load samtools/1.9

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

INPUT_DIR="/simone/pmax2023/out/06.rmdup_bams/10x"
OUTPUT_DIR="/simone/pmax2023/out/06.rmdup_bams/06x"

# Tab-separated file: sample_name  proportion_to_downsample
# Proportion computed as: target_depth / observed_depth (e.g. 6 / 10.4 = 0.577)
COVERAGE="/simone/pmax2023/out/file_lists/proportion_to_downsample_10x.txt"

mkdir -p "${OUTPUT_DIR}"

# ------------------------------------------------------------------------------
# Parse sample for this array task
# ------------------------------------------------------------------------------

BAM=$(ls -1 "${INPUT_DIR}"/*rmdup.bam)
THIS_BAM=$(echo "${BAM}" | sed -n "${SGE_TASK_ID}p")
PROP=$(sed -n "${SGE_TASK_ID}p" "${COVERAGE}" | awk '{print $2}')

echo "=================================================="
echo "Task ID:     ${SGE_TASK_ID}"
echo "Input BAM:   ${THIS_BAM}"
echo "Proportion:  ${PROP}"
echo "Host:        ${HOSTNAME}"
echo "Start time:  $(date)"
echo "=================================================="

# ------------------------------------------------------------------------------
# Downsample BAM
# ------------------------------------------------------------------------------

samtools view \
    -s "${PROP}" \
    -b "${THIS_BAM}" \
    > "${OUTPUT_DIR}/$(basename "${THIS_BAM%.bam}")_downsample_6x.bam"

echo "=================================================="
echo "Finished:   $(basename "${THIS_BAM}")"
echo "End time:   $(date)"
echo "=================================================="
