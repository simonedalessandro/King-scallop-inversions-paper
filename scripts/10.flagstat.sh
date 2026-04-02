#!/bin/sh
# ==============================================================================
# Script:       10.flagstat.sh
# Description:  Generate alignment statistics for final BAM files using SAMtools flagstat
#
# Input:        Deduplicated 6x BAM files from 06.MarkDup.sh / 08.Downsample.sh
# Output:       Per-sample flagstat text files in OUTPUT_DIR
#==============================================================================

# ------------------------------------------------------------------------------
# SGE array job configuration
# ------------------------------------------------------------------------------
#$ -N flagstat
#$ -cwd
#$ -l h_rt=2:00:00
#$ -l h_rss=4G
#$ -pe sharedmem 4
#$ -t 1-168
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

TARGET_DIR=/simone/pmax2023/out/06rmdup_bams/06x
OUTPUT_DIR=/simone/pmax2023/out/06rmdup_bams/06x/flagstat

mkdir -p "${OUTPUT_DIR}"

# ------------------------------------------------------------------------------
# Parse sample for this array task
# ------------------------------------------------------------------------------

BAM=$(ls -1 "${TARGET_DIR}"/*.bam)
THIS_BAM=$(echo "${BAM}" | sed -n "${SGE_TASK_ID}p")
BASE=$(basename "${THIS_BAM}" .bam)

echo "=================================================="
echo "Task ID:    ${SGE_TASK_ID}"
echo "Input BAM:  ${THIS_BAM}"
echo "Host:       ${HOSTNAME}"
echo "Start time: $(date)"
echo "=================================================="

# ------------------------------------------------------------------------------
# Run flagstat
# ------------------------------------------------------------------------------

samtools flagstat "${THIS_BAM}" > "${OUTPUT_DIR}/${BASE}_flagstat.txt"

echo "=================================================="
echo "Finished:   ${BASE}"
echo "End time:   $(date)"
echo "=================================================="
