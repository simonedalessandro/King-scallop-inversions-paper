#!/bin/sh
# ==============================================================================
# Script:       05.ReadGroups.sh
# Description:  Add read group tags to sorted BAM files using Picard AddOrReplaceReadGroups. #
# Input:        Coordinate-sorted BAM files from 04.SortSam.sh
# Output:       BAM files with read group tags in OUTPUT_DIR
==============================================================================

# ------------------------------------------------------------------------------
# SGE array job configuration
# ------------------------------------------------------------------------------
#$ -N addRG
#$ -cwd
#$ -l h_rt=08:00:00
#$ -l h_rss=4G
#$ -pe sharedmem 2
#$ -t 1-168
#$ -o o_files
#$ -e e_files

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load java

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

INPUT_DIR=/simone/pmax2023/out/04sorted_bams
OUTPUT_DIR=/simone/pmax2023/out/05group_bams

# Scratch space for Picard temporary files
SCRATCH=/scratch/sdaless

# ------------------------------------------------------------------------------
# Parse sample for this array task
# ------------------------------------------------------------------------------

BAM=$(ls -1 "${INPUT_DIR}"/*mapped_sorted.bam)
THIS_BAM=$(echo "${BAM}" | sed -n "${SGE_TASK_ID}p")
BASE=$(basename "${THIS_BAM}" .bam)

# Extract sample name from filename (used for all read group fields)
RGSM=$(echo "${BASE}" | cut -f 1 -d '_')
RGID="${RGSM}"
RGPU="${RGSM}"

echo "=================================================="
echo "Task ID:    ${SGE_TASK_ID}"
echo "Input BAM:  ${THIS_BAM}"
echo "Sample:     ${RGSM}"
echo "Host:       ${HOSTNAME}"
echo "Start time: $(date)"
echo "=================================================="

# ------------------------------------------------------------------------------
# Add read groups
# ------------------------------------------------------------------------------

java -Xmx4g -jar /software/picard/picard.jar AddOrReplaceReadGroups \
    I="${THIS_BAM}" \
    O="${OUTPUT_DIR}/${BASE}_RG.bam" \
    RGID="${RGID}" \
    RGPL=illumina \
    RGLB=lib1 \
    RGPU="${RGPU}" \
    RGSM="${RGSM}" \
    VALIDATION_STRINGENCY=SILENT \
    SORT_ORDER=coordinate \
    TMP_DIR="${SCRATCH}"

echo "=================================================="
echo "Finished: ${BASE}"
echo "End time: $(date)"
echo "=================================================="
