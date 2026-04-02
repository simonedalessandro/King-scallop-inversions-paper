#!/bin/sh
# ==============================================================================
# Script:       07.bamcov.sh
# Description:  Estimate depth and breadth of coverage per sample after duplicate removal using SAMtools depth and coverage.
#
# Input: Deduplicated BAM files from 06.MarkDup.sh
# Output: Per-sample depth.gz, coverage.gz, and summary TSV in OUTPUT_DIR
#
# Note: After all tasks complete, merge per-sample summary TSVs into a
#             single file using: awk 'FNR>1 || NR==1' OUTPUT_DIR/*_summary.txt
#             > all_samples_coverage_summary.txt
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE array job configuration
# ------------------------------------------------------------------------------
#$ -N bamcov
#$ -cwd
#$ -l h_rt=8:00:00
#$ -l h_rss=10G
#$ -pe sharedmem 4
#$ -t 1-168
#$ -o o_files
#$ -e e_files

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load samtools/1.10

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

TARGET_DIR=/simone/pmax2023/out/06rmdup_bams
OUTPUT_DIR=/simone/pmax2023/out/07coverage

mkdir -p "${OUTPUT_DIR}"

# ------------------------------------------------------------------------------
# Parse sample for this array task
# ------------------------------------------------------------------------------

BAM=$(ls -1 "${TARGET_DIR}"/*.bam)
THIS_BAM=$(echo "${BAM}" | sed -n "${SGE_TASK_ID}p")
BASE=$(basename "${THIS_BAM}" .bam | awk -F'_mapped' '{print $1}')

echo "=================================================="
echo "Task ID:    ${SGE_TASK_ID}"
echo "Sample:     ${BASE}"
echo "Input BAM:  ${THIS_BAM}"
echo "Host:       ${HOSTNAME}"
echo "Start time: $(date)"
echo "=================================================="

# ------------------------------------------------------------------------------
# Compute depth and coverage
# ------------------------------------------------------------------------------
# Key parameters:
#   samtools depth -a   Include all positions (even zero-coverage sites) for breadth calc
#   samtools coverage   Per-chromosome summary including breadth and mean depth

# Save full per-position depth and per-chromosome coverage (compressed)
samtools depth "${THIS_BAM}"    | gzip > "${OUTPUT_DIR}/${BASE}_depth.gz"
samtools coverage "${THIS_BAM}" | gzip > "${OUTPUT_DIR}/${BASE}_cov.gz"

# Breadth of coverage: proportion of reference positions covered by >= 1 read
coverage=$(samtools depth -a "${THIS_BAM}" | \
    awk '{if ($3 > 0) covered++} END {print covered/NR * 100}')

# Mean depth: average sequencing depth across all reference positions
depth=$(samtools depth "${THIS_BAM}" | \
    awk '{sum += $3} END {if (NR > 0) print sum/NR; else print 0}')

# Write per-sample summary TSV
echo -e "ID\tCoverage_Breadth\tDepth" > "${OUTPUT_DIR}/${BASE}_summary.txt"
echo -e "${BASE}\t${coverage}%\t${depth}" >> "${OUTPUT_DIR}/${BASE}_summary.txt"

echo "=================================================="
echo "Finished:          ${BASE}"
echo "Coverage breadth:  ${coverage}%"
echo "Mean depth:        ${depth}x"
echo "End time:          $(date)"
echo "=================================================="
