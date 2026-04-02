#!/bin/sh
# ==============================================================================
# Script:       06.MarkDup.sh
# Description:  Mark and remove PCR duplicates using Picard MarkDuplicates.
#               Extracts per-sample duplication statistics to a TSV file.
#
# Input:        Read-group-tagged BAM files from 05.ReadGroups.sh
# Output:       Deduplicated BAM files, Picard metrics, and per-sample stats TSV
# ==============================================================================
 
# ------------------------------------------------------------------------------
# SGE array job configuration
# ------------------------------------------------------------------------------
#$ -N MarkDups
#$ -cwd
#$ -l h_rt=4:00:00
#$ -l h_rss=8G
#$ -pe sharedmem 4
#$ -t 1-168
#$ -o o_files
#$ -e e_files
 
# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------
 
. /etc/profile.d/modules.sh
 
module load java
module load samtools/1.20
 
# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------
 
INPUT_DIR=/simone/pmax2023/out/05group_bams
OUTPUT_DIR=/simone/pmax2023/out/06rmdup_bams
METRICS_DIR=/simone/pmax2023/out/06rmdup_bams/metrics
STATS_DIR=/simone/pmax2023/out/06rmdup_bams/stats
 
# Job-specific scratch directory to avoid collisions between array tasks
SCRATCH=/scratch/sdaless/picard_tmp/${JOB_ID}_${SGE_TASK_ID}
 
PICARD=/software/picard/picard.jar
 
mkdir -p "${OUTPUT_DIR}" "${METRICS_DIR}" "${STATS_DIR}" "${SCRATCH}"
 
# ------------------------------------------------------------------------------
# Parse sample for this array task
# ------------------------------------------------------------------------------
 
BAM_LIST=$(ls -1 "${INPUT_DIR}"/*_RG.bam)
THIS_BAM=$(echo "${BAM_LIST}" | sed -n "${SGE_TASK_ID}p")
 
if [ -z "${THIS_BAM}" ]; then
    echo "ERROR: No BAM found for task ${SGE_TASK_ID}" >&2
    exit 1
fi
 
BASE=$(basename "${THIS_BAM}" .bam)
SAMPLE=$(basename "${THIS_BAM}" | cut -d '_' -f 1)
OUTPUT_BAM="${OUTPUT_DIR}/${BASE}_rmdup.bam"
METRICS="${METRICS_DIR}/${BASE}_rmdup.metrics"
 
echo "=================================================="
echo "Task ID:    ${SGE_TASK_ID}"
echo "Sample:     ${SAMPLE}"
echo "Input:      ${THIS_BAM}"
echo "Output:     ${OUTPUT_BAM}"
echo "Host:       ${HOSTNAME}"
echo "Start time: $(date)"
echo "=================================================="
 
# ------------------------------------------------------------------------------
# Mark and remove duplicates
# ------------------------------------------------------------------------------
 
java -Xmx14g -jar "${PICARD}" MarkDuplicates \
    I="${THIS_BAM}" \
    O="${OUTPUT_BAM}" \
    METRICS_FILE="${METRICS}" \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=true \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
    MAX_RECORDS_IN_RAM=500000 \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR="${SCRATCH}"
 
if [ $? -ne 0 ]; then
    echo "ERROR: MarkDuplicates failed for ${THIS_BAM}" >&2
    rm -rf "${SCRATCH}"
    exit 1
fi
 
# Clean up job-specific temp directory
rm -rf "${SCRATCH}"
 
# ------------------------------------------------------------------------------
# Extract duplication statistics from Picard metrics file
# ------------------------------------------------------------------------------
 
# Skip comment lines and blank lines; second data line contains library metrics
METRICS_LINE=$(grep -v "^#" "${METRICS}" | grep -v "^$" | awk 'NR==2')
 
unpaired_examined=$(echo "${METRICS_LINE}" | awk '{print $2}')
pairs_examined=$(echo    "${METRICS_LINE}" | awk '{print $3}')
unpaired_dups=$(echo     "${METRICS_LINE}" | awk '{print $6}')
pair_dups=$(echo         "${METRICS_LINE}" | awk '{print $7}')
optical_dups=$(echo      "${METRICS_LINE}" | awk '{print $8}')
pct_dup=$(echo           "${METRICS_LINE}" | awk '{printf "%.4f", $9}')
est_lib_size=$(echo      "${METRICS_LINE}" | awk '{print $10}')
 
# Count mapped non-duplicate reads in the output BAM (flags 2052: unmapped + duplicate)
reads_after_dedup=$(samtools view -c -F 2052 "${OUTPUT_BAM}")
 
# Derive total reads entering dedup from Picard counts
reads_before_dedup=$(awk -v pe="${pairs_examined}" -v up="${unpaired_examined}" \
    'BEGIN{print (pe*2) + up}')
 
# Write per-sample stats TSV
# Columns: sample | reads_before | reads_after | pct_dup | unpaired_dups | pair_dups | optical_dups | est_lib_size
echo -e "${SAMPLE}\t${reads_before_dedup}\t${reads_after_dedup}\t${pct_dup}\t${unpaired_dups}\t${pair_dups}\t${optical_dups}\t${est_lib_size}" \
    > "${STATS_DIR}/${SAMPLE}_dedup.stats.tsv"
 
echo "=================================================="
echo "Finished:            ${SAMPLE}"
echo "Reads before dedup:  ${reads_before_dedup}"
echo "Reads after dedup:   ${reads_after_dedup}"
echo "Duplication rate:    ${pct_dup}"
echo "Optical duplicates:  ${optical_dups}"
echo "Est. library size:   ${est_lib_size}"
echo "End time:            $(date)"
echo "=================================================="
