#!/bin/sh
# ==============================================================================
# Script:       03.bwa-mem.sh
# Description:  Align trimmed paired-end reads to the P. maximus reference
# genome using BWA-MEM. Filters unmapped reads with SAMtools.
#
# Input:        Trimmed FASTQ pairs from 02cutadapt (output of 01.TrimGalore.sh)
# Output:       Per-sample mapped BAM files
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE array job configuration
# ------------------------------------------------------------------------------
#$ -N bwa-mem
#$ -cwd
#$ -l h_rt=12:00:00
#$ -l h_rss=8G
#$ -t 1-168
#$ -pe sharedmem 6
#$ -o o_files
#$ -e e_files

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load bwa/0.7.17
module load samtools/1.9    

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

TARGET_DIR=/simone/pmax2023/out/02cutadapt
OUTPUT_DIR=/simone/pmax2023/out/03mapped_bams
SAMPLE_SHEET="/simone/pmax2023/out/file_lists/01.bwa-mem.sample.list.txt"

REFERENCE="/simone/reference/Pmax/GCA_902652985.1_xPecMax1.1_genomic.fna"

# ------------------------------------------------------------------------------
# Parse sample information for this array task
# ------------------------------------------------------------------------------

base=$(sed -n "${SGE_TASK_ID}p" "$SAMPLE_SHEET" | awk '{print $1}')
r1=$(sed -n "${SGE_TASK_ID}p"   "$SAMPLE_SHEET" | awk '{print $2}') 
r2=$(sed -n "${SGE_TASK_ID}p"   "$SAMPLE_SHEET" | awk '{print $3}')  

echo "=================================================="
echo "Task ID:    ${SGE_TASK_ID}"
echo "Sample:     ${base}"
echo "R1:         ${r1}"
echo "R2:         ${r2}"
echo "Host:       ${HOSTNAME}"
echo "Start time: $(date)"
echo "=================================================="

# ------------------------------------------------------------------------------
# Align reads to reference genome
# ------------------------------------------------------------------------------
# Key parameters:
#   -bF 4       Output BAM (-b), excluding unmapped reads (SAMtools flag 4)

bwa mem \
    -t 6 \
    "${REFERENCE}" \
    "${TARGET_DIR}/${r1}" \
    "${TARGET_DIR}/${r2}" \
    | samtools view -bF 4 - > "${OUTPUT_DIR}/${base}_mapped.bam"

echo "=================================================="
echo "Finished sample: ${base}"
echo "End time:        $(date)"
echo "=================================================="
