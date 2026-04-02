#!/bin/sh
# ==============================================================================
# Script:       01.TrimGalore.sh
# Quality and adapter trimming of paired-end raw reads using
#               TrimGalore (wrapper around Cutadapt + FastQC). Runs as an SGE
#               array job, one task per sample (168 samples total)
#
# Input:        Sample sheet with columns: [base] [R1_filename] [R2_filename]
# Output:       Trimmed FASTQ files and FastQC reports in OUTPUT_DIR
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE array job configuration
# ------------------------------------------------------------------------------
#$ -N TrimGalore          
#$ -cwd                  
#$ -l h_rt=16:00:00       
#$ -l h_rss=8G           
#$ -e e_files          
#$ -o o_files           
#$ -t 1-168              
#$ -pe sharedmem 4        

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

# Initialise the module system
. /etc/profile.d/modules.sh

module load cutadapt/4.6        
module load TrimGalore/0.6.6    
module load FastQC/0.11.9      
module load pigz/2.3.3         

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

# Directory containing raw paired-end FASTQ files
TARGET_DIR=/simone/pmax2023/raw

# Sample sheet: tab-separated, three columns — sample_name  R1_file  R2_file
SAMPLE_SHEET="/simone/pmax2023/out/file_lists/00.trimgalore.sample.list.txt"

# Output directory for trimmed reads and FastQC reports
OUTPUT_DIR=/simone/pmax2023/out/02cutadapt

# ------------------------------------------------------------------------------
# Parse sample information for this array task
# ------------------------------------------------------------------------------

# Each SGE task ID corresponds to one line in the sample sheet
base=$(sed -n "${SGE_TASK_ID}p" "$SAMPLE_SHEET" | awk '{print $1}')  # Sample name / base ID
r1=$(sed -n "${SGE_TASK_ID}p"   "$SAMPLE_SHEET" | awk '{print $2}')  # R1 filename
r2=$(sed -n "${SGE_TASK_ID}p"   "$SAMPLE_SHEET" | awk '{print $3}')  # R2 filename

echo "=================================================="
echo "Task ID:        ${SGE_TASK_ID}"
echo "Sample:         ${base}"
echo "R1:             ${r1}"
echo "R2:             ${r2}"
echo "Cutadapt path:  $(which cutadapt)"
echo "TrimGalore path:$(which trim_galore)"
echo "Start time:     $(date)"
echo "=================================================="

# ------------------------------------------------------------------------------
# Run TrimGalore
# ------------------------------------------------------------------------------
# Key parameters:
#   -q 30         Phred quality score threshold for 3' end trimming
#   --length 35   Discard reads shorter than 35 bp after trimming

trim_galore \
    --fastqc \
    -q 30 \
    -j 4 \
    --length 35 \
    --paired \
    -o "${OUTPUT_DIR}" \
    "${TARGET_DIR}/${r1}" \
    "${TARGET_DIR}/${r2}"

# Log completion
echo "=================================================="
echo "Finished sample: ${base}"
echo "End time:        $(date)"
echo "=================================================="
