#!/bin/sh
# ==============================================================================
# Script:       40b_SortSam_outgroup.sh
# Description:  Coordinate-sort and index the mapped outgroup BAMs.
# Called by:    40_polarise.sh (STEP 3)
# Output:       ${BAM_DIR}/<SRA>.sorted.bam (+ .bai)
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE job configuration
# ------------------------------------------------------------------------------
#$ -N sort_outgroup
#$ -cwd
#$ -l h_rt=16:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -t 1-15
#$ -o o_files
#$ -e e_files

# ------------------------------------------------------------------------------
# Modules
# ------------------------------------------------------------------------------
. /etc/profile.d/modules.sh

module load samtools/1.10

# ------------------------------------------------------------------------------
# Paths / parameters
# ------------------------------------------------------------------------------
BAM_DIR=/simone/pmax2023/out/30.ancestral/outgroup/bam
THREADS=4

# Pick this task's BAM
this_bam=$(ls -1 "${BAM_DIR}"/*.mapped.bam | sed -n "${SGE_TASK_ID}p")
out_bam="${this_bam%.mapped.bam}.sorted.bam"

echo "Sorting ${this_bam} on ${HOSTNAME}"

samtools sort -@ "${THREADS}" -o "${out_bam}" "${this_bam}"
samtools index "${out_bam}"

echo "Done: ${out_bam}"
