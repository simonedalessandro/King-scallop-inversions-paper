#!/bin/sh
# ==============================================================================
# Script:       40c_MergeSam_outgroup.sh
# Description:  Merge all sorted outgroup BAMs into a single file
# Called by:    40_polarise.sh (STEP 4)
# Output:       ${BAM_DIR}/outgroup_merged.bam (+ .bai)
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE job configuration
# ------------------------------------------------------------------------------
#$ -N merge_outgroup
#$ -cwd
#$ -l h_rt=6:00:00
#$ -l h_vmem=12G
#$ -pe sharedmem 4
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

echo "Merging sorted BAMs in ${BAM_DIR} on ${HOSTNAME}"

samtools merge -f -@ "${THREADS}" \
    "${BAM_DIR}/outgroup_merged.bam" \
    "${BAM_DIR}"/*.sorted.bam

samtools index "${BAM_DIR}/outgroup_merged.bam"

echo "Done: ${BAM_DIR}/outgroup_merged.bam"
