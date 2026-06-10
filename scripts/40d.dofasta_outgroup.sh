#!/bin/sh
# ==============================================================================
# Script:       40d_doFasta_outgroup.sh
# Description:  Build an outgroup consensus FASTA from the merged BAM with ANGSD
#               (-doFasta 2 = most-common base; -doCounts 1).
# Called by:    40_polarise.sh (STEP 5)
# Output:       ${OUT_DIR}/outgroup_consensus.fa
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE job configuration
# ------------------------------------------------------------------------------
#$ -N doFasta_outgroup
#$ -cwd
#$ -l h_rt=12:00:00
#$ -l h_vmem=12G
#$ -pe sharedmem 6
#$ -o o_files
#$ -e e_files

# ------------------------------------------------------------------------------
# Modules
# ------------------------------------------------------------------------------
. /etc/profile.d/modules.sh

# ------------------------------------------------------------------------------
# Paths / parameters  — edit to match your layout
# ------------------------------------------------------------------------------
ANGSD=angsd                                              
BAM_DIR=/simone/pmax2023/out/30.ancestral/outgroup/bam
OUT_DIR=/simone/pmax2023/out/30.ancestral/outgroup
THREADS=4

echo "Building consensus on ${HOSTNAME}"

"${ANGSD}" \
    -i "${BAM_DIR}/outgroup_merged.bam" \
    -P "${THREADS}" \
    -doFasta 2 \
    -doCounts 1 \
    -out "${OUT_DIR}/outgroup_consensus"

echo "Done: ${OUT_DIR}/outgroup_consensus.fa"
