#!/bin/sh
# ==============================================================================
# Script:       30.angsd_chr_no_inv.sh
# Description:  Re-run ANGSD per chromosome on the inversion-free dataset
#               after inversion removal (26.remove_inversion.sh).
#               Produces per-chromosome Beagle GL files and
#               PLINK tped/tfam files.
#
#               Run twice (see Note): first on the full 160-sample BAM list,
#               then — after QC outlier removal in 31.pruning.sh — on the
#               158-sample list. The BAM list is passed as a command-line
#               argument so the same script serves both passes.
#
# Usage:        qsub 30.angsd_chr_no_inv.sh <bam_list>
#                 Pass 1: qsub 30.angsd_chr_no_inv.sh \
#                           /simone/pmax2023/out/file_lists/angsd_bams_160.txt
#                 Pass 2: qsub 30.angsd_chr_no_inv.sh \
#                           /simone/pmax2023/out/file_lists/angsd_bams_158.txt
#
# Input:        BAM list (passed as $1)
#               Pmax_no_inversion.snplist — inversion-free SNP position list
#               from 26.remove_inversion.sh (PLINK --write-snplist), reformatted
#               to tab-separated "scaffold <TAB> position" and indexed.
#
# Output:       Per-chromosome Beagle GL files and PLINK tped/tfam files
#               in OUT_DIR.
#
# Note:         ANGSD -sites requires a binary index of the SNP list. Generate
#               once before running:
#                 angsd sites index Pmax_no_inversion.snplist
# Next step:    31.pruning.sh
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE array job configuration
# ------------------------------------------------------------------------------
#$ -N angsd_chr_no_inv
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=16G
#$ -t 1-19
#$ -pe sharedmem 4
#$ -o o_files
#$ -e e_files

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load java

# ------------------------------------------------------------------------------
# Command-line argument: BAM list
# ------------------------------------------------------------------------------

BAM_LIST="${1:?ERROR: provide a BAM list, e.g. qsub 30.angsd_chr_no_inv.sh /simone/pmax2023/out/file_lists/angsd_bams_160.txt}"

if [ ! -f "${BAM_LIST}" ]; then
    echo "ERROR: BAM list not found: ${BAM_LIST}" >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

ANGSD=/simone/pmax2023/software/angsd-v0.940/angsd

OUT_DIR=/simone/pmax2023/out/22.no_inversion
REFERENCE=/simone/reference/Pmax/GCA_902652985.1_xPecMax1.1_genomic.fna
SCAFF_LIST=/simone/pmax2023/out/file_lists/scaff_IDs.txt
SITES=/simone/pmax2023/out/file_lists/Pmax_no_inversion.snplist

mkdir -p "${OUT_DIR}"

# ------------------------------------------------------------------------------
# Parse chromosome for this array task
# ------------------------------------------------------------------------------

SCAFF=$(sed -n "${SGE_TASK_ID}p" "${SCAFF_LIST}" | awk '{print $1}')

echo "=================================================="
echo "Task ID:     ${SGE_TASK_ID}"
echo "Chromosome:  ${SCAFF}"
echo "BAM list:    ${BAM_LIST}"
echo "Sites file:  ${SITES}"
echo "Host:        ${HOSTNAME}"
echo "Start time:  $(date)"
echo "=================================================="

# ------------------------------------------------------------------------------
# Run ANGSD
# ------------------------------------------------------------------------------

"${ANGSD}" \
    -b          "${BAM_LIST}" \
    -ref        "${REFERENCE}" \
    -r          "${SCAFF}" \
    -sites      "${SITES}" \
    -out        "${OUT_DIR}/${SCAFF}" \
    -uniqueOnly 1 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -trim 0 \
    -C 50 \
    -baq 1 \
    -minMapQ 30 \
    -minQ 30 \
    -doCounts 1 \
    -GL 2 \
    -doGlf 2 \
    -doMajorMinor 4 \
    -doMaf 1 \
    -rmTriallelic 1e-6 \
    -SNP_pval 1e-6 \
    -doGeno 2 \
    -doPost 1 \
    -doPlink 2 \
    -P 4

echo "=================================================="
echo "Finished chromosome: ${SCAFF}"
echo "End time:            $(date)"
echo "=================================================="
