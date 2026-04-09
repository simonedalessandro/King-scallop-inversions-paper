#!/bin/sh
# ==============================================================================
# Script:       11a.angsd-chr.sh
# Description:  Compute genotype likelihoods (GLs) and allele frequencies per
#               chromosome using ANGSD.
#
# Input:        List of deduplicated BAM files (one path per line)
#                    List of scaffold/chromosome IDs (one ID per line)
# Output:     Per-chromosome Beagle GL files, MAF, and PLINK tped/tfam files
#
# Note:
# This script is run twice:
#
#    Step 1. full 168-sample dataset used for initial individual-level QC (12a/b.het_miss) and
#                 pairwise relatedness estimation (13a/b.relatedness).
#                 8 individuals identified as related are flagged for removal.
#
#      Step 2. cleaned 160-sample dataset re-run after removing the 8 flagged individuals 
#                 from the BAM list and re-running ANGSD (11.angsd-chr.sh). 
#
# ANGSD depth filters scale with sample count and -setMinDepth, -minInd and -setMaxDepth
# need to be updated accordingly
#
# Next steps:   11b.merge_ANGSD.sh ==============================================================================

# ------------------------------------------------------------------------------
# SGE array job configuration
# ------------------------------------------------------------------------------
#$ -N angsd_chr
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_rss=16G
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
# Path configuration
# ------------------------------------------------------------------------------

ANGSD=/software/angsd/angsd/angsd

# Update BAM_LIST between passes:
#   Pass 1: angsd_bams_168.txt  (full cohort)
#   Pass 2: angsd_bams_160.txt  (8 related individuals removed)
BAM_LIST=/simone/pmax2023/out/file_lists/angsd_bams_160.txt

OUTPUT_DIR=/simone/pmax2023/out/10.angsd
REFERENCE=/simone/reference/Pmax/GCA_902652985.1_xPecMax1.1_genomic.fna
SCAFF_LIST=/simone/pmax2023/out/file_lists/scaff_IDs.txt

# ------------------------------------------------------------------------------
# Parse chromosome for this array task
# ------------------------------------------------------------------------------

SCAFF=$(sed -n "${SGE_TASK_ID}p" "${SCAFF_LIST}" | awk '{print $1}')

echo "=================================================="
echo "Task ID:     ${SGE_TASK_ID}"
echo "Chromosome:  ${SCAFF}"
echo "BAM list:    ${BAM_LIST}"
echo "Host:        ${HOSTNAME}"
echo "Start time:  $(date)"
echo "=================================================="

# ------------------------------------------------------------------------------
# Run ANGSD
# ------------------------------------------------------------------------------
# Key parameters:
#   -GL 2                  GATK genotype likelihood model
#   -doGlf 2               Output Beagle format GLs
#   -doMajorMinor 4        Infer major/minor allele from reference
#   -doMaf 1               Estimate allele frequencies
#   -doPlink 2             Output tped/tfam for PLINK
#   -doGeno 2 / -doPost 1  Genotype calling using posterior probabilities
#   -rmTriallelic 1e-6     Remove triallelic sites
#   -SNP_pval 1e-6         Retain only statistically supported SNPs
#   -minMaf 0.05           Minor allele frequency filter
#   -minInd 128            Require data in >= 128 of 160 individuals (~80%)
#   -setMinDepth 480       Minimum global depth (0.5 * 160 * 6x)
#   -setMaxDepth 1920      Maximum global depth (2 * 160 * 6x)

"${ANGSD}" \
    -b "${BAM_LIST}" \
    -ref "${REFERENCE}" \
    -r "${SCAFF}" \
    -out "${OUTPUT_DIR}/${SCAFF}" \
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
    -P 4 \
    -setMinDepth 480 \
    -setMaxDepth 1920 \
    -minInd 128 \
    -minMaf 0.05

echo "=================================================="
echo "Finished chromosome: ${SCAFF}"
echo "End time:            $(date)"
echo "=================================================="
