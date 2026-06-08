#!/bin/sh
# ==============================================================================
# Script:       23.LD.sh
# Description:  Compute pairwise LD across full chromosomes associated with
#               putative inversions, using PLINK and VCFtools.
#
#               For each (CHR, LGC) pair in chrom_inversion_list.txt:
#                 1. SNP filtering (MAF > 0.05)
#                 2. VCF reheadering
#                 3. Subset to common-homozygous (AA) individuals for each inversion
#                 4. LD calculation within 100 kb for:
#                      - all individuals
#                      - AA-only individuals
#
# Input:        Chromosome-wide PLINK files from ANGSD -doPlink (11a.angsd-chr.sh)
#               Per-inversion AA sample lists: <LGC>_common_hom_individuals.txt
#               chrom_inversion_list.txt — two columns: CHR  LGC
#
# Output:       LD/<LGC>.ld     — LD for all individuals
#               LD/<LGC>.AA.ld  — LD for common-homozygous (AA) individuals
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE array job configuration
# ------------------------------------------------------------------------------
#$ -N LD
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_rss=30G
#$ -t 1-14
#$ -pe sharedmem 4
#$ -o o_files
#$ -e e_files

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load vcftools/0.1.16
module load plink/1.9
module load bcftools/1.20

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

PLINK_DIR=/simone/pmax2023/out/10.angsd/Pmax_160
LD_DIR=/simone/pmax2023/out/21.MDS_pca_het_LD/LD
TMP_DIR=${LD_DIR}/tmp
MAPFILE=/simone/pmax2023/out/file_lists/chrom_inversion_list.txt
REHEADER=/simone/pmax2023/out/file_lists/bcftools_reheader.txt

mkdir -p "${LD_DIR}" "${TMP_DIR}"

# ------------------------------------------------------------------------------
# Parse chromosome and inversion ID for this array task
# ------------------------------------------------------------------------------

LINE=$(sed -n "${SGE_TASK_ID}p" "${MAPFILE}")
CHR=$(echo "${LINE}" | awk '{print $1}')
LGC=$(echo "${LINE}" | awk '{print $2}')

AA_LIST=${LD_DIR}/${LGC}_common_hom_individuals.txt
TMP_PREFIX=${TMP_DIR}/${LGC}
OUT_PREFIX=${LD_DIR}/${LGC}

echo "=================================================="
echo "Task ID:     ${SGE_TASK_ID}"
echo "Chromosome:  ${CHR}"
echo "Inversion:   ${LGC}"
echo "Start time:  $(date)"
echo "=================================================="

# ==============================================================================
# STEP 1: SNP filtering and VCF conversion
# ==============================================================================
# Filter SNPs by MAF > 0.05, then convert to VCF for reheadering.

plink --bfile "${PLINK_DIR}/${CHR}" --allow-extra-chr --recode \
    --maf 0.05 --out "${TMP_PREFIX}.step1"

plink --file "${TMP_PREFIX}.step1" --allow-extra-chr --recode vcf \
    --make-bed --out "${TMP_PREFIX}.step2"

# ==============================================================================
# STEP 2: Reheader VCF with sample names
# ==============================================================================

bcftools reheader \
    -s "${REHEADER}" \
    "${TMP_PREFIX}.step2.vcf" > "${TMP_PREFIX}.step3.vcf"

# ==============================================================================
# STEP 3: Subset to common-homozygous (AA) individuals
# ==============================================================================

vcftools --vcf "${TMP_PREFIX}.step3.vcf" --keep "${AA_LIST}" \
    --recode --out "${TMP_PREFIX}.AA"

# ==============================================================================
# STEP 4: Convert VCF to PLINK binary (all individuals and AA-only)
# ==============================================================================

plink --vcf "${TMP_PREFIX}.step3.vcf" --make-bed --allow-extra-chr \
    --const-fid 0 --out "${TMP_PREFIX}.all"

plink --vcf "${TMP_PREFIX}.AA.recode.vcf" --make-bed --allow-extra-chr \
    --const-fid 0 --out "${TMP_PREFIX}.AA"

# ==============================================================================
# STEP 5: LD calculation within 100 kb
# ==============================================================================
# --ld-window-r2 0 retains all SNP pairs (no R2 threshold) for full heatmap

plink --bfile "${TMP_PREFIX}.all" --r2 --ld-window-kb 100000 --ld-window-r2 0 \
    --allow-extra-chr --out "${OUT_PREFIX}"

plink --bfile "${TMP_PREFIX}.AA" --r2 --ld-window-kb 100000 --ld-window-r2 0 \
    --allow-extra-chr --out "${OUT_PREFIX}.AA"

# ==============================================================================
# STEP 6: Remove intermediate files
# ==============================================================================

rm -f "${TMP_PREFIX}.step1"*
rm -f "${TMP_PREFIX}.step2"*
rm -f "${TMP_PREFIX}.step3"*
rm -f "${TMP_PREFIX}.all"*
rm -f "${TMP_PREFIX}.AA."*

echo "=================================================="
echo "Finished inversion: ${LGC} (chromosome ${CHR})"
echo "End time:           $(date)"
echo "=================================================="
