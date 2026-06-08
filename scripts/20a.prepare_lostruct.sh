#!/bin/sh
# ==============================================================================
# Script:       20a.prepare_lostruct.sh
# Description:  Prepare a merged, sorted, indexed VCF/BCF file for use with
#               lostruct. Converts the per-chromosome binary PLINK files
#               produced by 11b.merge_ANGSD.sh to VCF, merges all 19
#               chromosomes, and produces an indexed BCF and --012 genotype
#               matrix.
#
# Input:        Per-chromosome binary PLINK files (LR736838.1–LR736856.1)
#
# Output:       merged_sorted_chr.bcf  — indexed BCF for lostruct
#               merged_sorted_chr.012  — genotype matrix for global PCA
#
# Note:         Run as an interactive qlogin job — no SGE submission required.
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load vcftools/0.1.16
module load bcftools/1.20
module load plink/1.09
module load samtools/1.9

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

ANGSD_DIR=/simone/pmax2023/out/10.angsd
VCF_DIR=/simone/pmax2023/out/20.lostruct/vcf
OUT_DIR=/simone/pmax2023/out/20.lostruct

mkdir -p "${VCF_DIR}" "${OUT_DIR}"

echo "=================================================="
echo "Working directory: $(pwd)"
echo "Start time:        $(date)"
echo "=================================================="

# ==============================================================================
# STEP 1: Recode per-chromosome binary PLINK files to VCF
# ==============================================================================
# Binary PLINK files (.bed/.bim/.fam) were produced by 11b.merge_ANGSD.sh.
# Chromosomes span LR736838.1 (chr 1) to LR736856.1 (chr 19).

cd "${ANGSD_DIR}"

for chr in {838..856}; do

    BFILE="LR736${chr}.1"
    echo "[$(date)] Recoding: ${BFILE}"

    plink --bfile "${BFILE}" --recode vcf --allow-extra-chr --out "${VCF_DIR}/${BFILE}"

done

echo "[$(date)] All chromosomes recoded to VCF"

# ==============================================================================
# STEP 2: Merge all per-chromosome VCFs into a single file
# ==============================================================================

cd "${VCF_DIR}"

VCF_LIST=$(printf '%s ' LR736{838..856}.1.vcf)

bcftools concat -o "${OUT_DIR}/Pmax_160_merged_chr.vcf" ${VCF_LIST}

echo "[$(date)] Merged VCF written: ${OUT_DIR}/Pmax_160_merged_chr.vcf"

# ==============================================================================
# STEP 3: Sort, compress, and index merged VCF
# ==============================================================================

cd "${OUT_DIR}"

# Retain header lines first, then sort variant lines by chromosome then position
(grep "^#" Pmax_160_merged_chr.vcf; grep -v "^#" Pmax_160_merged_chr.vcf | sort -k1,1 -k2,2n) \
    > Pmax_160_merged_sorted_chr.vcf

echo "[$(date)] VCF sorted"

# Compress with bgzip (block gzip) and index with tabix for random access
bgzip -c Pmax_160_merged_sorted_chr.vcf > Pmax_160_merged_sorted_chr.vcf.gz
tabix -fp vcf Pmax_160_merged_sorted_chr.vcf.gz

echo "[$(date)] VCF compressed and indexed: Pmax_160_merged_sorted_chr.vcf.gz"

# ==============================================================================
# STEP 4: Convert to BCF and index
# ==============================================================================
# lostruct's vcf_windower requires an indexed BCF file as input.

bcftools convert -O b Pmax_160_merged_sorted_chr.vcf.gz > Pmax_160_merged_sorted_chr.bcf
bcftools index Pmax_160_merged_sorted_chr.bcf

echo "[$(date)] BCF written and indexed:"
ls -lh Pmax_160_merged_sorted_chr.bcf

# ==============================================================================
# STEP 5: Export genotype matrix in --012 format
# ==============================================================================
# The --012 matrix is used for the global PCA correlation step in 20b.lostruct.R.

vcftools --bcf  Pmax_160_merged_sorted_chr.bcf --012 --out  Pmax_160_merged_sorted_chr

echo "[$(date)] Genotype matrix written: Pmax_160_merged_sorted_chr.012"

echo ""
echo "[$(date)] Done."
echo "Outputs:"
echo "  BCF (lostruct input):  ${OUT_DIR}/Pmax_160_merged_sorted_chr.bcf"
echo "  012 matrix (PCA):      ${OUT_DIR}/Pmax_160_merged_sorted_chr.012"
echo ""
echo "Next step:"
echo "  Run 20b.lostruct.R to perform windowed PCA and MDS outlier detection."
echo "=================================================="
echo "End time: $(date)"
echo "=================================================="
