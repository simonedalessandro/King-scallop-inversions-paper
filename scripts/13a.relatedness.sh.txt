#!/bin/sh
# ==============================================================================
# Script:       13a.relatedness.sh
# Description:  Estimate pairwise relatedness using PLINK (IBD / PI_HAT) and ngsRelate (KING, R0, R1)
#
#               Pipeline:
#                 1. LD pruning of the merged ANGSD SNP dataset
#                 2. Additional SNP filtering (missingness, MAF, HWE)
#                 3. PLINK --genome to compute pairwise IBD (PI_HAT, Z0, Z1, Z2)
#                 4. ngsRelate to compute KING, R0, R1 from genotype calls
#
# Note:         Download .genome and .res to local machine and run
#               13b.relatedness.R for classification and visualisation.
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE job configuration
# ------------------------------------------------------------------------------
#$ -N relate
#$ -cwd
#$ -l h_rt=10:00:00
#$ -l h_rss=16G
#$ -pe sharedmem 1
#$ -o o_files
#$ -e e_files

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

set -euo pipefail

. /etc/profile.d/modules.sh

module /plink/1.90b7.2

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

ANGSD_DIR=/simone/pmax2023/out/10.angsd
OUT_DIR=/simone/pmax2023/out/11.relatedness

PREFIX=${ANGSD_DIR}/Pmax_168_merged_chr
BEAGLE=${ANGSD_DIR}/Pmax_168_merged_chr.beagle.gz
NGSRELATE=/software/ngsRelate/ngsRelate

mkdir -p "${OUT_DIR}"

echo "[$(date)] Host:         ${HOSTNAME}"
echo "[$(date)] Input prefix: ${PREFIX}"
echo "[$(date)] Beagle file:  ${BEAGLE}"
echo "[$(date)] Output dir:   ${OUT_DIR}"

# ==============================================================================
# STEP 1–2: LD pruning and SNP filtering
# ==============================================================================

echo ""
echo "[$(date)] === PLINK IBD ==="

# ------------------------------------------------------------------------------
# Step 1: LD pruning
# Generates prune.in (retained SNPs) and prune.out (removed SNPs)
# ------------------------------------------------------------------------------

plink \
    --bfile "${PREFIX}" \
    --allow-extra-chr \
    --make-founders \
    --indep-pairwise 50 5 0.5 \
    --out "${OUT_DIR}/Pmax_168_merged_chr"

echo "[$(date)] SNPs before pruning : $(wc -l < ${PREFIX}.bim)"
echo "[$(date)] SNPs in prune.out   : $(wc -l < ${OUT_DIR}/Pmax_168_merged_chr.prune.out)"
echo "[$(date)] SNPs in prune.in    : $(wc -l < ${OUT_DIR}/Pmax_168_merged_chr.prune.in)"

# ------------------------------------------------------------------------------
# Step 2: Apply pruning and additional SNP filters
# Filters applied:
#   --geno 0.01   exclude SNPs with > 1% missing genotypes
#   --maf 0.3     exclude SNPs with minor allele frequency < 0.3
#   --hwe 0.001   exclude SNPs deviating from HWE (p < 0.001)
# ------------------------------------------------------------------------------

plink \
    --bfile "${PREFIX}" \
    --allow-extra-chr \
    --make-founders \
    --geno 0.01 \
    --maf 0.3 \
    --hwe 0.001 \
    --extract "${OUT_DIR}/Pmax_168_merged_chr.prune.in" \
    --make-bed \
    --out "${OUT_DIR}/Pmax_168_merged_pruned"

echo "[$(date)] SNPs after pruning + filtering: $(wc -l < ${OUT_DIR}/Pmax_168_merged_pruned.bim)"

# ==============================================================================
# STEP 3: PLINK pairwise IBD (PI_HAT)
# Computes Z0, Z1, Z2 (probabilities of sharing 0, 1, 2 alleles IBD) and
# PI_HAT (overall proportion of alleles shared IBD) for all pairs.
# ==============================================================================

plink \
    --bfile "${OUT_DIR}/Pmax_168_merged_pruned" \
    --allow-extra-chr \
    --genome \
    --out "${OUT_DIR}/Pmax_168_merged_pruned"

# ==============================================================================
# STEP 4: ngsRelate 
# ==============================================================================

echo ""
echo "[$(date)] === ngsRelate ==="

# Extract individual IDs from pruned .fam file (column 2 = IID)
awk '{print $2}' "${OUT_DIR}/Pmax_168_merged_pruned.fam" \
    > "${OUT_DIR}/sample_ids.txt"

N_SAMPLES=$(wc -l < "${OUT_DIR}/sample_ids.txt")
echo "[$(date)] Sample count: ${N_SAMPLES}"

# Run ngsRelate on pruned genotype calls
"${NGSRELATE}" \
    -P "${OUT_DIR}/Pmax_168_merged_pruned" \
    -T GT \
    -O "${OUT_DIR}/Pmax_168_merged_pruned.res" \
    -c 1

echo ""
echo "[$(date)] Done."
echo "PLINK PI_HAT    : ${OUT_DIR}/Pmax_168_merged_pruned.genome"
echo "ngsRelate output: ${OUT_DIR}/Pmax_168_merged_pruned.res"
