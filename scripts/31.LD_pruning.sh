#!/bin/bash
# ==============================================================================
# Script:       31.LD_pruning.sh
# Description:  LD pruning of the inversion-free ANGSD dataset, followed by
#               individual QC (missingness and heterozygosity), relatedness
#               check, and preparation of a final 158-sample dataset for
#               downstream population genomic analyses.
#
#               Pipeline:
#                 1. Convert per-chromosome tped/tfam to binary PLINK format
#                 2. Merge all chromosomes into a single genome-wide PLINK file
#                 3. LD pruning (--indep-pairwise 50 5 0.5)
#                 4. Apply pruning and compute missingness and heterozygosity
#                 5. Flag het outliers (> 3 SD) and high-missingness individuals
#                 6. Relatedness check (PLINK --genome + ngsRelate)
#                 7. Remove confirmed outliers (WCO23_06_17; WCO23_06_04)
#                    and re-run ANGSD with the updated 158-sample BAM list
#                 8. Merge per-chromosome Beagle files into genome-wide file
#                 9. Run PCAngsd on the final 158-sample dataset
#
# Input:        Per-chromosome tped/tfam files from 30.angsd_chr_no_inv.sh
#
# Output:       Pmax_no_inv_pruned.bed/.bim/.fam — final PLINK dataset
#               Pmax_no_inv_pruned.beagle.gz     — merged Beagle GL file
#               Pmax_no_inv_pruned.cov           — PCAngsd covariance matrix
#
# Note:         Run as an interactive qlogin job — no SGE submission required.
#               Uses bash brace expansion — run with bash, not sh.
#               Outlier removal (step 7) requires re-running 30.angsd_chr_no_inv.sh
#               with the 158-sample BAM list before continuing from step 8.
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load plink/1.90b7.2
module load anaconda

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

ANGSD_DIR=/simone/pmax2023/out/22.no_inversion
OUT_DIR=/simone/pmax2023/out/22.no_inversion
PRUNE_DIR=${OUT_DIR}/pruned
NGSRELATE=/software/ngsRelate/ngsRelate

mkdir -p "${PRUNE_DIR}"

echo "=================================================="
echo "Start time:        $(date)"
echo "Working directory: ${ANGSD_DIR}"
echo "=================================================="

cd "${ANGSD_DIR}"

# ==============================================================================
# STEP 1: Convert per-chromosome tped/tfam to binary PLINK format
# ==============================================================================
# Chromosomes span LR736838.1 (chr 1) to LR736856.1 (chr 19).

echo ""
echo "[$(date)] === Converting tped/tfam to binary PLINK format ==="

for chr in {838..856}; do
    TFILE="LR736${chr}.1"
    echo "[$(date)] Converting: ${TFILE}"
    plink \
        --tfile       "${TFILE}" \
        --allow-extra-chr \
        --make-bed \
        --out         "${TFILE}"
done

echo "[$(date)] All chromosomes converted"

# ==============================================================================
# STEP 2: Merge all chromosomes into a single genome-wide PLINK file
# ==============================================================================

echo ""
echo "[$(date)] === Merging chromosomes ==="

printf '%s\n' LR736{839..856}.1 > merge_list.txt

plink \
    --bfile       LR736838.1 \
    --merge-list  merge_list.txt \
    --allow-extra-chr \
    --make-bed \
    --out         Pmax_no_inv

echo "[$(date)] Merged dataset: $(wc -l < Pmax_no_inv.bim) variants, $(wc -l < Pmax_no_inv.fam) individuals"

# ==============================================================================
# STEP 3: LD pruning
# ==============================================================================

echo ""
echo "[$(date)] === LD pruning ==="

plink \
    --bfile       Pmax_no_inv \
    --allow-extra-chr \
    --indep-pairwise 50 5 0.5 \
    --recode \
    --make-bed \
    --out         "${PRUNE_DIR}/Pmax_no_inv_LD"

echo "[$(date)] SNPs before pruning: $(wc -l < Pmax_no_inv.bim)"
echo "[$(date)] SNPs removed:        $(wc -l < ${PRUNE_DIR}/Pmax_no_inv_LD.prune.out)"
echo "[$(date)] SNPs retained:       $(wc -l < ${PRUNE_DIR}/Pmax_no_inv_LD.prune.in)"

# ==============================================================================
# STEP 4: Apply pruning and compute missingness and heterozygosity
# ==============================================================================

echo ""
echo "[$(date)] === Applying pruning and computing QC statistics ==="

cd "${PRUNE_DIR}"

plink \
    --file        Pmax_no_inv_LD \
    --allow-extra-chr \
    --exclude     Pmax_no_inv_LD.prune.out \
    --make-bed \
    --recode \
    --missing \
    --het \
    --out         Pmax_no_inv_pruned

echo "[$(date)] Pruned dataset: $(wc -l < Pmax_no_inv_pruned.bim) variants"

# ==============================================================================
# STEP 5: Flag QC outliers
# ==============================================================================
# Individuals with > 10% missing genotypes are flagged.
# Het outliers are defined as > 3 SD from the mean F statistic.

echo ""
echo "[$(date)] === Flagging QC outliers ==="

# High missingness (> 10%)
awk 'NR == 1 || $6 > 0.10' Pmax_no_inv_pruned.imiss > high_missing_pruned.txt
echo "[$(date)] High-missingness individuals: $(( $(wc -l < high_missing_pruned.txt) - 1 ))"

# Heterozygosity outliers (> 3 SD from mean F)
awk 'NR == 1 {print $0; next}
     NR > 1  {sum += $6; sumsq += $6*$6; vals[NR] = $6; lines[NR] = $0}
     END {
       mean = sum / (NR - 1)
       sd   = sqrt(sumsq / (NR - 1) - mean * mean)
       for (i in vals) {
         if (vals[i] < (mean - 3*sd) || vals[i] > (mean + 3*sd))
           print lines[i]
       }
     }' Pmax_no_inv_pruned.het > het_outliers_pruned.txt

echo "[$(date)] Het outliers: $(( $(wc -l < het_outliers_pruned.txt) - 1 ))"

# ==============================================================================
# STEP 6: Relatedness check
# ==============================================================================

echo ""
echo "[$(date)] === Relatedness check ==="

plink \
    --file        Pmax_no_inv_pruned \
    --geno 0.01 \
    --maf 0.3 \
    --hwe 0.001 \
    --allow-extra-chr \
    --genome \
    --make-bed \
    --recode \
    --missing \
    --out         Pmax_no_inv_pruned_relate

"${NGSRELATE}" \
    -P Pmax_no_inv_pruned_relate \
    -T GT \
    -O Pmax_no_inv_pruned_relate.res \
    -c 1

echo "[$(date)] Relatedness output: Pmax_no_inv_pruned_relate.res"

# ==============================================================================
# STEP 7: Remove outliers and re-run ANGSD
# ==============================================================================
# WCO23_06_17 (het outlier) and WCO23_06_04 are removed.
# Update the BAM list to exclude these two individuals, then re-run
# 30.angsd_chr_no_inv.sh with the 158-sample list (passed as an argument).
# Continue from STEP 8 once the new ANGSD run is complete.

echo ""
echo "[$(date)] === NOTE: Remove WCO23_06_17 and WCO23_06_04 from BAM list ==="
echo "[$(date)]   qsub 30.angsd_chr_no_inv.sh /simone/pmax2023/out/file_lists/angsd_bams_158.txt"
echo "[$(date)]   then continue from STEP 8."

# ==============================================================================
# STEP 8: Merge per-chromosome Beagle files (158-sample dataset)
# ==============================================================================
# Concatenate all 19 per-chromosome Beagle GL files, keeping the header from
# the first file only. zcat reads the .gz files without altering them; the
# merged file is line-counted before being compressed.

echo ""
echo "[$(date)] === Merging Beagle GL files (158 samples) ==="

cd "${ANGSD_DIR}"

MERGED=Pmax_no_inv_pruned.beagle
: > "${MERGED}"
first=1

for chr in {838..856}; do
    BEAGLE="LR736${chr}.1.beagle.gz"
    if [ ! -f "${BEAGLE}" ]; then
        echo "[$(date)] WARNING: missing ${BEAGLE}, skipping" >&2
        continue
    fi
    if [ "${first}" -eq 1 ]; then
        zcat "${BEAGLE}" >> "${MERGED}"            # keep header
        first=0
    else
        zcat "${BEAGLE}" | tail -n +2 >> "${MERGED}"   # drop header
    fi
done

echo "[$(date)] SNP count (excl. header): $(( $(wc -l < "${MERGED}") - 1 ))"

gzip -f "${MERGED}"

echo "[$(date)] Merged Beagle written: ${MERGED}.gz"

# ==============================================================================
# STEP 9: PCAngsd on the final 158-sample inversion-free dataset
# ==============================================================================

echo ""
echo "[$(date)] === Running PCAngsd ==="

source activate pcangsd

pcangsd \
    -b  Pmax_no_inv_pruned.beagle.gz \
    -o  Pmax_no_inv_pruned \
    -t  4 \
    --maf 0

echo ""
echo "[$(date)] Done."
echo "Outputs:"
echo "  PLINK dataset:   ${PRUNE_DIR}/Pmax_no_inv_pruned.bed/.bim/.fam"
echo "  Beagle GL file:  ${ANGSD_DIR}/Pmax_no_inv_pruned.beagle.gz"
echo "  PCAngsd cov:     ${ANGSD_DIR}/Pmax_no_inv_pruned.cov"
echo "=================================================="
echo "End time: $(date)"
echo "=================================================="
