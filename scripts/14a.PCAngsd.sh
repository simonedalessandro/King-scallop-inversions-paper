#!/bin/sh
# ==============================================================================
# Script:       14a.PCAngsd.sh
# Description:  Run PCAngsd on the merged genome-wide beagle (after LD pruning)
#               and on per-chromosome beagle files (unpruned).
#
#               Genome-wide PCA uses LD-pruned SNPs (PLINK --indep-pairwise
#               50 5 0.5, as for the collinear dataset in 31.Prune_filter.sh).
#               Per-chromosome PCAs are deliberately left unpruned: they are used
#               to localise the LD-driven clustering that motivates the local PCA
#               scan (20a-21b), where elevated LD is the diagnostic signal.
#
# Input:        Per-chromosome Beagle.gz produced by 11.angsd-chr.sh.
#               Merged Beagle GL file and PLINK files (.bed/.bim/.fam) from
#               11b.merge_ANGSD.sh
#
# Output:       Pmax_160_LD.prune.in / .prune.out  - PLINK pruning lists
#               Pmax_160_LD_pruned.beagle.gz        - LD-pruned Beagle file
#               Pmax_160_LD_pruned.cov              - genome-wide covariance (Fig. S3)
#               Per-chromosome covariance matrices (.cov) for Fig. 1B
#
# Note:         Run as an interactive qlogin job
#               Visualisation is performed locally in 14b.PCAngsd_plot.R.
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load plink/1.90b7.2
module load anaconda

source activate pcangsd

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

ANGSD_DIR=/anon_user/pmax2023/out/10.angsd
OUT_DIR=/anon_user/pmax2023/out/12.pcangsd

mkdir -p "${OUT_DIR}"

echo "=================================================="
echo "Working directory: $(pwd)"
echo "Start time:        $(date)"
echo "=================================================="

# ==============================================================================
# PART 1: Genome-wide PCA on LD-pruned SNPs
# ==============================================================================

echo ""
echo "[$(date)] === LD pruning (PLINK --indep-pairwise 50 5 0.5) ==="

plink \
    --bfile "${ANGSD_DIR}/Pmax_160_merged_chr" \
    --allow-extra-chr \
    --indep-pairwise 50 5 0.2 \
    --out "${OUT_DIR}/Pmax_160_LD"

# Keep pruned sites in the Beagle file (marker IDs match PLINK IDs: chr_pos)
zcat "${ANGSD_DIR}/Pmax_160_merged_chr.beagle.gz" \
| awk 'NR==FNR {keep[$1]; next} FNR==1 || ($1 in keep)' "${OUT_DIR}/Pmax_160_LD.prune.in" - \
| gzip > "${OUT_DIR}/Pmax_160_LD_pruned.beagle.gz"

echo "[$(date)] Sites kept: $(wc -l < "${OUT_DIR}/Pmax_160_LD.prune.in")"

echo ""
echo "[$(date)] === Genome-wide PCAngsd (LD-pruned) ==="

pcangsd \
    -b "${OUT_DIR}/Pmax_160_LD_pruned.beagle.gz" \
    -o "${OUT_DIR}/Pmax_160_LD_pruned" \
    -t 8 \
    --maf 0

echo "[$(date)] Genome-wide covariance matrix written: ${OUT_DIR}/Pmax_160_LD_pruned.cov"

# ==============================================================================
# PART 2: Per-chromosome PCA (one Beagle file per chromosome, unpruned)
# ==============================================================================
# Chromosomes span LR736838.1 (chr 1) to LR736856.1 (chr 19)

echo ""
echo "[$(date)] === Per-chromosome PCAngsd ==="

for chr in $(seq 838 856); do

    SCAFF="LR736${chr}.1"
    BEAGLE="${ANGSD_DIR}/${SCAFF}.beagle.gz"

    echo "[$(date)] Processing: ${SCAFF}"

    pcangsd \
        -b "${BEAGLE}" \
        -o "${OUT_DIR}/${SCAFF}" \
        -t 4

done

echo ""
echo "[$(date)] Done."
echo "Output directory: ${OUT_DIR}"
echo ""
echo "Next PART:"
echo "  Download ${OUT_DIR}/*.cov to local machine and run 14b.PCAngsd_plot.R"
echo "  for genome-wide and per-chromosome PCA plots."
echo "=================================================="
echo "End time: $(date)"
echo "=================================================="
