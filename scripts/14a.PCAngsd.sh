#!/bin/sh
# ==============================================================================
# Script:       14a.PCAngsd.sh
# Description:  Run PCAngsd on the merged genome-wide beagle and per-chr beagle files
#
# Input:        Per-chromosome Beagle.gz produced by 11.angsd-chr.sh.
#                    Merged Beagle GL file from 11a.merge_ANGSD.sh
#
# Output:    Per-chromosome covariance matrices (.cov) for chromosome-level PCA
#                    Covariance matrix (.cov) for genome-wide PCA
#                   
#
# Note:  Run as an interactive qlogin job
#              Visualisation is performed locally in 15b.pcangsd.R.
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load anaconda

source activate pcangsd

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

ANGSD_DIR=/simone/pmax2023/out/10.angsd
OUT_DIR=/simone/pmax2023/out/12.pcangsd

mkdir -p "${OUTPUT_DIR}"

echo "=================================================="
echo "Working directory: $(pwd)"
echo "Start time:        $(date)"
echo "=================================================="

# ==============================================================================
# PART 1: Genome-wide PCA (merged all-chromosome Beagle file)
# ==============================================================================

echo ""
echo "[$(date)] === Genome-wide PCAngsd ==="

pcangsd \
    -b "${ANGSD_DIR}/Pmax_160_merged_chr.beagle.gz" \
    -o "${OUT_DIR}/Pmax_160_merged_chr" \
    -t 4

echo "[$(date)] Genome-wide covariance matrix written: ${OUT_DIR}/Pmax_160_merged_chr.cov"

# ==============================================================================
# PART 2: Per-chromosome PCA (one Beagle file per chromosome)
# ==============================================================================
# Chromosomes span LR736838.1 (chr 1) to LR736856.1 (chr 19)

echo ""
echo "[$(date)] === Per-chromosome PCAngsd ==="

for chr in {838..856}; do

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
echo "  Download ${OUT_DIR}/*.cov to local machine and run 14b.PCAngsd.R"
echo "  for genome-wide and per-chromosome PCA plots."
echo "=================================================="
echo "End time: $(date)"
echo "=================================================="
