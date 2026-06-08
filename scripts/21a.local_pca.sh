#!/bin/bash
# ==============================================================================
# Script:       21a.local_pca.sh
# Description:  Extract each candidate inversion region identified from the
#               lostruct MDS (20b.lostruct.R), and prepare per-region genotype
#               data for local PCA and k-means clustering in 21b.local_pca.R.
#
#               For each candidate region:
#                 1. Extract region from merged VCF using vcftools
#                 2. Export genotype matrix in --012 format
#                 3. Convert to PLINK binary format for downstream analyses
#
# Input:        Pmax_160_merged_sorted_chr.vcf from 20a.prepare_lostruct.sh
#               Candidate regions from 20b.lostruct.R
#
# Output:       Per-region VCF, --012 genotype matrix, and PLINK binary files
#
# Note:         Uses bash associative arrays — run with bash, not sh.
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load vcftools/0.1.16
module load plink/1.9

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

VCF=/simone/pmax2023/out/20.lostruct/Pmax_160_merged_sorted_chr.vcf
OUT_DIR=/simone/pmax2023/out/20.lostruct/inversions
PLINK_DIR=${OUT_DIR}/plink

mkdir -p "${OUT_DIR}" "${PLINK_DIR}"

echo "=================================================="
echo "Start time: $(date)"
echo "VCF input:  ${VCF}"
echo "=================================================="

# ------------------------------------------------------------------------------
# Define candidate inversion regions: CODE -> "SCAFF START END"
# ------------------------------------------------------------------------------
# Scaffold IDs follow the P. maximus reference assembly (GCA_902652985.1).
# Regions were identified as MDS outlier clusters in 20b.lostruct.R.

declare -A REGIONS
REGIONS["LGC15"]="LR736852.1 53 15458540"
REGIONS["LGC02"]="LR736839.1 43890122 54839221"
REGIONS["LGC19.01"]="LR736856.1 2687622 10746375"
REGIONS["LGC17.01"]="LR736854.1 16491332 18526469"
REGIONS["LGC18.02"]="LR736855.1 7997250 21602004"
REGIONS["LGC18.03"]="LR736855.1 21602077 25574752"
REGIONS["LGC18.01"]="LR736855.1 14586 4763448"
REGIONS["LGC08.02"]="LR736845.1 19658887 29776768"
REGIONS["LGC13.01"]="LR736850.1 23743370 26676920"
REGIONS["LGC09"]="LR736846.1 3796103 9574175"
REGIONS["LGC13.02"]="LR736850.1 30396762 32578568"
REGIONS["LGC08.01"]="LR736845.1 4912295 6957646"
REGIONS["LGC11"]="LR736848.1 3301239 4978622"
REGIONS["LGC04"]="LR736841.1 24175137 25307425"
REGIONS["LGC18.04"]="LR736855.1 25575224 26402109"
REGIONS["LGC18.01_02"]="LR736855.1 4226706 6147332"
REGIONS["LGC17.02"]="LR736854.1 1099484 2076514"
REGIONS["LGC01"]="LR736838.1 41411476 42578240"
REGIONS["LGC06"]="LR736843.1 4381139 5156140"
REGIONS["LGC13.03"]="LR736850.1 22528454 23743247"
REGIONS["LGC19.02"]="LR736856.1 11093240 17180908"

# ==============================================================================
# Extract each region, export --012 matrix, and convert to PLINK binary
# ==============================================================================

for CODE in "${!REGIONS[@]}"; do

    read -r SCAFF START END <<< "${REGIONS[$CODE]}"

    echo ""
    echo "[$(date)] Processing: ${CODE} (${SCAFF}:${START}-${END})"

    # Extract region as VCF
    vcftools --vcf "${VCF}" --chr "${SCAFF}" --from-bp "${START}" --to-bp "${END}" \
        --recode --recode-INFO-all --out "${OUT_DIR}/${CODE}"

    # Rename recode output for clarity
    mv "${OUT_DIR}/${CODE}.recode.vcf" "${OUT_DIR}/${CODE}.vcf"

    # Export genotype matrix in --012 format (used for PCA in 21b.local_pca.R)
    vcftools --vcf "${OUT_DIR}/${CODE}.vcf" --012 --out "${OUT_DIR}/${CODE}"

    # Convert to PLINK binary format (used for downstream LD and Fst analyses)
    plink --vcf "${OUT_DIR}/${CODE}.vcf" --double-id --allow-extra-chr \
        --make-bed --out "${PLINK_DIR}/${CODE}"

    echo "[$(date)] Done: ${CODE}"

done

echo ""
echo "[$(date)] All regions extracted."
echo "Outputs: ${OUT_DIR}"
echo ""
echo "Next step:"
echo "  Run 21b.local_pca.R for local PCA and k-means clustering per region."
echo "=================================================="
echo "End time: $(date)"
echo "=================================================="
