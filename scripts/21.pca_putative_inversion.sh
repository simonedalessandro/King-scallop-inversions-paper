#!/bin/sh
#
# Run PCAngsd on putative inversion regions using PLINK files
#
# This script is designed for SGE job arrays on a computer cluster.
# Each array task reads one region name from putative_inversions_list.txt (see file_list)
# and runs PCAngsd on the corresponding PLINK (.bed/.bim/.fam) dataset.
#
# Input:
#   - PLINK files named <REGION>.bed/.bim/.fam under $INPUT
#   - A list of region names (one per line) in putative_inversions_list.txt

# Output:
#   PCA covariance matrices written to $OUTPUT/<REGION>.*
#

# ---------------------- Grid Engine options ---------------------- #
#$ -N local_pca
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=16G
#$ -t 1-14
#$ -pe sharedmem 4
#$ -o output_pmax
#$ -e error_pmax

# ------------------------- Environment --------------------------- #

# Initialise the Modules environment
. /etc/profile.d/modules.sh

# Load conda + PCAngsd environment
module load anaconda
conda activate pcangsd

# ----------------------------- Paths ----------------------------- #

INPUT=/local_path/pmax2023/out/12.inversion/inversion_plink
OUTPUT=/local_path/pmax2023/out/13.MDS_pca_het_LD/local_pca
LGC=/local path/pmax2023/out/file_lists/putative_inversions_list.txt
PCANGSD=/local_path/pcangsd

mkdir -p "$OUTPUT"

# ---------------------- Select region for task -------------------- #

REGION=$(sed -n "${SGE_TASK_ID}p" "$LGC")
PREFIX=${OUTPUT}/${REGION}

echo "[$(date)] Starting PCA for region: ${REGION}"
echo "Input PLINK prefix: ${INPUT}/${REGION}"

# ---------------------------- Run PCAngsd ------------------------- #

$PCANGSD \
    --plink "${INPUT}/${REGION}" \
    --threads 4 \
    --out "${PREFIX}"

echo "[$(date)] Finished PCA for region: ${REGION}"
