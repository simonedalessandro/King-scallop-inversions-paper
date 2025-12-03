#!/bin/sh
#
# het.sh — Compute observed heterozygosity (het) for putative inversion regions
#           using VCFtools.
#
# This script is designed for SGE job arrays on a computer cluster.
# Each array task reads one region name from putative_inversions_list.txt (see file_list)
# and runs `vcftools --het` on the corresponding VCF file.
#
#Input:
#   - VCF files named <REGION>.vcf under $INPUT
#   - A list of region names (one per line) in putative_inversions_list.txt
#
# Output:
#   For each region, a *.<REGION>.het* file written under $OUTPUT
#

# ---------------------- Grid Engine options ---------------------- #
#$ -N het
#$ -cwd
#$ -l h_rt=6:00:00
#$ -l h_vmem=16G
#$ -t 1-14
#$ -pe sharedmem 4
#$ -o output_pmax
#$ -e error_pmax

# ------------------------- Environment --------------------------- #

# Initialise the Modules environment
. /etc/profile.d/modules.sh

# Load VCFtools
module load vcftools/0.1.16

# ----------------------------- Paths ----------------------------- #

INPUT=/local_path/pmax2023/out/12.inversion/inversion_vcf
OUTPUT=/local_path/pmax2023/out/13.MDS_pca_het_LD/het
LGC=/local_path/pmax2023/out/file_lists/putative_inversions_list.txt

mkdir -p "$OUTPUT"

# ---------------------- Select region for task -------------------- #

REGION=$(sed -n "${SGE_TASK_ID}p" "$LGC")
PREFIX=${OUTPUT}/${REGION}

echo "[$(date)] Starting heterozygosity calculation for region: ${REGION}"
echo "Using VCF: ${INPUT}/${REGION}.vcf"

# ---------------------------- Run VCFtools ------------------------ #

vcftools \
    --vcf "${INPUT}/${REGION}.vcf" \
    --het \
    --out "${PREFIX}"

echo "[$(date)] Finished heterozygosity for region: ${REGION}"
