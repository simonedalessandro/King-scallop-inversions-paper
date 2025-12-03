#!/bin/sh
#
# Compute pairwise LD for full chromosomes associated
# with putative inversions, using PLINK and VCFtools.
#
# This script is designed for SGE job arrays on a computer cluster.
# Each array task reads one (CHR, LGC) pair from chrom_inversion_list.txt (see file_list)
# and performs:
#    1. SNP filtering (MAF > 0.05)
#    2. VCF reheadering
#    3. Subsetting to the "common homozygous group" individuals (i.e., AA) for the inversion
#       for subsequent plot
#    4. LD calculation within 100 kb for:
#         - all individuals
#         - common hom individuals
#
# Inputs:
#   - Chromosome-wide PLINK files from angsd -doPlink in $PLINK
#   - Per-inversion AA lists: <LGC>_common_hom_individuals.txt generated during
#     [INSERT PIPELINE STEP]
#   - chrom_inversion_list.txt
#
# Output:
#   LD/<LGC>.ld       — LD for all individuals
#   LD/<LGC>.AA.ld    — LD for common-hom individuals
#

# ---------------------- Grid Engine options ---------------------- #
#$ -N LD
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=30G
#$ -t 1-14
#$ -pe sharedmem 4
#$ -o output_pmax
#$ -e error_pmax

# ------------------------- Environment --------------------------- #

# Initialise the Modules environment
. /etc/profile.d/modules.sh

# Load necessary software
module load vcftools/0.1.16
module load plink/1.90b7.2
module load bcftools/1.20

# ----------------------------- Paths ----------------------------- #

BASE=/local_path/pmax2023/out
PLINK=${BASE}/10.angsd/rm8_160_06x
LD=${BASE}/13.MDS_pca_het_LD/LD
TMP=${LD}/tmp
MAPFILE=${BASE}/file_lists/chrom_inversion_list.txt
REHEADER=${BASE}/file_lists/bcftools_reheader.txt

mkdir -p "$LD" "$TMP"

# -------------------------- Select region ------------------------ #

LINE=$(sed -n "${SGE_TASK_ID}p" "$MAPFILE")
CHR=$(echo "$LINE" | awk '{print $1}')
LGC=$(echo "$LINE" | awk '{print $2}')

echo "[$(date)] Starting LD calculation for inversion ${LGC} on chromosome ${CHR}"

AA_LIST=${LD}/${LGC}_common_hom_individuals.txt
TMP_PREFIX=${TMP}/${LGC}
OUT_PREFIX=${LD}/${LGC}

# ------------------------------ Step 1 ---------------------------- #
# Temporary filtered PLINK + VCF
# ------------------------------------------------------------------ #

plink --bfile ${PLINK}/${CHR} \
      --allow-extra-chr \
      --recode \
      --maf 0.05 \
      --out ${TMP_PREFIX}.step1

plink --file ${TMP_PREFIX}.step1 \
      --allow-extra-chr \
      --recode vcf --make-bed \
      --out ${TMP_PREFIX}.step2

# ------------------------------ Step 2 ---------------------------- #
# Reheader VCF (optional)
# ------------------------------------------------------------------ #

bcftools reheader -s ${REHEADER} \
    ${TMP_PREFIX}.step2.vcf > ${TMP_PREFIX}.step3.vcf

# ------------------------------ Step 3 ---------------------------- #
# Extract common-homozygote individuals for this inversion
# ------------------------------------------------------------------ #

vcftools --vcf ${TMP_PREFIX}.step3.vcf \
         --keep ${AA_LIST} \
         --recode --out ${TMP_PREFIX}.AA

# ------------------------------ Step 4 ---------------------------- #
# Convert VCF to PLINK (all individuals and AA-only)
# ------------------------------------------------------------------ #

plink --vcf ${TMP_PREFIX}.step3.vcf \
      --make-bed --allow-extra-chr --const-fid 0 \
      --out ${TMP_PREFIX}.all

plink --vcf ${TMP_PREFIX}.AA.recode.vcf \
      --make-bed --allow-extra-chr --const-fid 0 \
      --out ${TMP_PREFIX}.AA

# ------------------------------ Step 5 ---------------------------- #
# LD calculation (within 100 kb)
# ------------------------------------------------------------------ #

plink --bfile ${TMP_PREFIX}.all \
      --r2 --ld-window-kb 100000 --ld-window-r2 0 --allow-extra-chr \
      --out ${OUT_PREFIX}

plink --bfile ${TMP_PREFIX}.AA \
      --r2 --ld-window-kb 100000 --ld-window-r2 0 --allow-extra-chr \
      --out ${OUT_PREFIX}.AA

# ------------------------------ Step 6 ---------------------------- #
# Cleanup temporary files
# ------------------------------------------------------------------ #

rm -f ${TMP_PREFIX}.step1*
rm -f ${TMP_PREFIX}.step2*
rm -f ${TMP_PREFIX}.step3*
rm -f ${TMP_PREFIX}.all*
rm -f ${TMP_PREFIX}.AA.*

echo "[$(date)] Finished LD for inversion ${LGC} (chromosome ${CHR})"
