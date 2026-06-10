#!/bin/sh
# ===================================================================================================================
# Script:       40_polarise.sh
# Description:  Polarise alleles in the P. maximus inversion dataset using a
#               consensus sequence built from two outgroup species mapped to
#               the P. maximus reference genome:
#                 - Mizuhopecten yessoensis   (10 samples; subfamily Pedinae)
#                 - Placopecten magellanicus  ( 5 samples; subfamily Palliolinae)
#
#               Pipeline:
#                 1. (Interactive) Download outgroup reads from SRA
#                 2. (SGE job)     Map reads to reference   -> 40a_align_outgroup.sh
#                 3. (SGE job)     Coordinate-sort BAMs     -> 40b_SortSam_outgroup.sh
#                 4. (Interactive) Per-sample coverage
#                    (SGE job)     Merge BAMs               -> 40c_MergeSam_outgroup.sh
#                 5. (SGE job)     Consensus FASTA (ANGSD)  -> 40d_doFasta_outgroup.sh
#                    (Interactive) Extract ancestral allele per SNP
#                 6. (Interactive) Rotate REF to ancestral allele (PLINK2)
#
# Input:        Outgroup SRA SRAessions (listed in STEP 1)
#               P. maximus reference genome (FASTA)
#               Inversion SNPs dataset (VCF + PLINK .map) of 160 samples
#
# Output:       ${POL_DIR}/dataset_inversion_polarised.vcf
#               ${POL_DIR}/<INV>/<INV>.traw - per-inversion ancestral-allele
#                   dosage matrix for downstream R load analysis
#               ${OUT_DIR}/ancestral_alleles.txt - SNP_ID -> ancestral allele
#
# Usage:        Run interactively (qlogin). Steps 2, 3, the merge in step 4, and
#               step 5 are submitted as separate SGE jobs (qsub) and must finish
#               before the next step is run. After polarisation: 0/0 = ancestral, 0/1 = heterozygote, 1/1 = derived.
# ===================================================================================================================

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load sratoolkit/2.10.8
module load samtools/1.10
module load bedtools/2.29.2
module load bcftools/1.16
module load vcftools/0.1.16
module load plink/2.0.0
module load R/4.5.0

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

BASE=/simone/pmax2023
RAW_DIR=${BASE}/raw/outgroup                       
INV_DIR=${BASE}/out/20.lostruct/inversions          # inversion SNP dataset
OUT_DIR=${BASE}/out/30.ancestral/outgroup           
BAM_DIR=${OUT_DIR}/bam                              
COV_DIR=${BAM_DIR}/cov                             
POL_DIR=${OUT_DIR}/polarised                        # final polarised VCF

REFERENCE=${BASE}/reference/GCA_902652985.1_xPecMax1.1_genomic.fna
VCF_FILE=${INV_DIR}/dataset_inversion.vcf
MAP_FILE=${INV_DIR}/dataset_inversion.map


mkdir -p "${RAW_DIR}" "${BAM_DIR}" "${COV_DIR}" "${POL_DIR}"

echo "=================================================="
echo "Start time: $(date)"
echo "=================================================="

# ==============================================================================
# STEP 1: Download outgroup reads from SRA
# ==============================================================================
# fasterq-dump --split-files writes paired reads as ${SRA}_1.fastq / ${SRA}_2.fastq
# -e 6 uses 6 threads.

echo "[$(date)] === STEP 1: Downloading SRA reads ==="

cd "${RAW_DIR}"

# Mizuhopecten yessoensis (10 samples)
for SRA in SRR18361775 SRR18361774 SRR18361772 SRR18361771 SRR18361770 \
           SRR18361769 SRR18361773 SRR18361777 SRR18361766 SRR18361826; do
    echo "[$(date)] Downloading ${SRA}"
    fasterq-dump --split-files -e 6 "${SRA}"
done

# Placopecten magellanicus (5 samples)
for SRA in SRR18361812 SRR18361820 SRR18361815 SRR18361919 SRR18361818; do
    echo "[$(date)] Downloading ${SRA}"
    fasterq-dump --split-files -e 6 "${SRA}"
done

echo "[$(date)] SRA downloads complete"

# ==============================================================================
# STEP 2: Align outgroup reads to the P. maximus reference  (SGE job)
# ==============================================================================
# Array job over the 15 outgroup accessions. Make sure the reference is indexed
# first ( bwa index "${REFERENCE}" ). Submit and wait for completion:
#
#   qsub 40a_align_outgroup.sh
#
# Output: ${BAM_DIR}/<SRA>.mapped.bam  (q20, mapped reads only)

echo "[$(date)] === STEP 2: submit 40a_align_outgroup.sh and wait ==="

# ==============================================================================
# STEP 3: Coordinate-sort BAMs  (SGE job)
# ==============================================================================
#   qsub 40b_SortSam_outgroup.sh
#
# Output: ${BAM_DIR}/<SRA>.sorted.bam (+ .bai)

echo "[$(date)] === STEP 3: submit 40b_SortSam_outgroup.sh and wait ==="

# ==============================================================================
# STEP 4: Per-sample coverage, then merge BAMs
# ==============================================================================

echo "[$(date)] === STEP 4: Coverage and BAM merge ==="

cd "${BAM_DIR}"

# Per-sample mean depth (samtools coverage column 7 = meandepth; skip header row)
MYESS="SRR18361775 SRR18361774 SRR18361772 SRR18361771 SRR18361770
       SRR18361769 SRR18361773 SRR18361777 SRR18361766 SRR18361826"
PMAG="SRR18361812 SRR18361820 SRR18361815 SRR18361919 SRR18361818"

i=1
for SRA in ${MYESS}; do
    samtools coverage "${SRA}.sorted.bam" \
        | awk 'NR > 1 {print $7}' > "${COV_DIR}/myessoensis${i}_cov.txt"
    i=$((i + 1))
done

i=1
for SRA in ${PMAG}; do
    samtools coverage "${SRA}.sorted.bam" \
        | awk 'NR > 1 {print $7}' > "${COV_DIR}/pmagellanicus${i}_cov.txt"
    i=$((i + 1))
done

echo "[$(date)] Coverage files written to ${COV_DIR}"

# Merge all sorted outgroup BAMs into a single file for consensus calling (SGE job):
#   qsub 40c_MergeSam_outgroup.sh
# Output: ${BAM_DIR}/outgroup_merged.bam

echo "[$(date)] === STEP 4: submit 40c_MergeSam_outgroup.sh and wait ==="

# ==============================================================================
# STEP 5: Outgroup consensus, then extract ancestral alleles at SNP positions
# ==============================================================================
# Consensus called with ANGSD -doFasta 2 -doCounts 1 (most-common base) (SGE job):
#   qsub 40d_doFasta_outgroup.sh
# Output: ${OUT_DIR}/outgroup_consensus.fa

echo "[$(date)] === STEP 5: submit 40d_doFasta_outgroup.sh and wait ==="

cd "${OUT_DIR}"

# ── 5a: Build a BED of SNP positions from the PLINK .map ──────────────────────
# BED uses 0-based, half-open coordinates; PLINK .map positions are 1-based.

Rscript - --args "${MAP_FILE}" "${OUT_DIR}/dataset_inversion.bed" <<'REOF'
args     <- commandArgs(trailingOnly = TRUE)
MAP_FILE <- args[1]
BED_FILE <- args[2]

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# PLINK .map columns: V1 = chr, V2 = SNP ID, V3 = cM, V4 = bp (1-based)
map <- fread(MAP_FILE)

bed <- map %>%
  transmute(chr   = V1,
            start = V4 - 1,   # 1-based -> 0-based
            end   = V4,
            name  = "SNP")

fwrite(bed, BED_FILE, quote = FALSE, row.names = FALSE,
       col.names = FALSE, sep = "\t")

cat("BED file written:", BED_FILE, "\n")
cat("SNP count:", nrow(bed), "\n")
REOF

# ── 5b: Pull the consensus base at each SNP position from the outgroup FASTA ──
bedtools getfasta \
    -fi  "${OUT_DIR}/outgroup_consensus.fa" \
    -bed "${OUT_DIR}/dataset_inversion.bed" \
    -fo  "${OUT_DIR}/ancestral_alleles.out"

echo "[$(date)] ancestral_alleles.out written"

# ── 5c: Reformat to a two-column PLINK table: SNP_ID  ancestral_base ──────────
# SNPs where the outgroup base is 'N' (no coverage / ambiguous) are dropped.

Rscript - --args "${OUT_DIR}/ancestral_alleles.out" "${MAP_FILE}" \
                 "${OUT_DIR}/ancestral_alleles.txt" <<'REOF'
args     <- commandArgs(trailingOnly = TRUE)
ANC_FILE <- args[1]
MAP_FILE <- args[2]
OUT_FILE <- args[3]

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

anc <- fread(ANC_FILE, header = FALSE)
map <- fread(MAP_FILE)                       # V2 = SNP ID

# getfasta interleaves ">chr:start-end" headers with single-base allele lines,
# in BED order. The BED was built from the .map (same order), so the i-th allele
# corresponds to the i-th .map row.
alleles <- toupper(anc$V1[!grepl("^>", anc$V1)])
stopifnot(length(alleles) == nrow(map))

out <- data.frame(SNP_ID = map$V2, ANC = alleles, stringsAsFactors = FALSE)

n_N <- sum(out$ANC == "N")
cat("SNPs with N (no outgroup coverage, excluded):", n_N, "\n")

out <- out[out$ANC != "N", ]

fwrite(out, OUT_FILE, quote = FALSE, col.names = FALSE,
       row.names = FALSE, sep = " ")

cat("ancestral_alleles.txt written:", OUT_FILE, "\n")
cat("Informative SNPs written:", nrow(out), "\n")
REOF

# ==============================================================================
# STEP 6: Polarise alleles with PLINK2
# ==============================================================================
# --ref-allele rotates allele codes so the ancestral allele is REF (0).
# Assumes the .map SNP IDs match the variant IDs in the VCF (same PLINK dataset);
# --set-missing-var-ids only fills IDs for variants whose ID is '.'.

echo "[$(date)] === STEP 6: Polarising alleles ==="

cd "${OUT_DIR}"

# Ensure every variant has an ID (chr:pos for any that are missing)
plink2 \
    --vcf "${VCF_FILE}" \
    --set-missing-var-ids @:# \
    --recode vcf \
    --allow-extra-chr \
    --out temp

echo "[$(date)] Variant IDs ensured: temp.vcf"

# Keep SNPs with a known ancestral allele and set it as REF
awk '{print $1}' ancestral_alleles.txt > ancestral_positions.txt

plink2 \
    --vcf temp.vcf \
    --extract ancestral_positions.txt \
    --ref-allele force ancestral_alleles.txt 2 1 \
    --recode vcf \
    --allow-extra-chr \
    --out temp2

echo "[$(date)] Alleles rotated: temp2.vcf"

# Collect variants where the ancestral allele matched neither REF nor ALT
# (uninformative) from the PLINK log, so they can be removed
grep 'Warning' temp2.log \
    | awk '{print $7}' \
    | sed "s/^'//;s/'\.$//" \
    > mismatches.txt

echo "[$(date)] Mismatch variants to exclude: $(wc -l < mismatches.txt)"

# Write the final polarised VCF
plink2 \
    --vcf temp2.vcf \
    --exclude mismatches.txt \
    --allow-extra-chr \
    --export vcf-4.2 \
    --threads 1 \
    --out "${POL_DIR}/dataset_inversion_polarised"

echo "[$(date)] Polarised VCF written: ${POL_DIR}/dataset_inversion_polarised.vcf"
echo "[$(date)] Final SNP count: $(grep -vc '^#' "${POL_DIR}/dataset_inversion_polarised.vcf")"

# Tidy intermediate VCFs
rm -f temp.vcf temp2.vcf

# ------------------------------------------------------------------------------
# Allele coding after polarisation  (REF forced to ancestral in STEP 6)
# ------------------------------------------------------------------------------
# VCF (dataset_inversion_polarised.vcf):
#   0/0 = ancestral homozygous (REF/REF)
#   0/1 = heterozygous         (REF/ALT)
#   1/1 = derived homozygous   (ALT/ALT)
#
# .traw (PLINK2 --export A-transpose): with REF = ancestral, COUNTED = REF, so
# each value = number of copies of the ANCESTRAL allele per individual:
#   traw 2 -> 0/0 -> ancestral homozygous (A/A)
#   traw 1 -> 0/1 -> heterozygous         (A/D)
#   traw 0 -> 1/1 -> derived homozygous   (D/D)
#
# ==============================================================================
# STEP 7: Per-inversion ancestral-allele .traw export
# ==============================================================================
# For each inversion, subset the polarised VCF to its region and export an
# ancestral-allele dosage matrix (.traw) for visualisation/load analysis in R.
# Add one "<name> <chr> <from_bp> <to_bp>" line per inversion below.

echo "[$(date)] === STEP 8: Per-inversion .traw export ==="

INVERSIONS="LGC15 LR736852.1 53 15458540"
# e.g. add more, one per line:
# INVERSIONS="${INVERSIONS}
# LGC06 LR736843.1 12345 6789012"

# Loop via a temp file (not a pipe) so 'exit 1' below halts the script rather
# than just a pipe subshell.
inv_file=$(mktemp)
printf '%s\n' "${INVERSIONS}" > "${inv_file}"

while read -r INV CHR FROM TO; do
    [ -z "${INV}" ] && continue
    echo "[$(date)] ${INV}: ${CHR}:${FROM}-${TO}"

    INV_OUT=${POL_DIR}/${INV}
    mkdir -p "${INV_OUT}"

    # Subset to the inversion region
    vcftools \
        --vcf dataset_inversion_polarised.vcf \
        --chr "${CHR}" --from-bp "${FROM}" --to-bp "${TO}" \
        --recode --out "${INV_OUT}/${INV}"

    # Export ancestral-allele dosage matrix
    plink2 \
        --vcf "${INV_OUT}/${INV}.recode.vcf" \
        --export A-transpose \
        --allow-extra-chr \
        --threads 1 \
        --out "${INV_OUT}/${INV}"

    # Sanity check
    n_mismatch=$(awk '
        FNR == NR { if ($0 !~ /^#/) ref[++i] = $4; next }   # region VCF REF, in order
        FNR > 1   { j++; if ($5 != ref[j]) c++ }            # .traw COUNTED
        END       { print c + 0 }
    ' "${INV_OUT}/${INV}.recode.vcf" "${INV_OUT}/${INV}.traw")

    if [ "${n_mismatch}" -ne 0 ]; then
        echo "[$(date)] ERROR: ${INV}: ${n_mismatch} SNPs where COUNTED != REF."
        echo "[$(date)] .traw counts the wrong allele — do NOT use for load counts."
        exit 1
    fi

    echo "[$(date)] ${INV}: COUNTED == REF — .traw counts the ANCESTRAL allele"
    echo "[$(date)] ${INV}: written ${INV_OUT}/${INV}.traw"
done < "${inv_file}"

rm -f "${inv_file}"

# Per-karyotype individual lists (AA / BB homozygotes), used by the R viz, are
# derived from the inversion-clustering output
# (/simone/pmax2023/out/20.lostruct/inversions/local_pca/<INV>/)
# e.g. for LGC15:
#   awk 'NR>1{$1=""; gsub("\"",""); print}' AA.list.txt | sed 's/^[ \t]*//' \
#       > LGC15_AA_individuals.txt
#   awk 'NR>1{$1=""; gsub("\"",""); print}' BB.list.txt | sed 's/^[ \t]*//' \
#       > LGC15_BB_individuals.txt

echo "=================================================="
echo "Done."
echo "Output: ${POL_DIR}/dataset_inversion_polarised.vcf"
echo "        ${POL_DIR}/<INV>/<INV>.traw  (per inversion)"
echo "End time: $(date)"
echo "=================================================="
