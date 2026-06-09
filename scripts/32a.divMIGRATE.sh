#!/bin/bash
# ==============================================================================
# Script:       32a_divMigrate.sh
# Description:  Estimate relative migration rates between populations using
#               diveRsity::divMigrate on a 50 000-SNP subsample of the inversion-free filtered dataset.
#
#                 1.  Subsample 50 000 SNPs from the filtered VCF and write a binary PLINK fileset (PLINK)
#                 2.  PLINK -> GDS (SNPRelate) -> genlight (adegenet) -> Genepop file (dartR) -> divMigrate
#
# Input:        Pmax_no_inv_pruned_filtered.vcf: filtered VCF from 31_prune_filter.sh
#               ID_pop.txt: sample ID / population assignment, tab-separated, header row: ID  Location
#
# Output:       <DATASET>.bed/.bim/.fam
#                divMigrate_results.rds
#                dRelMig.csv / gRelMig.csv / nmRelMig.csv

# ==============================================================================

# ------------------------------------------------------------------------------
# SGE job configuration
# ------------------------------------------------------------------------------
#$ -N divMigrate
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_rss=32G
#$ -pe sharedmem 8
#$ -o output_pmax
#$ -e error_pmax

# ------------------------------------------------------------------------------
# Path / parameter configuration
# ------------------------------------------------------------------------------

BASE_DIR=/simone/pmax2023/out/23.het_DAPC_IBD
OUT_DIR=${BASE_DIR}/divmigrate
VCF_FILE=${BASE_DIR}/Pmax_no_inv_pruned_filtered.vcf
ID_POP=${BASE_DIR}/ID_pop.txt

DATASET=Pmax_no_inv_pruned_filtered_subset
NSUBSAMPLE=50000                                     
BOOTS=100                                             

VCF_BASE=$(basename "${VCF_FILE}" .vcf)                

mkdir -p "${OUT_DIR}"

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load plink/1.9
module load R/4.5

# ==============================================================================
# STEP 1: Subsample 50 000 SNPs from the filtered VCF -> binary PLINK
# ==============================================================================

echo "[$(date)] === Subsampling ${NSUBSAMPLE} SNPs ==="

# Write the full SNP list
plink \
    --vcf              "${VCF_FILE}" \
    --double-id \
    --allow-extra-chr \
    --set-all-var-ids  '@_#' \
    --write-snplist \
    --out              "${OUT_DIR}/${VCF_BASE}"

# Draw a random 50k-SNP subsample of SNP IDs
shuf "${OUT_DIR}/${VCF_BASE}.snplist" | head -n "${NSUBSAMPLE}" \
    > "${OUT_DIR}/${VCF_BASE}_50k.snps"

# Extract the subset and write a binary PLINK fileset
plink \
    --vcf "${VCF_FILE}" \
    --double-id \
    --allow-extra-chr \
    --set-all-var-ids  '@_#' \
    --extract "${OUT_DIR}/${VCF_BASE}_50k.snps" \
    --make-bed \
    --out "${OUT_DIR}/${DATASET}"

# Population map for reference (FID_IID  FID) which assigns populations from ID_pop.txt
awk '{print $1"_"$2, $1}' "${OUT_DIR}/${DATASET}.fam" > "${OUT_DIR}/popmap.txt"

# Clean up intermediate SNP-list files
rm -f "${OUT_DIR}"/*.snps "${OUT_DIR}"/*.snplist

echo "[$(date)] STEP 1 complete: ${OUT_DIR}/${DATASET}.bed/.bim/.fam"

# ==============================================================================
# STEP 2: divMigrate
# ==============================================================================

echo "[$(date)] ===  divMigrate ==="

Rscript - \
  --args \
  "${OUT_DIR}" \
  "${DATASET}" \
  "${ID_POP}" \
  "${BOOTS}" \
<<'EOF'

args    <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- args[1]                 # working directory 
DATASET <- args[2]                 # PLINK fileset basename
ID_POP  <- args[3]                 # ID_pop.txt (header: ID  Location)
BOOTS   <- as.integer(args[4])     # number of  bootstraps

suppressPackageStartupMessages({
  library(adegenet)
  library(dartR)
  library(SNPRelate)
  library(diveRsity)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

bed     <- file.path(OUT_DIR, paste0(DATASET, ".bed"))
bim     <- file.path(OUT_DIR, paste0(DATASET, ".bim"))
famfile <- file.path(OUT_DIR, paste0(DATASET, ".fam"))
gdsfile <- file.path(OUT_DIR, paste0(DATASET, ".gds"))
genfile <- file.path(OUT_DIR, "pmax_50k.gen")

# ==============================================================================
# Convert PLINK to GDS (skip if already exists)
# ==============================================================================

if (!file.exists(gdsfile)) {
  message("Converting PLINK to GDS ...")
  snpgdsBED2GDS(
    bed.fn    = bed,
    bim.fn    = bim,
    fam.fn    = famfile,
    out.gdsfn = gdsfile)
} else {
  message("GDS exists — skipping conversion: ", gdsfile)
}

# ==============================================================================
# Read genotypes and build genlight object
# ==============================================================================
# SNPRelate codes missing genotypes as 3 - convert to NA.

message("Reading genotypes from GDS ...")
gds  <- snpgdsOpen(gdsfile)
geno <- snpgdsGetGeno(gds, with.id = TRUE)
snpgdsClose(gds)

X <- geno$genotype
X[X == 3] <- NA

message("Building genlight object ...")
gl <- new("genlight", X)
indNames(gl) <- geno$sample.id
locNames(gl) <- geno$snp.id

# dartR gl2genepop requires locus names without dots - replace with underscores
locNames(gl) <- gsub("\\.", "_", locNames(gl))

# ==============================================================================
# Assign population labels from ID_pop.txt
# ==============================================================================
# The VCF-based subset does not encode population in FID, so populations are
# taken from ID_pop.txt, matched by sample ID. GDS sample order matches the .fam row order, 
# so assignment is positional.

message("Assigning population labels from ID_pop.txt ...")
fam    <- read.table(famfile,
                     col.names = c("FID","IID","PAT","MAT","SEX","PHENO"))
id_pop <- read.delim(ID_POP, header = TRUE)   # columns: ID, Location

fam$Location <- id_pop$Location[match(fam$IID, id_pop$ID)]

if (any(is.na(fam$Location))) {
  warning("Some samples have no population in ID_pop.txt — check that the .fam ",
          "IID column matches the ID column of ID_pop.txt.")
}

stopifnot(length(indNames(gl)) == nrow(fam))
pop(gl) <- fam$Location

message("Population sizes:")
print(table(pop(gl)))

# ==============================================================================
# Export Genepop file (skip if already exists)
# ==============================================================================

if (!file.exists(genfile)) {
  message("Writing Genepop file (may be slow): ", genfile)
  dartR::gl2genepop(
    gl,
    outfile = basename(genfile),
    outpath = dirname(genfile))
} else {
  message("Genepop exists — skipping export: ", genfile)
}

# ==============================================================================
# Run divMigrate
# ==============================================================================

setwd(OUT_DIR)

message("Running divMigrate ...")
dm <- diveRsity::divMigrate(infile = genfile, stat = "all", outfile = NULL, boots = BOOTS,
filter_threshold = 0, plot_network = TRUE, para = FALSE)

# ==============================================================================
# PART 6: Save outputs
# ==============================================================================

saveRDS(dm, file = "divMigrate_results.rds")

if (!is.null(dm$dRelMig))  write.csv(dm$dRelMig,  "dRelMig.csv")
if (!is.null(dm$gRelMig))  write.csv(dm$gRelMig,  "gRelMig.csv")
if (!is.null(dm$nmRelMig)) write.csv(dm$nmRelMig, "nmRelMig.csv")

message("Done. Outputs written to: ", OUT_DIR)

EOF

echo "[$(date)] === All done ==="
