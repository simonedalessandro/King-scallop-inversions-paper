#!/bin/sh
# ==============================================================================
# Script:       51_AnnotationForge.sh
# Description:  Build and install a custom OrgDb package for Pecten maximus
#               (org.Pmaximus.eg.db) using AnnotationForge. The OrgDb packages
#               the genome-wide GeneID -> GO term mapping into a format usable
#               by clusterProfiler::enrichGO() in 52_enrichGO.sh.
#
# Input:        AnnotationForge/gene2go.tsv from 50a_GO_prepare.sh
#               (columns: GeneID <tab> GO:xxxxxxx <tab> evidence_code)
# Output:       org.Pmaximus.eg.db built in ${OUTDIR} and installed into the
#               user R library (${RLIB})
#
# Note:         AnnotationForge packages an existing GeneID -> GO mapping; it does
#               not generate or infer GO terms. Run once before 52_enrichGO.sh.
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE job configuration
# ------------------------------------------------------------------------------
#$ -N AnnotationForge
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_rss=32G
#$ -pe sharedmem 8
#$ -o output_pmax
#$ -e error_pmax

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load R/4.5

# ------------------------------------------------------------------------------
# Path configuration   
# ------------------------------------------------------------------------------

BASE=/simone/pmax2023/out
GO="${BASE}/40.GO"
ANN="${GO}/AnnotationForge"

GENE2GO="${ANN}/gene2go.tsv"
QC_INV="${ANN}/LGC08.02_GO.genes.txt"   # any per-inversion GO gene list for QC
OUTDIR="${GO}"                          

RLIB=/local/sdaless/R/4.5

mkdir -p "${OUTDIR}" "${RLIB}"

# ------------------------------------------------------------------------------
# Validate input
# ------------------------------------------------------------------------------

test -s "${GENE2GO}" || { echo "ERROR: missing/empty ${GENE2GO}"; exit 1; }

echo "=================================================="
echo "Start time: $(date)"
echo "gene2go:    ${GENE2GO}"
echo "OrgDb out:  ${OUTDIR}"
echo "R library:  ${RLIB}"
echo "=================================================="

# ------------------------------------------------------------------------------
# Run R   (paths passed in via --args)
# ------------------------------------------------------------------------------

Rscript - --args "${GENE2GO}" "${QC_INV}" "${OUTDIR}" "${RLIB}" <<'EOF'

args         <- commandArgs(trailingOnly = TRUE)
GENE2GO_FILE <- args[1]
QC_INV_FILE  <- args[2]
OUTDIR       <- args[3]
RLIB         <- args[4]

# Prepend the user library so install.packages()/library() resolve there
.libPaths(c(RLIB, .libPaths()))

suppressPackageStartupMessages({
  library(AnnotationForge)
  library(AnnotationDbi)
})

# ==============================================================================
# PART 1: Read and clean GeneID -> GO mapping
# ==============================================================================

gene2go <- read.table(GENE2GO_FILE,
                      sep              = "\t",
                      header           = FALSE,
                      stringsAsFactors = FALSE,
                      quote            = "",
                      comment.char     = "")
colnames(gene2go) <- c("GID", "GO", "EVIDENCE")

# Drop duplicates and rows with missing or malformed values
gene2go <- unique(gene2go)
gene2go <- gene2go[!is.na(gene2go$GID)      & gene2go$GID      != "", ]
gene2go <- gene2go[!is.na(gene2go$GO)       & gene2go$GO       != "", ]
gene2go <- gene2go[!is.na(gene2go$EVIDENCE) & gene2go$EVIDENCE != "", ]
gene2go <- gene2go[grepl("^GO:[0-9]{7}$", gene2go$GO), ]   # keep valid GO IDs only

cat("Gene-GO mappings after cleaning:", nrow(gene2go), "\n")
stopifnot(nrow(gene2go) > 0)

# ==============================================================================
# PART 2: Build minimal gene_info table
# ==============================================================================
# SYMBOL is required by makeOrgPackage; reuse GID when gene symbols are unavailable

all_gids <- sort(unique(gene2go$GID))

gene_info <- data.frame(
  GID    = all_gids,
  SYMBOL = all_gids,
  stringsAsFactors = FALSE)

cat("Unique GeneIDs:", nrow(gene_info), "\n")

# ==============================================================================
# PART 3: Create OrgDb source package
# ==============================================================================
# Remove any previous build so makeOrgPackage (which errors on an existing
# package directory) can be re-run cleanly.

pkg_dir <- file.path(OUTDIR, "org.Pmaximus.eg.db")
if (dir.exists(pkg_dir)) unlink(pkg_dir, recursive = TRUE)

makeOrgPackage(
  gene_info  = gene_info,
  go         = gene2go,
  version    = "0.1",
  maintainer = "Simone Dalessandro <your.email@ed.ac.uk>",
  author     = "Simone Dalessandro",
  outputDir  = OUTDIR,
  tax_id     = "6576", 
  genus      = "Pecten",
  species    = "maximus",
  goTable    = "go")

cat("\n[OK] OrgDb source package created:", pkg_dir, "\n")

# ==============================================================================
# PART 4: Install OrgDb into the user library
# ==============================================================================

install.packages(pkg_dir, repos = NULL, type = "source")
library(org.Pmaximus.eg.db)

cat("\n[OK] org.Pmaximus.eg.db installed and loaded.\n")

# ==============================================================================
# PART 5: Sanity checks
# ==============================================================================

cat("\nKeytypes available:\n")
print(keytypes(org.Pmaximus.eg.db))

cat("\nFirst 6 GeneIDs in OrgDb:\n")
print(head(keys(org.Pmaximus.eg.db, keytype = "GID")))

cat("\nExample GeneID -> GO mapping (first 5 GeneIDs):\n")
print(AnnotationDbi::select(org.Pmaximus.eg.db,
                            keys    = head(keys(org.Pmaximus.eg.db, keytype = "GID"), 5),
                            columns = c("GO", "EVIDENCE"),
                            keytype = "GID"))

# ==============================================================================
# PART 6: QC 
# ==============================================================================
# Confirm inversion genes map to the OrgDb

if (file.exists(QC_INV_FILE)) {

  inv_label <- basename(QC_INV_FILE)
  inv       <- scan(QC_INV_FILE, what = "character", quiet = TRUE)
  db_keys   <- keys(org.Pmaximus.eg.db, keytype = "GID")

  inv_in_db <- sum(inv %in% db_keys)
  cat(sprintf("\n[QC] %s genes in OrgDb: %d / %d (%.1f%%)\n",
              inv_label, inv_in_db, length(inv),
              100 * inv_in_db / max(1, length(inv))))

  missing <- setdiff(inv, db_keys)
  if (length(missing) > 0) {
    cat("[QC] Example missing GeneIDs (first 10):",
        paste(head(missing, 10), collapse = ", "), "\n")
  }

} else {
  cat("\n[WARN] QC file not found:", QC_INV_FILE, "- skipping QC.\n")
}

cat("\n[DONE] OrgDb ready. Run 52_enrichGO.sh for GO enrichment.\n")

EOF

echo "=================================================="
echo "End time: $(date)"
echo "=================================================="
