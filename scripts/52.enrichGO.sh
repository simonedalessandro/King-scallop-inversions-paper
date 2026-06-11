#!/bin/sh
# ==============================================================================
# Script:       52_enrichGO.sh
# Description:  Run GO enrichment analysis (clusterProfiler::enrichGO) for each
#               inversion region. One SGE task per inversion; uses the custom
#               P. maximus OrgDb built in 51_AnnotationForge.sh.
#
# Input:        Per-inversion  gene lists (<LGC>_GO.genes.txt) from
#               50b_GO.sh; background universe (all_GO.genes.txt) from
#               50a_GO_prepare.sh; org.Pmaximus.eg.db installed in
#               51_AnnotationForge.sh
# Output:       Per-inversion enrichGO results TSV in enrichGO/
#
# Note:         Gene IDs must match the OrgDb keyType "GID" (numeric NCBI GeneIDs).
#               p-values adjusted by Benjamini-Hochberg (BH).
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE array job configuration
# ------------------------------------------------------------------------------
#$ -N enrichGO
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_rss=16G
#$ -t 1-14        #number of inversions
#$ -pe sharedmem 4
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
MAPFILE="${BASE}/file_lists/chrom_inversion_list.txt"

ALL_GO_GENES="${ANN}/all_GO.genes.txt"
OUTDIR="${GO}/enrichGO"

# User R library holding org.Pmaximus.eg.db (edit to match your account)
RLIB=/local/sdaless/R/4.5

mkdir -p "${OUTDIR}"

# ------------------------------------------------------------------------------
# Parse inversion for this array task
# ------------------------------------------------------------------------------
# chrom_inversion_list.txt columns: index  chr  LGC_name  start  end

test -s "${MAPFILE}" || { echo "ERROR: missing/empty ${MAPFILE}"; exit 1; }

read -r _ _ LGC _ _ <<EOF
$(sed -n "${SGE_TASK_ID}p" "${MAPFILE}")
EOF

if [ -z "${LGC}" ]; then
    echo "ERROR: no inversion at line ${SGE_TASK_ID} of ${MAPFILE}"
    exit 1
fi

FG="${ANN}/${LGC}_GO.genes.txt"
OUT="${OUTDIR}/${LGC}_enrichGO.tsv"

echo "=================================================="
echo "Task ID:    ${SGE_TASK_ID}"
echo "Inversion:  ${LGC}"
echo "Start time: $(date)"
echo "=================================================="

# ------------------------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------------------------

test -s "${ALL_GO_GENES}" || { echo "ERROR: missing/empty ${ALL_GO_GENES}"; exit 1; }
test -s "${FG}"           || { echo "ERROR: missing/empty ${FG}";           exit 1; }

# ------------------------------------------------------------------------------
# Run GO enrichment in R
# ------------------------------------------------------------------------------

Rscript - --args "${FG}" "${ALL_GO_GENES}" "${OUT}" "${RLIB}" <<'EOF'

args     <- commandArgs(trailingOnly = TRUE)
FG_FILE  <- args[1]
BG_FILE  <- args[2]
OUT_FILE <- args[3]
RLIB     <- args[4]

# Prepend the user library so org.Pmaximus.eg.db resolves there
.libPaths(c(RLIB, .libPaths()))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Pmaximus.eg.db)
})

# ── Read foreground and background gene lists ──────────────────────────────────
fg <- unique(scan(FG_FILE, what = "character", quiet = TRUE))
bg <- unique(scan(BG_FILE, what = "character", quiet = TRUE))

cat("[R] Foreground genes:", length(fg), "\n")
cat("[R] Background genes:", length(bg), "\n")

# ── Run enrichment (Biological Process) ───────────────────────────────────────
GO_result <- enrichGO(
  gene          = fg,
  OrgDb         = org.Pmaximus.eg.db,
  keyType       = "GID",
  ont           = "BP",   # biological process
  universe      = bg,
  pAdjustMethod = "BH",   # adj. with Benjamini-Hochberg
  pvalueCutoff  = 0.1,
  qvalueCutoff  = 0.1,
  readable      = FALSE)

# enrichGO returns NULL when nothing maps / nothing is enriched
if (is.null(GO_result)) {
  cat("[R] enrichGO returned NULL (no mappable genes or no enrichment)\n")
  res <- data.frame()
} else {
  res <- as.data.frame(GO_result)
}

# ── Save results ───────────────────────────────────────────────────────────────
write.table(res,
            file      = OUT_FILE,
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

cat("[R] Significant GO terms:", nrow(res), "\n")
cat("[R] Wrote:", OUT_FILE, "\n")

EOF

echo "=================================================="
echo "Finished: ${LGC}"
echo "End time: $(date)"
echo "=================================================="
