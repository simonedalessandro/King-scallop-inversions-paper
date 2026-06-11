#!/bin/sh
# ==============================================================================
# Script:       53b_eggNOG-mapper.sh
# Description:  Functionally annotate enriched inversion genes using eggNOG-mapper.
#               For each inversion (one SGE task):
#                 1. Extract enriched GeneIDs (FDR <= 0.05) from enrichGO output
#                 2. Map GeneIDs to RefSeq protein accessions
#                 3. Extract protein FASTA sequences with seqtk
#                 4. Run eggNOG-mapper on extracted proteins
#                 5. Join protein -> GeneID back onto eggNOG annotations
#
# Input:        Per-inversion enrichGO TSV from 52_enrichGO.sh;
#               GeneID2protein.tsv from 53a_geneID2protein.sh;
#               GCF_902652985.1_xPecMax1.1_protein.faa.gz (from 50a_GO_prepare.sh)
# Output:       Per-inversion eggNOG annotation TSV in eggNOG-mapper/
#
# Note:         Requires the eggnog-mapper conda environment and the eggNOG
#               database directory (${DBDIR}).
#               Manually combine eggNOG annotation TSV and per-inversion enrichGO
#               results TSV from 52_enrichGO.sh
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE array job configuration
# ------------------------------------------------------------------------------
#$ -N eggNOG
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=32G
#$ -t 1-14            # number of inversions
#$ -pe sharedmem 16
#$ -o output_pmax
#$ -e error_pmax

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load anaconda
module load seqtk/1.4-r122

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

BASE=/simone/pmax2023/out
GO="${BASE}/40.GO"
MAPFILE="${BASE}/file_lists/chrom_inversion_list.txt"

ENRICH="${GO}/enrichGO"
EMAP="${GO}/eggNOG-mapper"

PROT="${GO}/GCF_902652985.1_xPecMax1.1_protein.faa.gz"
G2P="${EMAP}/GeneID2protein.tsv"
DBDIR=/simone/software/eggnog-db          # eggNOG database directory

CPU=${NSLOTS:-16}                      

mkdir -p "${EMAP}"

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

echo "=================================================="
echo "Task ID:    ${SGE_TASK_ID}"
echo "Inversion:  ${LGC}"
echo "Start time: $(date)"
echo "=================================================="

# ------------------------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------------------------

test -s "${ENRICH}/${LGC}_enrichGO.tsv" || { echo "ERROR: missing/empty ${ENRICH}/${LGC}_enrichGO.tsv"; exit 1; }
test -s "${G2P}"                        || { echo "ERROR: missing/empty ${G2P}";                        exit 1; }
test -s "${PROT}"                       || { echo "ERROR: missing/empty ${PROT}";                       exit 1; }
test -d "${DBDIR}"                      || { echo "ERROR: missing eggNOG database ${DBDIR}";            exit 1; }

# ==============================================================================
# STEP 1: Extract enriched GeneIDs from enrichGO output (FDR <= 0.05)
# ==============================================================================
# enrichGO TSV (clusterProfiler 4.x layout): column 9 = p.adjust,
# column 11 = geneID (slash-separated GeneIDs).

awk -F'\t' 'NR > 1 && $9 != "" && $9 <= 0.05 {print $11}' "${ENRICH}/${LGC}_enrichGO.tsv" \
    | tr '/' '\n' \
    | awk 'NF' \
    | sort -u > "${EMAP}/${LGC}_geneIDs.enrichGO.txt"

echo "[$(date)] ${LGC}: enriched GeneIDs = $(wc -l < "${EMAP}/${LGC}_geneIDs.enrichGO.txt")"

if [ ! -s "${EMAP}/${LGC}_geneIDs.enrichGO.txt" ]; then
    echo "[$(date)] ${LGC}: no enriched GeneIDs (p.adjust <= 0.05) — skipping eggNOG."
    exit 0
fi

# ==============================================================================
# STEP 2: Map GeneIDs to RefSeq protein accessions
# ==============================================================================
# Exact match on GeneID2protein column 1 (GeneID), emitting column 2 (protein).

awk -F'\t' 'NR==FNR { ids[$1]; next } ($1 in ids) { print $2 }' \
    "${EMAP}/${LGC}_geneIDs.enrichGO.txt" \
    "${G2P}" \
    | sort -u > "${EMAP}/${LGC}_enriched.genes.proteins.txt"

echo "[$(date)] ${LGC}: proteins = $(wc -l < "${EMAP}/${LGC}_enriched.genes.proteins.txt")"

if [ ! -s "${EMAP}/${LGC}_enriched.genes.proteins.txt" ]; then
    echo "[$(date)] ${LGC}: no proteins found for enriched GeneIDs — skipping eggNOG."
    exit 0
fi

# ==============================================================================
# STEP 3: Extract protein FASTA sequences
# ==============================================================================

seqtk subseq \
    "${PROT}" \
    "${EMAP}/${LGC}_enriched.genes.proteins.txt" \
    > "${EMAP}/${LGC}_enriched.genes.proteins.faa"

test -s "${EMAP}/${LGC}_enriched.genes.proteins.faa" \
    || { echo "ERROR: FASTA extraction failed for ${LGC}"; exit 1; }

echo "[$(date)] ${LGC}: FASTA sequences extracted"

# ==============================================================================
# STEP 4: Run eggNOG-mapper
# ==============================================================================
# Activate the eggnog-mapper conda environment

CONDA_BASE=$(conda info --base)
. "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate eggnog-mapper

emapper.py \
    -i "${EMAP}/${LGC}_enriched.genes.proteins.faa" \
    -o "${EMAP}/${LGC}_enriched.genes.eggnog" \
    --itype proteins \
    --tax_scope eukaryota \
    --target_orthologs all \
    --cpu "${CPU}" \
    --go_evidence all \
    --data_dir "${DBDIR}" \
    --override

test -s "${EMAP}/${LGC}_enriched.genes.eggnog.emapper.annotations" \
    || { echo "ERROR: missing eggNOG annotations for ${LGC}"; exit 1; }

echo "[$(date)] ${LGC}: eggNOG-mapper complete"

# ==============================================================================
# STEP 5: Join protein accession -> GeneID back onto eggNOG annotations
# ==============================================================================
# GeneID2protein.tsv has columns: GeneID  protein_accession. Invert it
# (protein -> GeneID) and join on the annotations' first column (protein).

TAB=$(printf '\t')
INV_TMP="${EMAP}/${LGC}.prot2gene.tmp"
ANN_TMP="${EMAP}/${LGC}.annot.tmp"

awk -F'\t' 'BEGIN{OFS="\t"} {print $2, $1}' "${G2P}" \
    | LC_ALL=C sort -t "${TAB}" -k1,1 > "${INV_TMP}"

grep -v '^#' "${EMAP}/${LGC}_enriched.genes.eggnog.emapper.annotations" \
    | LC_ALL=C sort -t "${TAB}" -k1,1 > "${ANN_TMP}"

LC_ALL=C join -t "${TAB}" -1 1 -2 1 "${INV_TMP}" "${ANN_TMP}" \
    > "${EMAP}/${LGC}_enrichGO_eggnog_annotated.tsv"

rm -f "${INV_TMP}" "${ANN_TMP}"

echo "=================================================="
echo "Finished: ${LGC}"
echo "Output:   ${EMAP}/${LGC}_enrichGO_eggnog_annotated.tsv"
echo "End time: $(date)"
echo "=================================================="
