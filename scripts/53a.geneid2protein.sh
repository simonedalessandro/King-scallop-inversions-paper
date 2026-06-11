#!/bin/sh
# ==============================================================================
# Script:       53a_geneID2protein.sh
# Description:  Build a mapping between NCBI GeneIDs and RefSeq protein
#               accessions from the P. maximus genomic GFF3 annotation.
#               This mapping is used in 53b_eggNOG-mapper.sh to retrieve
#               protein sequences for enriched genes.
#
# Input:        GCF_902652985.1_xPecMax1.1_genomic.gff.gz (from 50a_GO_prepare.sh)
# Output:       eggNOG-mapper/GeneID2protein.tsv
#               (two columns: GeneID <tab> RefSeq_protein_accession)
#
# Note:         Run once interactively before 53b_eggNOG-mapper.sh.
# ==============================================================================

. /etc/profile.d/modules.sh

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

BASE=/simone/pmax2023/out
GO="${BASE}/40.GO"

GFF="${GO}/GCF_902652985.1_xPecMax1.1_genomic.gff.gz"
OUTDIR="${GO}/eggNOG-mapper"
OUT="${OUTDIR}/GeneID2protein.tsv"

mkdir -p "${OUTDIR}"

echo "=================================================="
echo "Start time: $(date)"
echo "=================================================="

test -s "${GFF}" || { echo "ERROR: missing/empty ${GFF}"; exit 1; }

# ==============================================================================
# Extract GeneID -> protein accession mapping from CDS features
# ==============================================================================
# For each CDS entry:
#   - GeneID is taken from   Dbxref=GeneID:########
#   - Protein accession from Name=XP_########.#
# Only lines where both fields are found are written. Uses 2-argument match()
# + substr (POSIX awk) so it works under gawk, mawk and BSD awk alike.

echo "[$(date)] Parsing GFF3 CDS features..."

zcat "${GFF}" \
| awk -F'\t' 'BEGIN{OFS="\t"}
$3 == "CDS" {
    geneid = ""; prot = "";
    n = split($9, a, ";");
    for (i = 1; i <= n; i++) {
        if (match(a[i], /GeneID:[0-9]+/))
            geneid = substr(a[i], RSTART + 7, RLENGTH - 7);   # digits after "GeneID:"
        if (a[i] ~ /^Name=XP_/)
            prot = substr(a[i], 6);                           # strip "Name="
    }
    if (geneid != "" && prot != "")
        print geneid, prot;
}' \
| sort -u > "${OUT}"

echo "[$(date)] GeneID->protein mappings: $(wc -l < "${OUT}")"

echo "=================================================="
echo "Done."
echo "Output: ${OUT}"
echo "Next step: run 53b_eggNOG-mapper.sh"
echo "End time: $(date)"
echo "=================================================="
