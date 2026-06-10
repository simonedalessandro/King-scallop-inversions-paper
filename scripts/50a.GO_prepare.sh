#!/bin/sh
# ==============================================================================
# Script:       50a_GO_prepare.sh
# Description:  Download NCBI annotation files for P. maximus (GCF_902652985.1)
#               and build three shared inputs for the inversion GO pipeline:
#                 1. all.genes.bed
#                    Genome-wide gene coordinates (BED, 0-based start) with
#                    gene IDs extracted from the GFF3.
#                 2. AnnotationForge/gene2go.tsv
#                    GeneID -> GO term mapping extracted from the GAF file.
#                 3. AnnotationForge/all_GO.genes.txt
#                    Unique numeric GeneIDs with >= 1 GO annotation (background
#                    universe for enrichGO).
#
# Input:        NCBI FTP (downloaded here)
# Output:       all.genes.bed, AnnotationForge/gene2go.tsv,
#               AnnotationForge/all_GO.genes.txt
#
# Note:         Run once interactively before 50b_GO.sh,
#               51_AnnotationForge.sh, and 52_enrichGO.sh.
# ==============================================================================

. /etc/profile.d/modules.sh

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

GO=/simone/pmax2023/out/40.GO          

mkdir -p "${GO}/AnnotationForge"
cd "${GO}"

echo "=================================================="
echo "Working directory: $(pwd)"
echo "Start time:        $(date)"
echo "=================================================="

# ==============================================================================
# STEP 1: Download annotation files from NCBI FTP
# ==============================================================================

echo "[$(date)] === Downloading annotation files ==="

NCBI_FTP="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/652/985/GCF_902652985.1_xPecMax1.1"

wget -nc "${NCBI_FTP}/GCF_902652985.1_xPecMax1.1_genomic.gff.gz"
wget -nc "${NCBI_FTP}/GCF_902652985.1_xPecMax1.1_protein.faa.gz"
wget -nc "${NCBI_FTP}/GCF_902652985.1_xPecMax1.1_gene_ontology.gaf.gz"

echo "[$(date)] Downloads complete"

# ==============================================================================
# STEP 2: Build genome-wide gene BED file
# ==============================================================================
# Extract "gene" features from the GFF3 and output BED (0-based start).
# Gene ID is the ID= attribute value (e.g. ID=gene-LOC117332084).

echo "[$(date)] === Building all.genes.bed ==="

GFF="GCF_902652985.1_xPecMax1.1_genomic.gff.gz"

zcat "${GFF}" \
| awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" {
    match($9, /ID=[^;]+/)
    if (RSTART > 0) {
        id = substr($9, RSTART + 3, RLENGTH - 3)   # strip leading "ID="
        print $1, $4 - 1, $5, id                   # 1-based GFF start -> 0-based BED
    }
}' > all.genes.bed

echo "[$(date)] Genes in all.genes.bed: $(wc -l < all.genes.bed)"

# ==============================================================================
# STEP 3: Build GO background universe from GAF file
# ==============================================================================
# GAF format: column 2 = GeneID, column 5 = GO term, column 7 = evidence code.
# Lines beginning with '!' are comments and are skipped.

echo "[$(date)] === Building GO universe ==="

GAF="GCF_902652985.1_xPecMax1.1_gene_ontology.gaf.gz"

zcat "${GAF}" \
| awk -F'\t' 'BEGIN{OFS="\t"} !/^!/ && $2 != "" {print $2, $5, $7}' \
| sort -u > AnnotationForge/gene2go.tsv

echo "[$(date)] GeneID->GO mappings: $(wc -l < AnnotationForge/gene2go.tsv)"

# Unique numeric GeneIDs with >= 1 GO annotation (enrichGO background universe)
cut -f1 AnnotationForge/gene2go.tsv \
| sort -u > AnnotationForge/all_GO.genes.txt

echo "[$(date)] GO-annotated genes (background universe): $(wc -l < AnnotationForge/all_GO.genes.txt)"

echo "=================================================="
echo "Done."
echo "  all.genes.bed    : ${GO}/all.genes.bed"
echo "  gene2go.tsv      : ${GO}/AnnotationForge/gene2go.tsv"
echo "  all_GO.genes.txt : ${GO}/AnnotationForge/all_GO.genes.txt"
echo ""
echo "Next step: run 50b_GO.sh to extract per-inversion foreground gene sets."
echo "End time: $(date)"
echo "=================================================="
