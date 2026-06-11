#!/bin/sh
# ==============================================================================
# Script:       50b_GO.sh
# Description:  Extract per-inversion gene sets for GO enrichment.
#               For each inversion region (one SGE task per inversion):
#                 1. Build a BED interval for the inversion region
#                 2. Restrict the genome-wide gene BED to the chromosome
#                 3. Intersect to identify genes overlapping the inversion
#                 4. Retain genes with >= 1 GO annotation
#
# Input:        all.genes.bed and AnnotationForge/all_GO.genes.txt from
#               50_GO_prepare.sh; chrom_inversion_list.txt (one inversion per
#               line: index  chr  LGC_name  start  end)
# Output:       Per-inversion foreground GO gene lists in AnnotationForge/
#
# Note:         Run after 50a_GO_prepare.sh.
#               Gene IDs in all.genes.bed have the form "gene-LOC########";
#               these are stripped to numeric IDs to match all_GO.genes.txt.
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE array job configuration
# ------------------------------------------------------------------------------
#$ -N GO
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=16G
#$ -t 1-14           #number of inversions
#$ -pe sharedmem 4
#$ -o output_pmax
#$ -e error_pmax

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load BEDTools/2.31.1

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

BASE=/simone/pmax2023/out              
GO="${BASE}/40.GO"
ANN="${GO}/AnnotationForge"
MAPFILE="${BASE}/file_lists/chrom_inversion_list.txt"

ALL_GENES="${GO}/all.genes.bed"
ALL_GO_GENES="${ANN}/all_GO.genes.txt"

mkdir -p "${GO}" "${ANN}"

# ------------------------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------------------------

test -s "${ALL_GENES}"    || { echo "ERROR: missing/empty ${ALL_GENES}";    exit 1; }
test -s "${ALL_GO_GENES}" || { echo "ERROR: missing/empty ${ALL_GO_GENES}"; exit 1; }
test -s "${MAPFILE}"      || { echo "ERROR: missing/empty ${MAPFILE}";      exit 1; }

# ------------------------------------------------------------------------------
# Parse inversion for this array task
# ------------------------------------------------------------------------------
# chrom_inversion_list.txt columns: index  chr  LGC_name  start  end

read -r _ CHR LGC START END <<EOF
$(sed -n "${SGE_TASK_ID}p" "${MAPFILE}")
EOF

if [ -z "${CHR}" ] || [ -z "${START}" ] || [ -z "${END}" ]; then
    echo "ERROR: no inversion at line ${SGE_TASK_ID} of ${MAPFILE}"
    exit 1
fi

echo "=================================================="
echo "Task ID:    ${SGE_TASK_ID}"
echo "Inversion:  ${LGC}  |  ${CHR}:${START}-${END}"
echo "Start time: $(date)"
echo "=================================================="

# ==============================================================================
# STEP 1: Build inversion BED interval (0-based start)
# ==============================================================================

BED_START=$((START - 1))
printf "%s\t%s\t%s\n" "${CHR}" "${BED_START}" "${END}" > "${GO}/${LGC}.bed"

# ==============================================================================
# STEP 2: Restrict genome-wide gene BED to this chromosome
# ==============================================================================
# Cached per chromosome to avoid redundant awk passes across array tasks.
# Written via a task-unique temp + atomic rename so concurrent tasks sharing a
# chromosome can't read a half-written cache file.

CHR_GENES="${GO}/${CHR}.genes.bed"

if [ ! -s "${CHR_GENES}" ]; then
    tmp="${CHR_GENES}.tmp.${SGE_TASK_ID}.$$"
    awk -v chr="${CHR}" '$1==chr' "${ALL_GENES}" > "${tmp}"
    mv -f "${tmp}" "${CHR_GENES}"
fi

# ==============================================================================
# STEP 3: Intersect chromosome gene BED with inversion interval
# ==============================================================================
# -wa: report the full gene record from the gene BED when an overlap is found.
# Strip the "gene-LOC" prefix to produce numeric IDs matching the GO universe.

bedtools intersect \
    -a "${CHR_GENES}" \
    -b "${GO}/${LGC}.bed" \
    -wa \
    | cut -f4 \
    | sed 's/^gene-LOC//' \
    | sort -u > "${GO}/${LGC}.genes.txt"

echo "[$(date)] ${LGC}: overlapping genes = $(wc -l < "${GO}/${LGC}.genes.txt")"

# ==============================================================================
# STEP 4: GO gene set
# ==============================================================================
# Keep inversion genes that also appear in the GO background universe.
# grep -Fxf: lines of the inversion gene list that exactly match a whole line
# in the universe  (a portable set intersection with no sort/collation
# dependency and no process substitution (unlike comm <(...) <(...)).

grep -Fxf "${ALL_GO_GENES}" "${GO}/${LGC}.genes.txt" \
    > "${ANN}/${LGC}_GO.genes.txt" || true

echo "[$(date)] ${LGC}: foreground GO genes = $(wc -l < "${ANN}/${LGC}_GO.genes.txt")"

echo "=================================================="
echo "Finished: ${LGC}"
echo "End time: $(date)"
echo "=================================================="
