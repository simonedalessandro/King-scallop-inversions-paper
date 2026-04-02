#!/bin/sh
# ==============================================================================
# Script:       02.IndexRef.sh
# Description:  Index the P. maximus reference genome

# Output: .fai index (SAMtools) and BWA index files alongside the FASTA
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE job configuration
# ------------------------------------------------------------------------------
#$ -N IndexRef
#$ -cwd
#$ -l h_rt=01:00:00
#$ -l h_rss=16G
#$ -pe sharedmem 4
#$ -o o_files
#$ -e e_files

# ------------------------------------------------------------------------------
# Load required modules
# ------------------------------------------------------------------------------

. /etc/profile.d/modules.sh

module load bwa/0.7.18    
module load samtools/1.20    

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

REFERENCE="/simone/reference/Pmax/GCA_902652985.1_xPecMax1.1_genomic.fna"

# ------------------------------------------------------------------------------
# Index reference genome
# ------------------------------------------------------------------------------

echo "=================================================="
echo "Reference:  ${REFERENCE}"
echo "Start time: $(date)"
echo "=================================================="

# Generate .fai index
samtools faidx "${REFERENCE}"

# Generate BWA index files (.amb .ann .bwt .pac .sa) required for alignment
bwa index "${REFERENCE}"

echo "=================================================="
echo "Indexing complete"
echo "End time: $(date)"
echo "=================================================="
