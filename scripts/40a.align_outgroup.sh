#!/bin/sh
# ==============================================================================
# Script:       40a_align_outgroup.sh
# Description:  Map outgroup reads (M. yessoensis + P. magellanicus) to the
#               P. maximus reference with bwa mem, keeping only mapped reads
#               (q >= 20).
# Called by:    40_polarise.sh (STEP 2)
# Output:       ${BAM_DIR}/<SRA>.mapped.bam
# Note:         Index the reference first:  bwa index "${REFERENCE}"
# ==============================================================================

# ------------------------------------------------------------------------------
# SGE job configuration
# ------------------------------------------------------------------------------
#$ -N align_outgroup
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=2G
#$ -pe sharedmem 6
#$ -t 1-15
#$ -o o_files
#$ -e e_files

# ------------------------------------------------------------------------------
# Modules
# ------------------------------------------------------------------------------
. /etc/profile.d/modules.sh

module load bwa/0.7.17
module load samtools/1.10

# ------------------------------------------------------------------------------
# Paths / parameters
# ------------------------------------------------------------------------------
RAW_DIR=/simone/pmax2023/raw/outgroup
BAM_DIR=/simone/pmax2023/out/30.ancestral/outgroup/bam
REFERENCE=/simone/pmax2023/reference/GCA_902652985.1_xPecMax1.1_genomic.fna
THREADS=6

# 15 outgroup accessions (10 M. yessoensis + 5 P. magellanicus)
SRAS="SRR18361775 SRR18361774 SRR18361772 SRR18361771 SRR18361770 \
      SRR18361769 SRR18361773 SRR18361777 SRR18361766 SRR18361826 \
      SRR18361812 SRR18361820 SRR18361815 SRR18361919 SRR18361818"

# Pick this task's accession (awk treats runs of whitespace as one separator)
SRA=$(echo "${SRAS}" | awk -v i="${SGE_TASK_ID}" '{print $i}')

echo "Processing ${SRA} on ${HOSTNAME}"

bwa mem -t "${THREADS}" "${REFERENCE}" \
        "${RAW_DIR}/${SRA}_1.fastq" "${RAW_DIR}/${SRA}_2.fastq" \
  | samtools view -q 20 -bF 4 - > "${BAM_DIR}/${SRA}.mapped.bam"

echo "Done: ${BAM_DIR}/${SRA}.mapped.bam"
