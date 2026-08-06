# Workflow to the manuscript "Multiple chromosomal inversions shape the genetic structure of a commercial bivalve"

## 1. Summary

This repository contains the full bioinformatic workflow used to identify and characterise
putative chromosomal inversions in the king scallop (*Pecten maximus*) from whole-genome
resequencing data. Starting from raw paired-end short reads, the pipeline:

1. Performs quality trimming, reference alignment, and standard BAM processing/QC.
2. Calls genotype likelihoods and estimates allele frequencies genome-wide using ANGSD.
3. Screens for population structure and relatedness, and produces a quality-filtered sample set.
4. Detects candidate inversion regions genome-wide using windowed local PCA / MDS outlier
   detection (`lostruct`).
5. Characterises each candidate region via local PCA, k-means clustering of karyotypes,
   heterozygosity, and linkage disequilibrium.
6. Produces an inversion-free genome-wide dataset for downstream population-genomic analyses
   (PCA, DAPC, isolation-by-distance, migration rate estimation).
7. Polarises inversion genotypes into ancestral/derived states using outgroup species.
8. Annotates genes within inversion regions and runs GO enrichment and functional annotation
   (eggNOG) on candidate genes.

## 2. Code and package versions

Software/module versions as referenced in the pipeline scripts:

| Tool | Version(s) used |
|---|---|
| Cutadapt | 4.6 |
| TrimGalore | 0.6.6 |
| FastQC | 0.11.9 |
| pigz | 2.3.3 |
| BWA | 0.7.17 / 0.7.18 |
| SAMtools | 1.9 / 1.10 / 1.20 |
| Picard (SortSam, AddOrReplaceReadGroups, MarkDuplicates) | via module system |
| ANGSD | 0.940 |
| PLINK | 1.9 / 2.0.0 |
| BCFtools | 1.16 / 1.20 |
| BEDTools | 2.29.2 / 2.31.1 |
| VCFtools | 0.1.16 |
| ngsRelate | 2 |
| seqtk | 1.4-r122 |
| SRA Toolkit | 2.10.8 |
| R | 4.5 / 4.5.0 |

R packages used across the `.R` scripts: `data.table`, `lostruct`, `ggplot2`, `patchwork`,
`tidyverse`, `dplyr`, `purrr`, `cowplot`, `ggrepel`, `ggthemes`, `wesanderson`, plus
`clusterProfiler` and `AnnotationForge`/`org.Pmaximus.eg.db` (custom OrgDb) for the GO
enrichment stage.

Jobs are configured for a Sun Grid Engine (SGE) high-performance computing cluster (`#$ -N`,
`#$ -t`, `#$ -pe sharedmem`, etc. directives). Users on other schedulers (e.g. Slurm) will need
to translate these headers accordingly.

## 3. Repository structure

```
.
├── README.md
├── file_list/     # Sample sheets and ID lists consumed by the scripts (empty templates
│                   # in this anonymised repository; see Section 5)
└── scripts/        # Numbered pipeline scripts, run in ascending numeric order
```

### `scripts/`

Scripts are numbered to reflect pipeline order. Paired `a`/`b` scripts indicate a shell
wrapper (`a`) plus its associated R analysis/plotting script (`b`), where applicable.

| Stage | Scripts | Description |
|---|---|---|
| Read QC & mapping | `01`–`10` | Adapter/quality trimming (TrimGalore), reference indexing, BWA-MEM alignment, BAM sorting, read-group tagging, duplicate marking, coverage estimation, downsampling to uniform depth, BAM indexing, alignment statistics (flagstat) |
| Genotype likelihoods | `11a`–`11b` | Per-chromosome genotype likelihood and allele frequency estimation (ANGSD), merged into genome-wide Beagle/PLINK datasets |
| Individual QC & relatedness | `12a`–`13b` | Per-individual heterozygosity/missingness (QC), pairwise relatedness (PLINK IBD, ngsRelate) to flag and remove related individuals |
| Population structure | `14a`–`14b` | PCAngsd covariance matrices and PCA visualisation |
| Inversion detection | `20a`–`21b` | VCF/BCF preparation, windowed local PCA + MDS outlier scan (`lostruct`) to flag candidate inversion regions, local PCA/k-means clustering to resolve arrangement karyotypes per region |
| Inversion characterisation | `22`–`25` | Per-region heterozygosity, linkage disequilibrium, and composite MDS/PCA/heterozygosity/LD figures per inversion (including a dedicated script for one double-inversion region) |
| Inversion-free dataset | `26`, `30`–`32a` | Removal of inversion regions from the genome-wide VCF, re-run of genotype calling, LD pruning, QC, and population-genetic summaries (maps, heterozygosity, DAPC, isolation-by-distance, migration rate estimation) on the inversion-free dataset |
| Ancestral state polarisation | `40`–`41` | Outgroup read mapping, consensus building, and polarisation of inversion genotypes to ancestral/derived state using two outgroup scallop species |
| Visualisation | `42` | Alluvial plot of inversion genotype combinations across individuals |
| Functional annotation | `50a`–`53b` | Extraction of genes within inversion regions, GO enrichment (`clusterProfiler`, custom OrgDb via `AnnotationForge`), and functional annotation of enriched genes (eggNOG-mapper) |
