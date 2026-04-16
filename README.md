# Kcnq1ot1/KCNQ1OT1 in Metabolic Liver Disease: Integrative Transcriptome Re-Analysis

Analysis scripts for the manuscript:

> **Model-dependent regulation of lncRNA Kcnq1ot1/KCNQ1OT1 in metabolic liver disease models and human NAFLD: inverse correlation with oxidative phosphorylation gene expression**
>
> Hidenori Tani
> Yokohama University of Pharmacy

## Overview

This repository contains all R scripts used for the integrative re-analysis of public RNA-seq datasets (GSE246221 and GSE207855) to investigate the role of lncRNA Kcnq1ot1/KCNQ1OT1 in metabolic liver disease progression.

### Key Findings

- Kcnq1ot1/KCNQ1OT1 regulation is **model-dependent**: unchanged in chemical hepatotoxicity (CCl4), upregulated in metabolic disease (STZ+HFD, human NAFLD)
- Human KCNQ1OT1 is significantly upregulated in NAFLD (log2FC = +0.88, padj = 6.4 x 10^-8)
- Strong inverse correlation with **mitochondrial oxidative phosphorylation** genes (KEGG adjusted p = 3.4 x 10^-35)
- Cross-species conservation of the Kcnq1ot1-mitochondrial axis

## Data Sources

Raw count matrices were downloaded from NCBI GEO:

| Dataset | Description | Samples |
|---------|-------------|---------|
| [GSE246221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246221) | Mouse STZ+HFD model + Human NAFLD | 52 mouse + 32 human |
| [GSE207855](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207855) | Mouse CCl4 fibrosis model | 12 |

To reproduce the analyses, download the following files from GEO and place them in a `data/` directory:
- `GSE246221_rawcounts_allsamples_mouse.txt.gz`
- `GSE246221_rawcounts_allsamples_human.txt.gz`
- `GSE207855_count.tsv.gz`

## Scripts

Run scripts in numerical order:

| Script | Description |
|--------|-------------|
| `scripts/01_data_overview.R` | Data loading, metadata parsing, TMM normalization, initial expression survey |
| `scripts/02_Kcnq1ot1_comprehensive.R` | DESeq2 differential expression (5 comparisons), genome-wide co-expression, GO enrichment |
| `scripts/03_KEGG_and_network.R` | KEGG pathway enrichment, mitochondrial complex gene analysis |
| `scripts/04_human_NAFLD_severity.R` | Human NAFLD DESeq2 analysis, severity marker correlations |
| `scripts/05_publication_figures.R` | Individual publication-quality figures |
| `scripts/06_main_figures.R` | Multi-panel main figures for the manuscript (Fig. 1-6) |

## Requirements

R (>= 4.0) with the following packages:

```r
# Bioconductor
BiocManager::install(c("DESeq2", "edgeR", "clusterProfiler",
                       "org.Mm.eg.db", "EnhancedVolcano"))

# CRAN
install.packages(c("ggplot2", "pheatmap", "patchwork", "cowplot"))
```

## Output

- `results/` — CSV files with DEG tables, correlation data, pathway enrichment results
- `figures/` — Individual figure files (PDF + PNG)
- `figures/main/` — Multi-panel main figures for the manuscript

## License

MIT License

## Citation

If you use these scripts, please cite:

> Tani, H. Model-dependent regulation of lncRNA Kcnq1ot1/KCNQ1OT1 in metabolic liver disease models and human NAFLD: inverse correlation with oxidative phosphorylation gene expression. *Hepatology Research* (submitted).

## Contact

Hidenori Tani — Yokohama University of Pharmacy — hidenori.tani@yok.hamayaku.ac.jp
