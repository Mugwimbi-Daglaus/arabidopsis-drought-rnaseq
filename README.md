# Arabidopsis Drought Stress RNA-seq Analysis

## Overview
This project performs a complete RNA-seq differential expression analysis 
of Arabidopsis thaliana drought stress response using both Galaxy 
and Linux (Ubuntu WSL2) pipelines.

## Dataset
- **GEO Accession:** GSE119382
- **Organism:** Arabidopsis thaliana (TAIR10)
- **Conditions:** Drought stress vs Control
- **Replicates:** 3 biological replicates per condition

## Pipeline Steps
1. Quality Control — FastQC
2. Trimming — Trimmomatic
3. Alignment — HISAT2 (TAIR10 genome)
4. Read Counting — featureCounts
5. Differential Expression — DESeq2
6. GO Enrichment Analysis

## Key Results
- 25,726 genes analyzed
- 2,625 upregulated in drought
- 1,845 downregulated in drought
- 1,425 highly significant drought-responsive genes

## Tools Used
- FastQC v0.12.1
- Trimmomatic v0.40
- HISAT2 v2.2.2
- SAMtools v1.19.2
- featureCounts v2.1.1
- DESeq2
- R v4.3.1

## Repository Structure
├── Scripts/     # Pipeline scripts
├── Results/     # Output files and plots
├── README.md    # Project documentation
└── LICENSE      # MIT License

## Author
**Daglaus Mugwimbi**
- GitHub: Mugwimbi-Daglaus
- Email: mugwimbi.daglaus@gmail.com

## Reproducibility
This analysis was run on both Galaxy platform and Linux 
command line confirming reproducibility of results across 
both platforms.
