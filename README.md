# Arabidopsis Drought Stress RNA-seq Analysis

## Overview
This project performs a complete RNA-seq transcriptomics analysis of 
*Arabidopsis thaliana* drought stress response using both Galaxy and 
Linux (Ubuntu WSL2) pipelines, confirming reproducibility across platforms.

## Dataset
- **GEO Accession:** GSE119382
- **BioProject:** PRJNA490054
- **Organism:** *Arabidopsis thaliana* (TAIR10)
- **Conditions:** Drought stress vs Control
- **Replicates:** 3 biological replicates per condition
- **Sequencing:** Illumina HiSeq, Single-end

## Pipeline Steps
1. **Download** - SRA Toolkit (fasterq-dump)
2. **Quality Control & Trimming** - fastp
3. **Genome Indexing** - HISAT2 (TAIR10 genome)
4. **Alignment** - HISAT2
5. **Read Counting** - featureCounts
6. **Differential Expression** - DESeq2
7. **Visualization** - ggplot2, pheatmap
8. **GO Enrichment** - clusterProfiler, org.At.tair.db
9. **Co-expression Network** - WGCNA
10. **Hub Gene Network** - ggraph

## Key Results
- 25,727 genes analyzed
- 2,616 upregulated in drought
- 1,844 downregulated in drought
- 1,421 highly significant drought-responsive genes (padj < 0.05, |LFC| > 1)
- 35 WGCNA co-expression modules identified
- Key modules: Brown (downregulated, r=-0.97) and Turquoise (upregulated, r=+0.89)

## Functional Enrichment
  Top GO biological processes:
- Regulation of transcription, DNA-templated (GO:0006355)
- Response to abscisic acid (GO:0009737)
- Response to water deprivation (GO:0009414)
- Response to cold (GO:0009409)
- Response to salt stress (GO:0009651)

## Tools Used
| Tool | Version | Purpose |
|------|---------|---------|
| fasterq-dump | 3.4.1 | Download raw data |
| fastp | 1.3.6 | QC + Trimming |
| HISAT2 | 2.2.2 | Alignment |
| SAMtools | 1.19.2 | BAM processing |
| featureCounts | 2.1.1 | Read counting |
| DESeq2 | Bioconductor | Differential expression |
| WGCNA | CRAN | Co-expression network |
| ggraph | CRAN | Network visualization |
| clusterProfiler | Bioconductor | GO/KEGG enrichment |
| R | 4.5.3 | Statistical analysis |

## Repository Structure
## How to Run
```bash
# Activate environment
conda activate plantgenomics

# Run pipeline step by step
bash scripts/01_download.sh
bash scripts/02_fastp_trimming.sh
bash scripts/03_genome_index.sh
bash scripts/04_alignment.sh
bash scripts/05_featurecounts.sh
Rscript scripts/06_deseq2_wgcna_analysis.R
```

## Reproducibility
This analysis was performed using both Galaxy and Linux command line
(Ubuntu WSL2) pipelines, confirming reproducibility across platforms.

All six RNA-seq samples showed high alignment efficiency using HISAT2:

- SRR7779219: 98.60%
- SRR7779220: 99.04%
- SRR7779221: 99.03%
- SRR7779225: 98.91%
- SRR7779226: 99.14%
- SRR7779227: 99.12%

Overall alignment rate ranged from 98.60% to 99.14%.

## Author
**Daglaus Mugwimbi**
- GitHub: [Mugwimbi-Daglaus](https://github.com/Mugwimbi-Daglaus)
- Institution: Institute of Botany-CAS, UCAS
- Research: Plant stress genomics and transcriptomics

## License
MIT License
