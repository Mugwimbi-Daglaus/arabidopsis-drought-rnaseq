#!/bin/bash
# Step 3: Download TAIR10 genome and build HISAT2 index
# Environment: conda activate plantgenomics

# Download genome
cd /mnt/e/plantgenomics/genome
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

# Download annotation
cd /mnt/e/plantgenomics/annotation
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.57.gtf.gz
gunzip Arabidopsis_thaliana.TAIR10.57.gtf.gz

# Build HISAT2 index
hisat2-build \
    /mnt/e/plantgenomics/genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
    /mnt/e/plantgenomics/genome/tair10_index

echo "Genome index built!"
