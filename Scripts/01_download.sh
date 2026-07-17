#!/bin/bash
# Arabidopsis Drought Stress RNA-seq Pipeline
# Step 1: Download samples using fasterq-dump
# Dataset: GSE119382
# Conditions: Drought vs Control
# Environment: conda activate plantgenomics

cd /mnt/e/plantgenomics/raw_data

for sample in SRR7779219 SRR7779220 SRR7779221 SRR7779225 SRR7779226 SRR7779227; do
    prefetch ${sample}
    fasterq-dump --split-files ${sample}
    pigz ${sample}_1.fastq ${sample}_2.fastq
    echo "Downloaded: ${sample}"
done

echo "All downloads complete!"
