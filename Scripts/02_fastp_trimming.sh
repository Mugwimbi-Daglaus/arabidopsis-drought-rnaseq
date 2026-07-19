#!/bin/bash
# Arabidopsis Drought Stress RNA-seq Pipeline
# Step 2: Quality control and trimming with fastp
# Single-end reads!
# Dataset: GSE119382
# Environment: conda activate plantgenomics

mkdir -p /mnt/e/arab/results/fastp_reports
mkdir -p /mnt/e/arab/trimmed

cd /mnt/e/arab/raw_data

for sample in SRR7779219 SRR7779220 SRR7779221 SRR7779225 SRR7779226 SRR7779227; do
    fastp \
        -i ${sample}.fastq.gz \
        -o /mnt/e/arab/trimmed/${sample}_trimmed.fastq.gz \
        -h /mnt/e/arab/results/fastp_reports/${sample}_fastp.html \
        -j /mnt/e/arab/results/fastp_reports/${sample}_fastp.json \
        --thread 4 \
        --detect_adapter_for_se \
        --qualified_quality_phred 15 \
        --length_required 36
    echo "Done: ${sample}"
done

echo "All trimming complete!"