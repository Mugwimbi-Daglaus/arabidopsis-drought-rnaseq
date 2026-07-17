#!/bin/bash
# Step 2: Quality control and trimming with fastp
# fastp replaces both FastQC AND Trimmomatic
# Environment: conda activate plantgenomics

mkdir -p /mnt/e/plantgenomics/results/fastp_reports
mkdir -p /mnt/e/plantgenomics/trimmed

cd /mnt/e/plantgenomics/raw_data

for sample in SRR7779219 SRR7779220 SRR7779221 SRR7779225 SRR7779226 SRR7779227; do
    fastp \
        -i ${sample}_1.fastq.gz \
        -I ${sample}_2.fastq.gz \
        -o /mnt/e/plantgenomics/trimmed/${sample}_1_trimmed.fastq.gz \
        -O /mnt/e/plantgenomics/trimmed/${sample}_2_trimmed.fastq.gz \
        -h /mnt/e/plantgenomics/results/fastp_reports/${sample}_fastp.html \
        -j /mnt/e/plantgenomics/results/fastp_reports/${sample}_fastp.json \
        --thread 4 \
        --detect_adapter_for_pe \
        --qualified_quality_phred 15 \
        --length_required 36
    echo "Done: ${sample}"
done

echo "All trimming complete!"
