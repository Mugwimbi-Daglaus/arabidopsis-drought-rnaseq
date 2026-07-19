#!/bin/bash
# Arabidopsis Drought Stress RNA-seq Pipeline
# Step 4: Align samples to TAIR10 genome
# Single-end reads — uses -U flag!
# Processes one sample at a time to save disk space
# Environment: conda activate plantgenomics

for sample in SRR7779219 SRR7779220 SRR7779221 SRR7779225 SRR7779226 SRR7779227; do
    # Align with HISAT2 (single-end)
    hisat2 \
        -x /mnt/e/arab/genome/tair10_index \
        -U /mnt/e/arab/trimmed/${sample}_trimmed.fastq.gz \
        -S /mnt/e/arab/aligned/${sample}.sam \
        --summary-file /mnt/e/arab/aligned/${sample}_stats.txt \
        -p 4

    # Convert SAM to BAM immediately
    /usr/bin/samtools sort \
        /mnt/e/arab/aligned/${sample}.sam \
        -o /mnt/e/arab/aligned/${sample}.bam
    /usr/bin/samtools index \
        /mnt/e/arab/aligned/${sample}.bam

    # Delete SAM to save space
    rm /mnt/e/arab/aligned/${sample}.sam

    # Delete trimmed files to save space
    rm /mnt/e/arab/trimmed/${sample}_trimmed.fastq.gz

    echo "Done: ${sample}"
done

echo "All samples aligned!"