#!/bin/bash
# Step 4: Align samples to TAIR10 genome
# Processes one sample at a time to save disk space
# Environment: conda activate plantgenomics

for sample in SRR7779219 SRR7779220 SRR7779221 SRR7779225 SRR7779226 SRR7779227; do
    # Align with HISAT2
    hisat2 -x /mnt/e/plantgenomics/genome/tair10_index \
        -1 /mnt/e/plantgenomics/trimmed/${sample}_1_trimmed.fastq.gz \
        -2 /mnt/e/plantgenomics/trimmed/${sample}_2_trimmed.fastq.gz \
        -S /mnt/e/plantgenomics/aligned/${sample}.sam \
        --summary-file /mnt/e/plantgenomics/aligned/${sample}_stats.txt \
        -p 4

    # Convert SAM to BAM immediately
    /usr/bin/samtools sort /mnt/e/plantgenomics/aligned/${sample}.sam \
        -o /mnt/e/plantgenomics/aligned/${sample}.bam
    /usr/bin/samtools index /mnt/e/plantgenomics/aligned/${sample}.bam

    # Delete SAM to save space
    rm /mnt/e/plantgenomics/aligned/${sample}.sam

    # Delete trimmed files to save space
    rm /mnt/e/plantgenomics/trimmed/${sample}*.fastq.gz

    echo "Done: ${sample}"
done

echo "All samples aligned!"
