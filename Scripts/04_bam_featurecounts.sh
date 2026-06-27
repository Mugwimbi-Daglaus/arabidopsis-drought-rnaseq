#!/bin/bash
# Step 4: Convert SAM to BAM and run featureCounts

# Convert SAM to BAM
for sample in SRR7779219 SRR7779220 SRR7779221 SRR7779225 SRR7779226 SRR7779227; do
    /usr/bin/samtools sort /mnt/e/plantgenomics/aligned/${sample}.sam \
        -o /mnt/e/plantgenomics/aligned/${sample}.bam
    /usr/bin/samtools index /mnt/e/plantgenomics/aligned/${sample}.bam
    rm /mnt/e/plantgenomics/aligned/${sample}.sam
    echo "Done converting: ${sample}"
done

# Run featureCounts
featureCounts \
    -a /mnt/e/plantgenomics/annotation/Arabidopsis_thaliana.TAIR10.57.gtf \
    -o /mnt/e/plantgenomics/counts/all_samples_counts.txt \
    -t exon \
    -g gene_id \
    -T 4 \
    /mnt/e/plantgenomics/aligned/SRR7779219.bam \
    /mnt/e/plantgenomics/aligned/SRR7779220.bam \
    /mnt/e/plantgenomics/aligned/SRR7779221.bam \
    /mnt/e/plantgenomics/aligned/SRR7779225.bam \
    /mnt/e/plantgenomics/aligned/SRR7779226.bam \
    /mnt/e/plantgenomics/aligned/SRR7779227.bam

echo "featureCounts complete!"
