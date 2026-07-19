#!/bin/bash
# Step 5: Count reads per gene
# Environment: conda activate plantgenomics

featureCounts \
    -a /mnt/e/arab/annotation/Arabidopsis_thaliana.TAIR10.57.gtf \
    -o /mnt/e/arab/counts/all_samples_counts.txt \
    -t exon \
    -g gene_id \
    -T 4 \
    /mnt/e/arab/aligned/SRR7779219.bam \
    /mnt/e/arab/aligned/SRR7779220.bam \
    /mnt/e/arab/aligned/SRR7779221.bam \
    /mnt/e/arab/aligned/SRR7779225.bam \
    /mnt/e/arab/aligned/SRR7779226.bam \
    /mnt/e/arab/aligned/SRR7779227.bam

echo "featureCounts complete!"
