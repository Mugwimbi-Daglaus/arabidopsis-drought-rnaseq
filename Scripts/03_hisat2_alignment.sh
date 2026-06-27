#!/bin/bash
# Step 3: Build index and align

# Build genome index
hisat2-build \
    /mnt/e/plantgenomics/genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
    /mnt/e/plantgenomics/genome/tair10_index

# Align all samples
for sample in SRR7779219 SRR7779220 SRR7779221 SRR7779225 SRR7779226 SRR7779227; do
    hisat2 -x /mnt/e/plantgenomics/genome/tair10_index \
        -U /mnt/e/plantgenomics/trimmed/${sample}_trimmed.fastq.gz \
        -S /mnt/e/plantgenomics/aligned/${sample}.sam \
        --summary-file /mnt/e/plantgenomics/aligned/${sample}_stats.txt \
        -p 4
    echo "Done aligning: ${sample}"
done

echo "Alignment complete!"
