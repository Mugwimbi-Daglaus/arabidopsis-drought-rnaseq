#!/bin/bash
# Step 1: Download samples and run FastQC

cd /mnt/e/plantgenomics/raw_data

# Download all samples
for sample in SRR7779219 SRR7779220 SRR7779221 SRR7779225 SRR7779226 SRR7779227; do
    prefetch ${sample}
    fastq-dump --gzip ${sample}
    echo "Downloaded: ${sample}"
done

# FastQC on raw data
mkdir -p /mnt/e/plantgenomics/results/fastqc_raw
fastqc /mnt/e/plantgenomics/raw_data/*.fastq.gz \
       -o /mnt/e/plantgenomics/results/fastqc_raw/ \
       -t 6

echo "FastQC on raw data complete!"
