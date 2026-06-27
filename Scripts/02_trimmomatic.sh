#!/bin/bash
# Step 2: Trim all samples

for sample in SRR7779219 SRR7779220 SRR7779221 SRR7779225 SRR7779226 SRR7779227; do
    trimmomatic SE \
        /mnt/e/plantgenomics/raw_data/${sample}.fastq.gz \
        /mnt/e/plantgenomics/trimmed/${sample}_trimmed.fastq.gz \
        ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 MINLEN:36
    echo "Done: ${sample}"
done

# FastQC on trimmed data
mkdir -p /mnt/e/plantgenomics/results/fastqc_trimmed
fastqc /mnt/e/plantgenomics/trimmed/*_trimmed.fastq.gz \
       -o /mnt/e/plantgenomics/results/fastqc_trimmed/ \
       -t 6

echo "Trimming complete!"
