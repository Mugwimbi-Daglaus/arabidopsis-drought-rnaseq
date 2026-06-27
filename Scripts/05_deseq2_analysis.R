# DESeq2 Differential Expression Analysis
# Arabidopsis Drought Stress

library(DESeq2)

# Read counts
counts <- read.table("/mnt/e/plantgenomics/counts/all_samples_counts.txt",
                     header=TRUE, skip=1, row.names=1)
counts <- counts[,6:11]
colnames(counts) <- c("Drought1","Drought2","Drought3",
                      "Control1","Control2","Control3")

# Sample info
coldata <- data.frame(
  condition = factor(c("Drought","Drought","Drought",
                       "Control","Control","Control")),
  row.names = colnames(counts))

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=coldata,
                              design=~condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","Drought","Control"))

# Summary
summary(res)

# Save results
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(res),
          file="/mnt/e/plantgenomics/results/DESeq2_drought_results.csv")
write.csv(as.data.frame(sig_genes),
          file="/mnt/e/plantgenomics/results/DESeq2_significant_genes.csv")

cat("Significant genes:", nrow(sig_genes), "\n")

# Volcano plot
library(ggplot2)
res_df <- as.data.frame(res)
res_df$significant <- ifelse(res_df$padj < 0.05 &
                             abs(res_df$log2FoldChange) > 1,
                             "Significant", "Not significant")
ggplot(res_df, aes(x=log2FoldChange, y=-log10(pvalue),
                   color=significant)) +
    geom_point(alpha=0.5, size=0.8) +
    scale_color_manual(values=c("grey","red")) +
    geom_vline(xintercept=c(-1,1), linetype="dashed") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    labs(title="Volcano Plot: Arabidopsis Drought vs Control",
         x="Log2 Fold Change", y="-Log10 P-value") +
    theme_minimal()
ggsave("/mnt/e/plantgenomics/results/volcano_plot.png",
       width=8, height=6, dpi=300)

# PCA plot
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition") +
    labs(title="PCA Plot: Drought vs Control") +
    theme_minimal()
ggsave("/mnt/e/plantgenomics/results/PCA_plot.png",
       width=8, height=6, dpi=300)

# Heatmap
library(pheatmap)
top_genes <- sig_genes[order(sig_genes$padj),]
top50_names <- rownames(head(as.data.frame(top_genes), 50))
vsd_matrix <- assay(vsd)[top50_names,]
pheatmap(vsd_matrix,
         cluster_rows=TRUE, cluster_cols=TRUE,
         show_rownames=TRUE,
         annotation_col=data.frame(condition=coldata$condition,
                                   row.names=rownames(coldata)),
         main="Top 50 Drought-Responsive Genes",
         filename="/mnt/e/plantgenomics/results/heatmap_top50.png",
         width=10, height=12)

# GO enrichment
library(org.At.tair.db)
up_genes <- rownames(subset(res, padj < 0.05 & log2FoldChange > 1))
go_annotations <- select(org.At.tair.db,
                         keys = up_genes,
                         columns = c("GO","ONTOLOGY","GENENAME"),
                         keytype = "TAIR")
bp_terms <- go_annotations[!is.na(go_annotations$ONTOLOGY) &
                           go_annotations$ONTOLOGY == "BP",]
top_terms <- sort(table(bp_terms$GO), decreasing=TRUE)
head(top_terms, 20)
write.csv(as.data.frame(top_terms),
          file="/mnt/e/plantgenomics/results/GO_enrichment_results.csv")

cat("Analysis complete! All results saved.\n")
