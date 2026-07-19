# ================================================
# Arabidopsis thaliana Drought Stress RNA-seq Analysis
# Dataset: GSE119382
# Conditions: Drought vs Control
# Environment: conda activate plantgenomics → R
# ================================================

# ---- Load libraries ----
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(WGCNA)
library(ggraph)
library(igraph)
library(clusterProfiler)
library(org.At.tair.db)

# Fix cor conflict
cor <- WGCNA::cor

# ---- SECTION 1: DESeq2 ----
counts <- read.table("/mnt/e/arab/counts/all_samples_counts.txt",
                     header=TRUE, skip=1, row.names=1)
counts <- counts[,6:11]
colnames(counts) <- c("Drought1","Drought2","Drought3",
                      "Control1","Control2","Control3")

coldata <- data.frame(
  condition = factor(c("Drought","Drought","Drought",
                       "Control","Control","Control")),
  row.names = colnames(counts))

dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=coldata,
                              design=~condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","Drought","Control"))
summary(res)

# Save results
write.csv(as.data.frame(res),
          "/mnt/e/arab/results/DESeq2_drought_results.csv")
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(sig_genes),
          "/mnt/e/arab/results/DESeq2_significant_genes.csv")
cat("Significant genes:", nrow(sig_genes), "\n")
cat("Upregulated:", nrow(subset(sig_genes, log2FoldChange > 1)), "\n")
cat("Downregulated:", nrow(subset(sig_genes, log2FoldChange < -1)), "\n")

# ---- SECTION 2: Visualization ----
vsd <- vst(dds, blind=FALSE)

# Volcano plot
res_df <- as.data.frame(res)
res_df$significant <- ifelse(res_df$padj < 0.05 &
                             abs(res_df$log2FoldChange) > 1,
                             "Significant","Not significant")
ggplot(res_df, aes(x=log2FoldChange, y=-log10(pvalue),
                   color=significant)) +
    geom_point(alpha=0.5, size=0.8) +
    scale_color_manual(values=c("grey","red")) +
    geom_vline(xintercept=c(-1,1), linetype="dashed") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    labs(title="Volcano Plot: Arabidopsis Drought vs Control",
         x="Log2 Fold Change", y="-Log10 P-value") +
    theme_minimal()
ggsave("/mnt/e/arab/results/volcano_plot.png",
       width=8, height=6, dpi=300)

# PCA plot
plotPCA(vsd, intgroup="condition") +
    labs(title="PCA Plot: Drought vs Control") +
    theme_minimal()
ggsave("/mnt/e/arab/results/PCA_plot.png",
       width=8, height=6, dpi=300)

# Heatmap top 50# ================================================
# Arabidopsis thaliana Drought Stress RNA-seq Analysis
# Dataset: GSE119382
# Conditions: Drought vs Control
# Environment: conda activate plantgenomics → R
# ================================================

# ---- Load libraries ----
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(WGCNA)
library(ggraph)
library(igraph)
library(clusterProfiler)
library(org.At.tair.db)

# Fix cor conflict
cor <- WGCNA::cor

# ---- SECTION 1: DESeq2 ----
counts <- read.table("/mnt/e/arab/counts/all_samples_counts.txt",
                     header=TRUE, skip=1, row.names=1)
counts <- counts[,6:11]
colnames(counts) <- c("Drought1","Drought2","Drought3",
                      "Control1","Control2","Control3")

coldata <- data.frame(
  condition = factor(c("Drought","Drought","Drought",
                       "Control","Control","Control")),
  row.names = colnames(counts))

dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=coldata,
                              design=~condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","Drought","Control"))
summary(res)

# Save results
write.csv(as.data.frame(res),
          "/mnt/e/arab/results/DESeq2_drought_results.csv")
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(sig_genes),
          "/mnt/e/arab/results/DESeq2_significant_genes.csv")
cat("Significant genes:", nrow(sig_genes), "\n")
cat("Upregulated:", nrow(subset(sig_genes, log2FoldChange > 1)), "\n")
cat("Downregulated:", nrow(subset(sig_genes, log2FoldChange < -1)), "\n")

# ---- SECTION 2: Visualization ----
vsd <- vst(dds, blind=FALSE)

# Volcano plot
res_df <- as.data.frame(res)
res_df$significant <- ifelse(res_df$padj < 0.05 &
                               abs(res_df$log2FoldChange) > 1,
                             "Significant","Not significant")
ggplot(res_df, aes(x=log2FoldChange, y=-log10(pvalue),
                   color=significant)) +
  geom_point(alpha=0.5, size=0.8) +
  scale_color_manual(values=c("grey","red")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  labs(title="Volcano Plot: Arabidopsis Drought vs Control",
       x="Log2 Fold Change", y="-Log10 P-value") +
  theme_minimal()
ggsave("/mnt/e/arab/results/volcano_plot.png",
       width=8, height=6, dpi=300)

# PCA plot
plotPCA(vsd, intgroup="condition") +
  labs(title="PCA Plot: Drought vs Control") +
  theme_minimal()
ggsave("/mnt/e/arab/results/PCA_plot.png",
       width=8, height=6, dpi=300)

# Heatmap top 50
top50 <- rownames(head(sig_genes[order(sig_genes$padj),], 50))
vsd_matrix <- assay(vsd)[top50,]
pheatmap(vsd_matrix,
         cluster_rows=TRUE, cluster_cols=TRUE,
         scale="row", show_rownames=TRUE,
         annotation_col=data.frame(condition=coldata$condition,
                                   row.names=rownames(coldata)),
         main="Top 50 Drought-Responsive Genes",
         filename="/mnt/e/arab/results/heatmap_top50.png",
         width=10, height=12)

# ---- SECTION 3: GO Enrichment ----
up_genes <- rownames(subset(res, padj < 0.05 & log2FoldChange > 1))
go_annotations <- select(org.At.tair.db,
                         keys=up_genes,
                         columns=c("GO","ONTOLOGY","GENENAME"),
                         keytype="TAIR")
bp_terms <- go_annotations[!is.na(go_annotations$ONTOLOGY) &
                             go_annotations$ONTOLOGY=="BP",]
top_terms <- sort(table(bp_terms$GO), decreasing=TRUE)
write.csv(as.data.frame(top_terms),
          "/mnt/e/arab/results/GO_enrichment_results.csv")

# ---- SECTION 4: WGCNA ----
allowWGCNAThreads()
datExpr <- t(assay(vsd))
gsg <- goodSamplesGenes(datExpr, verbose=3)
if(!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

net <- blockwiseModules(datExpr, power=18,
                        TOMType="unsigned", minModuleSize=30,
                        reassignThreshold=0, mergeCutHeight=0.25,
                        numericLabels=TRUE, pamRespectsDendro=FALSE,
                        maxBlockSize=5000, saveTOMs=FALSE, verbose=3)

moduleColors <- labels2colors(net$colors)
cat("Modules found:", length(table(moduleColors)), "\n")

# Module-trait correlation
traitData <- data.frame(Drought=c(1,1,1,0,0,0),
                        row.names=rownames(datExpr))
MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
moduleTraitCor <- cor(MEs, traitData, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Hub genes
brown_genes <- names(net$colors)[moduleColors=="brown"]
turquoise_genes <- names(net$colors)[moduleColors=="turquoise"]
MMbrown <- cor(datExpr, MEs[,"MEbrown"], use="p")
MMturquoise <- cor(datExpr, MEs[,"MEturquoise"], use="p")
top_brown <- names(sort(abs(MMbrown[brown_genes,]), decreasing=TRUE))[1:10]
top_turquoise <- names(sort(abs(MMturquoise[turquoise_genes,]), decreasing=TRUE))[1:10]

# Save WGCNA results
write.csv(data.frame(gene=top_brown),
          "/mnt/e/arab/results/WGCNA_brown_hub_genes.csv")
write.csv(data.frame(gene=top_turquoise),
          "/mnt/e/arab/results/WGCNA_turquoise_hub_genes.csv")
write.csv(as.data.frame(moduleTraitCor),
          "/mnt/e/arab/results/WGCNA_module_trait_correlations.csv")
write.csv(data.frame(gene=names(moduleColors), module=moduleColors),
          "/mnt/e/arab/results/WGCNA_module_assignments.csv")

# Module-trait heatmap
png("/mnt/e/arab/results/WGCNA_module_trait_heatmap.png",
    width=800, height=1200, res=150)
pheatmap(as.matrix(moduleTraitCor),
         main="Module-Trait Relationships",
         display_numbers=round(moduleTraitCor,2),
         color=colorRampPalette(c("blue","white","red"))(50),
         cluster_rows=TRUE, cluster_cols=FALSE)
dev.off()

# Hub gene heatmap
all_hub_genes <- c(top_brown, top_turquoise)
hub_matrix <- assay(vsd)[all_hub_genes,]
pheatmap(hub_matrix,
         cluster_rows=FALSE, cluster_cols=FALSE, scale="row",
         show_rownames=TRUE,
         annotation_col=data.frame(condition=coldata$condition,
                                   row.names=rownames(coldata)),
         main="Hub Gene Expression: Brown & Turquoise Modules",
         filename="/mnt/e/arab/results/WGCNA_hub_gene_heatmap.png",
         width=10, height=8)

# ---- SECTION 5: Hub gene network with ggraph ----
hub_expr <- t(assay(vsd)[all_hub_genes,])
adjacency_hub <- adjacency(hub_expr, power=18, type="unsigned")
adj_matrix <- as.matrix(adjacency_hub)
diag(adj_matrix) <- 0
adj_matrix[adj_matrix < 0.1] <- 0

g <- graph_from_adjacency_matrix(adj_matrix,
                                 mode="undirected",
                                 weighted=TRUE,
                                 diag=FALSE)

V(g)$module <- ifelse(V(g)$name %in% top_brown,
                      "Brown module", "Turquoise module")
V(g)$degree <- degree(g)

p <- ggraph(g, layout="fr") +
  geom_edge_link(aes(alpha=weight), color="grey50",
                 width=1, show.legend=FALSE) +
  scale_edge_alpha(range=c(0.3, 0.9)) +
  geom_node_point(aes(color=module, size=degree)) +
  scale_color_manual(values=c("Brown module"="#8B4513",
                              "Turquoise module"="#40E0D0"),
                     name="WGCNA Module") +
  scale_size_continuous(range=c(8,16), name="Node Connectivity") +
  geom_node_text(aes(label=name), size=2.8, fontface="bold",
                 repel=TRUE, max.overlaps=30,
                 bg.color="white", bg.r=0.1) +
  theme_graph(base_family="sans") +
  theme(legend.position="bottom",
        plot.title=element_text(face="bold", size=15, hjust=0.5),
        plot.subtitle=element_text(size=11, hjust=0.5, color="grey40")) +
  labs(title="Hub Gene Co-expression Network",
       subtitle="Arabidopsis thaliana Drought Stress Response (WGCNA)")

ggsave("/mnt/e/arab/results/hub_gene_network.png",
       plot=p, width=12, height=10, dpi=300, bg="white")

cat("====== Analysis Complete! ======\n")
cat("All results saved to: /mnt/e/arab/results/\n")

top50 <- rownames(head(sig_genes[order(sig_genes$padj),], 50))
vsd_matrix <- assay(vsd)[top50,]
pheatmap(vsd_matrix,
         cluster_rows=TRUE, cluster_cols=TRUE,
         scale="row", show_rownames=TRUE,
         annotation_col=data.frame(condition=coldata$condition,
                                   row.names=rownames(coldata)),
         main="Top 50 Drought-Responsive Genes",
         filename="/mnt/e/arab/results/heatmap_top50.png",
         width=10, height=12)

# ---- SECTION 3: GO Enrichment ----
up_genes <- rownames(subset(res, padj < 0.05 & log2FoldChange > 1))
go_annotations <- select(org.At.tair.db,
                         keys=up_genes,
                         columns=c("GO","ONTOLOGY","GENENAME"),
                         keytype="TAIR")
bp_terms <- go_annotations[!is.na(go_annotations$ONTOLOGY) &
                           go_annotations$ONTOLOGY=="BP",]
top_terms <- sort(table(bp_terms$GO), decreasing=TRUE)
write.csv(as.data.frame(top_terms),
          "/mnt/e/arab/results/GO_enrichment_results.csv")

# ---- SECTION 4: WGCNA ----
allowWGCNAThreads()
datExpr <- t(assay(vsd))
gsg <- goodSamplesGenes(datExpr, verbose=3)
if(!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

net <- blockwiseModules(datExpr, power=18,
                        TOMType="unsigned", minModuleSize=30,
                        reassignThreshold=0, mergeCutHeight=0.25,
                        numericLabels=TRUE, pamRespectsDendro=FALSE,
                        maxBlockSize=5000, saveTOMs=FALSE, verbose=3)

moduleColors <- labels2colors(net$colors)
cat("Modules found:", length(table(moduleColors)), "\n")

# Module-trait correlation
traitData <- data.frame(Drought=c(1,1,1,0,0,0),
                        row.names=rownames(datExpr))
MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
moduleTraitCor <- cor(MEs, traitData, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Hub genes
brown_genes <- names(net$colors)[moduleColors=="brown"]
turquoise_genes <- names(net$colors)[moduleColors=="turquoise"]
MMbrown <- cor(datExpr, MEs[,"MEbrown"], use="p")
MMturquoise <- cor(datExpr, MEs[,"MEturquoise"], use="p")
top_brown <- names(sort(abs(MMbrown[brown_genes,]), decreasing=TRUE))[1:10]
top_turquoise <- names(sort(abs(MMturquoise[turquoise_genes,]), decreasing=TRUE))[1:10]

# Save WGCNA results
write.csv(data.frame(gene=top_brown),
          "/mnt/e/arab/results/WGCNA_brown_hub_genes.csv")
write.csv(data.frame(gene=top_turquoise),
          "/mnt/e/arab/results/WGCNA_turquoise_hub_genes.csv")
write.csv(as.data.frame(moduleTraitCor),
          "/mnt/e/arab/results/WGCNA_module_trait_correlations.csv")
write.csv(data.frame(gene=names(moduleColors), module=moduleColors),
          "/mnt/e/arab/results/WGCNA_module_assignments.csv")

# Module-trait heatmap
png("/mnt/e/arab/results/WGCNA_module_trait_heatmap.png",
    width=800, height=1200, res=150)
pheatmap(as.matrix(moduleTraitCor),
         main="Module-Trait Relationships",
         display_numbers=round(moduleTraitCor,2),
         color=colorRampPalette(c("blue","white","red"))(50),
         cluster_rows=TRUE, cluster_cols=FALSE)
dev.off()

# Hub gene heatmap
all_hub_genes <- c(top_brown, top_turquoise)
hub_matrix <- assay(vsd)[all_hub_genes,]
pheatmap(hub_matrix,
         cluster_rows=FALSE, cluster_cols=FALSE, scale="row",
         show_rownames=TRUE,
         annotation_col=data.frame(condition=coldata$condition,
                                   row.names=rownames(coldata)),
         main="Hub Gene Expression: Brown & Turquoise Modules",
         filename="/mnt/e/arab/results/WGCNA_hub_gene_heatmap.png",
         width=10, height=8)

# ---- SECTION 5: Hub gene network with ggraph ----
hub_expr <- t(assay(vsd)[all_hub_genes,])
adjacency_hub <- adjacency(hub_expr, power=18, type="unsigned")
adj_matrix <- as.matrix(adjacency_hub)
diag(adj_matrix) <- 0
adj_matrix[adj_matrix < 0.1] <- 0

g <- graph_from_adjacency_matrix(adj_matrix,
                                  mode="undirected",
                                  weighted=TRUE,
                                  diag=FALSE)

V(g)$module <- ifelse(V(g)$name %in% top_brown,
                      "Brown module", "Turquoise module")
V(g)$degree <- degree(g)

p <- ggraph(g, layout="fr") +
    geom_edge_link(aes(alpha=weight), color="grey50",
                   width=1, show.legend=FALSE) +
    scale_edge_alpha(range=c(0.3, 0.9)) +
    geom_node_point(aes(color=module, size=degree)) +
    scale_color_manual(values=c("Brown module"="#8B4513",
                                "Turquoise module"="#40E0D0"),
                       name="WGCNA Module") +
    scale_size_continuous(range=c(8,16), name="Node Connectivity") +
    geom_node_text(aes(label=name), size=2.8, fontface="bold",
                   repel=TRUE, max.overlaps=30,
                   bg.color="white", bg.r=0.1) +
    theme_graph(base_family="sans") +
    theme(legend.position="bottom",
          plot.title=element_text(face="bold", size=15, hjust=0.5),
          plot.subtitle=element_text(size=11, hjust=0.5, color="grey40")) +
    labs(title="Hub Gene Co-expression Network",
         subtitle="Arabidopsis thaliana Drought Stress Response (WGCNA)")

ggsave("/mnt/e/arab/results/hub_gene_network.png",
       plot=p, width=12, height=10, dpi=300, bg="white")

cat("====== Analysis Complete! ======\n")
cat("All results saved to: /mnt/e/arab/results/\n")
