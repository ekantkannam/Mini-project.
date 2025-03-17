library(DESeq2)
library(tidyverse)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(apeglm)
library(RColorBrewer)
library(ggplot2)
library(readxl)
library(writexl)
library(airway)
setwd('C:/Users/OneDrive/Documents/ngs')
library(readr)
count_blood <- read_tsv("GSE154881_raw_counts_GRCh38.p13_NCBI.tsv", show_col_types = FALSE)
# Read the TSV file
data <- read.delim("GSE154881_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t", header = TRUE)

# Write to CSV file
write.csv(data, "count_blood.csv", row.names = FALSE)
head(data)
dim(data)

# metadata data frame
metadata <- data.frame(
  Accession = c("GSM4681806", "GSM4681807", "GSM4681808", "GSM4681809", "GSM4681810", "GSM4681811"),
  Phenotypes = c("Healthy", "Healthy", "Healthy", "T2D", "T2D", "T2D")
)

# Write to a TSV file
write.table(metadata, "metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
sample_information <- read.csv('metadata.tsv', sep = '\t', row.names = 1)
sample_information
dim(sample_information)

#processing data
count_blood<-data[1:7]
view(count_blood)
count_blood <- count_blood[rowSums(count_blood[-1]!=0)>0,]
view(count_blood)

dim(count_blood)
count_blood <- as.data.frame(count_blood)
rownames(count_blood) <- count_blood[[1]]
count_blood <- count_blood[,-1]
view(count_blood)

factors <- factor(sample_information$Phenotypes)
Phenotypes <- unique(sample_information$Phenotypes)
sample_information$Phenotypes
Phenotypes

sample_information$Phenotypes <- factors
sample_information$Phenotypes
dds <- DESeqDataSetFromMatrix(countData = count_blood, colData = sample_information, design = ~Phenotypes)
dds$Phenotypes <- relevel(dds$Phenotypes, ref = "Healthy")

deg <- DESeq(dds)
res <- results(deg)
res
res <- as.data.frame(res)
class(res)
head(res)
dim(res)
names(res)
res$GeneName <- row.names(res)
names(res)
head(res)

res <- subset(res,
              select = c("GeneName", "padj", "pvalue", "lfcSE", "stat", "log2FoldChange", "baseMean"))
names(res)
write.table(res, file = "result.all.tsv", row.names = FALSE, sep = '\t')

summary(res)
d <- subset(res, padj<0.05 & abs(log2FoldChange) >= 1)
dim(d)
dim(res)

keep <- rowSums(counts(dds) >= 10) >= 3  # At least 10 counts in 3 samples
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)


#DispersionPlot
plotDispEsts(dds, main = "DISPERSION PLOT OF BLOOD")

#volcanoPlot
library(ggplot2)
ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1)) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(padj)") +
  theme(legend.position = "none")

#or
ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), size = 1.5) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(padj)")

#Enhanced volcano plot
ressig <- subset(res, pvalue < 0.05)
summary(ressig)
head(ressig[order(ressig$log2FoldChange),])
head(ressig[order(-ressig$log2FoldChange),])
topGene <- rownames(res)[which.max(res$log2FoldChange)]
plotCounts(dds, gene = topGene, intgroup = c("Phenotypes"))
EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue', pCutoff = 10e-2, FCcutoff = 1.00, xlim = c(-10, 10), ylim = c(0, 30), pointSize = 1.3, labSize = 2.6, max.overlaps = 100, title = "Healthy vs T2D", subtitle = "Differential Expression Analysis", caption = "log2fc cut cutoff=1.3333; pValue cutoff = 10e-5", legendPosition =  "right", legendLabels =  c("NS", "Log2FC", "NA", "p-value and log2FC"), col = c('gray', 'green','blue', 'red'), colAlpha =  0.6, drawConnectors =  FALSE, hline = c(10e-8), widthConnectors = 0.5)

#MA plot
plotMA(res)

# Order results by adjusted p-value (padj)
res_sorted <- res[order(res$padj), ]

# View top 10 most significant genes
head(res_sorted, 10)
sig_genes <- subset(res_sorted, padj < 0.05)

# Filter upregulated genes (logFC > 0, padj < 0.05)
upregulated_genes <- res[which(res$log2FoldChange > 0 & res$padj < 0.05), ]


# Filter downregulated genes (logFC < 0, padj < 0.05)
downregulated_genes <- res[which(res$log2FoldChange < 0 & res$padj < 0.05), ]

# View upregulated genes
upregulated_genes
upregulated_genes <- as.data.frame(upregulated_genes)
write.csv(upregulated_genes, "upregulated_genes.csv", row.names = TRUE)

# View downregulated genes
downregulated_genes
downregulated_genes <- as.data.frame(downregulated_genes)
write.csv(downregulated_genes, "downregulated_genes.csv", row.names = TRUE)

#normalization
norm <- vst(dds)

#heatmap
sampleDist <- dist(t(assay(norm)))
sampleDistMatrix <- as.matrix(sampleDist)
rownames(sampleDistMatrix) <- paste(norm$Phenotypes)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDist, clustering_distance_cols = sampleDist, col = colors)

#PCA plot
plotPCA(norm, intgroup = c("Phenotypes"))

res_df <- as.data.frame(res)
library(AnnotationDbi)
library(org.Hs.eg.db)
res_df$ensembl_id <- mapIds(org.Hs.eg.db,
                            keys = rownames(res_df),
                            column = "ENSEMBL",      # Map to Ensembl IDs
                            keytype = "ENTREZID",      # Replace 'SYMBOL' with your gene ID format
                            multiVals = "first")

res_df$gene_name <- mapIds(org.Hs.eg.db,
                           keys = rownames(res_df),
                           column = "GENENAME",     # Map to Gene Names (descriptions)
                           keytype = "ENTREZID",      # Replace 'SYMBOL' with your input ID type
                           multiVals = "first")

res_df$gene_symbol <- mapIds(org.Hs.eg.db,
                             keys = rownames(res_df),
                             column = "SYMBOL",      # Map to Gene Symbols
                             keytype = "ENTREZID",     # Replace 'SYMBOL' with the current gene ID type
                             multiVals = "first")

head(res_df)
write.csv(res_df, "res.csv", row.names = TRUE)


library(clusterProfiler)
library(enrichplot)  # For visualization
library(org.Hs.eg.db)  # Annotation database for humans
library(ggplot2)

gene_list <- rownames(res_df)

# Biological Process (BP)
ego_bp <- enrichGO(gene = gene_list, 
                   OrgDb = org.Hs.eg.db, 
                   keyType = "ENTREZID", 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.05)

# Molecular Function (MF)
ego_mf <- enrichGO(gene = gene_list, 
                   OrgDb = org.Hs.eg.db, 
                   keyType = "ENTREZID", 
                   ont = "MF", 
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.05)

# Cellular Component (CC)
ego_cc <- enrichGO(gene = gene_list, 
                   OrgDb = org.Hs.eg.db, 
                   keyType = "ENTREZID", 
                   ont = "CC", 
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.05)

#dot plots
dotplot(ego_bp, showCategory = 20) + 
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  ggtitle("Biological Process Enrichment")

dotplot(ego_mf, showCategory = 20) + 
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  ggtitle("Molecular Function")

dotplot(ego_cc, showCategory = 20) + 
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  ggtitle("Cellular Component")


#or
# Biological Process Dot Plot
dotplot(ego_bp, showCategory = 15, font.size = 12, label_format = 50) +
  ggtitle("GO Enrichment - Biological Process")

# Molecular Function Dot Plot
dotplot(ego_mf, showCategory = 15, font.size = 12, label_format = 50) +
  ggtitle("GO Enrichment - Molecular Function")

# Cellular Component Dot Plot
dotplot(ego_cc, showCategory = 15, font.size = 12, label_format = 50) +
  ggtitle("GO Enrichment - Cellular Component")
