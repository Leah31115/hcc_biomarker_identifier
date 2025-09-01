library("DESeq2")
library(ggplot2)
library(ggrepel)
library("vsn")
library("RColorBrewer")
library("pheatmap")
library(tidyverse)

# GO ontology
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Clear environment
rm(list=ls())

count_data <- read.csv("C:/Users/leahe/Documents/LIFE4137_independent_project/xena/clean_script/s1s3/stage1_stage3_rawcounts.csv", row.names = 1)
condition <- read.csv("C:/Users/leahe/Documents/LIFE4137_independent_project/xena/clean_script/s1s3/s1s3_raw_col_data.csv", header = TRUE, row.names = 1)

# Set the row names in condition as the column names in count_data
colnames(count_data) <- row.names(condition)
# Make values integers
count_data <- round(count_data, 0)
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = condition,
                              design = ~ condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

# Differential expression analysis
res <- results(dds, name="condition_Stage3_HCC_vs_Stage1_HCC") # y experiment (stage3) against x control (stage1)

# P-values and adj-pvalues
resOrdered <- res[order(res$padj),]
summary(res) # Summary of adjusted pvalues FOR ALL SAMPLES

# padj < 0.01
sum(res$padj < 0.01, na.rm=TRUE) # 10
# padj < 0.05
sum(res$padj < 0.05, na.rm=TRUE) # 100

# MA plot 
plotMA(res, alpha = 0.01, ylim=c(-5, 5))
plotMA(res, alpha = 0.05, ylim=c(-3, 3))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="limegreen", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=1, col="limegreen", cex = 0.8)
})

# Log counts of the most significant gene
clean_d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                      returnData=TRUE)

ggplot(clean_d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10() +
  ylab("Log count") +
  xlab("Condition") +
  theme_bw()

# Top differentially expressed gene
which.min(res$padj) # 2368 - ENSG00000100652.5 gene
stage1_gene_count <- clean_d %>%
  filter(condition=="Stage1_HCC")
mean(stage1_gene_count[["count"]]) # 12.02866
median(stage1_gene_count[["count"]])  # 12.60543

stage3_gene_count <- clean_d %>%
  filter(condition=="Stage3_HCC")
mean(stage3_gene_count[["count"]]) # 9.644948
median(stage3_gene_count[["count"]]) # 9.787448


# Heatmap
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
head(assay(vsd), 3)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(assay(vsd)))
colnames(sampleDistMatrix) <- paste(colnames(assay(vsd)))
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# PCA
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()

# Volcano plot
# Subset data
volcano <- data.frame(res@rownames, res$padj, res$log2FoldChange)
# Rename columns
colnames(volcano) <- c("Gene", "padj", "log2FoldChange")
volcano <- na.omit(volcano) # remove rows with NAs
volcano$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
volcano$diffexpressed[volcano$log2FoldChange > 0.6 & volcano$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
volcano$diffexpressed[volcano$log2FoldChange < -0.6 & volcano$padj < 0.05] <- "DOWN"


# Label the top 5 up and downregulated DEGs
top_genes <- c() # empty vector for gene annotations to go into
# Get the top 5 most significant up-regulated genes (smallest p-values)
upreg <- volcano %>%
  filter(diffexpressed == "UP")
sig_upreg_genes <- upreg %>%
  arrange(desc(-padj))
for (idx in 1:5) {
  row <- sig_upreg_genes[idx, "Gene"]
  print(row)
  top_genes <- append(top_genes, row)
}
top_genes

# Get the top 5 most significant down-regulated genes (smallest p-values)
downreg <- volcano %>%
  filter(diffexpressed == "DOWN")
sig_downreg_genes <- downreg %>%
  arrange(desc(-padj))
for (idx in 1:5) {
  row <- sig_downreg_genes[idx, "Gene"]
  print(row)
  top_genes <- append(top_genes, row)
}
top_genes

# Label the top 10 DEGs
volcano$label <- ifelse(volcano$Gene %in% top_genes, volcano$Gene, NA)


# Volcano plot
ggplot(data=volcano, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=label)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("cornflowerblue", "grey", "firebrick2"),
                     labels=c("Downregulated", "No", "Upregulated")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="grey25") +
  geom_hline(yintercept=-log10(0.05), col="grey25") +
  ylim(-1, 3) +
  xlim(-3, 3) +
  geom_text_repel(max.overlaps = 14, size = 3, force = 20, col = "black")

# Density plot
ggplot(volcano, aes(x=log2FoldChange, fill=diffexpressed)) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values = c("cornflowerblue", "grey", "firebrick2"), name="Differential Gene Expression") +
  theme_classic() +
  labs(title="Density plot of Log2 Fold Changes", x="Log2 Fold Change", y="Density", legend="He") +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black",
                                fill = NA, size = 0.8)
  ) +
  xlim(-2, 2)


# GO ontology 
# Upregulated pathways - NA
upreg_genes <- as.vector(upreg$Gene)
# Remove gene versions
upreg_genes <- substr(upreg_genes, 1, 15)

upreg_GO <- enrichGO(gene = upreg_genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
as.data.frame(upreg_GO)
fit <- plot(barplot(upreg_GO, showCategory = 10))
fit
write.csv(as.data.frame(upreg_GO),file="C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/GO_analysis/upreg/s1s3_upreg.csv")

# Downregulated pathways - NA
downreg_genes <- as.vector(downreg$Gene)
# Remove gene versions
downreg_genes <- substr(downreg_genes, 1, 15)

downreg_GO <- enrichGO(gene = downreg_genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
as.data.frame(downreg_GO)
fit <- plot(barplot(downreg_GO, showCategory = 10))
fit
write.csv(as.data.frame(downreg_GO),file="C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/GO_analysis/downreg/s1s3_downreg.csv")



# Upregulated, genes where padj < 0.05
upreg_05 <- upreg %>%
  arrange(desc(log2FoldChange))
write.table(upreg_05["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/clean_script/s1s3/s1s3_upreg_05_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)

# Top 10 upregulated significant genes (padj < 0.05)
upreg_10 <- upreg %>%
  arrange(desc(-padj))
upreg_10 <- upreg_10[1:10, ]
write.table(upreg_10["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/clean_script/s1s3/s1s3_upreg_10_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)


# Top 50 upregulated significant genes (padj < 0.05)
upreg_50 <- upreg %>%
  arrange(desc(-padj))
upreg_50 <- upreg_50[1:50, ]
write.table(upreg_50["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/clean_script/s1s3/s1s3_upreg_50_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)


# Summary of the upregulated significant genes (padj < 0.05)
upreg_sum <- upreg_10 %>%
  dpylr::select(Gene, padj, log2FoldChange)
write.table(upreg_sum,
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/clean_script/s1s3/s1s3_upreg_summary.txt",
            sep = ",",
            row.names = FALSE,
            col.names = FALSE
)


# downregulated, genes where padj < 0.05
downreg_sum <- downreg %>%
  arrange(desc(-padj)) %>%
  dpylr::select(Gene, padj, log2FoldChange)
write.table(downreg_sum,
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/clean_script/s1s3/s1s3_downreg_summary.txt",
            sep = ",",
            row.names = FALSE,
            col.names = FALSE
)




# Pvalues for the top 5 up and down significant genes across all DGEA comparisons (including stage vs stage and paired)
top_genes <- c("ENSG00000254772.10",
               "ENSG00000280614.1",
               "ENSG00000281181.1",
               "ENSG00000136305.11",
               "ENSG00000178922.17",
               "ENSG00000210196.2",
               "ENSG00000199753.1",
               "ENSG00000256393.1",
               "ENSG00000231131.8",
               "ENSG00000250982.2",
               "ENSG00000240801.1",
               "ENSG00000213793.5",
               "ENSG00000099977.14",
               "ENSG00000237039.1",
               "ENSG00000227063.5",
               "ENSG00000243176.6",
               "ENSG00000227077.3",
               "ENSG00000136149.6",
               "ENSG00000011052.21",
               "ENSG00000058673.17",
               "ENSG00000105583.11",
               "ENSG00000129965.15",
               "ENSG00000287928.1",
               "ENSG00000288048.1",
               "ENSG00000234389.1",
               "ENSG00000273209.1",
               "ENSG00000176040.13",
               "ENSG00000213867.4",
               "ENSG00000236824.2",
               "ENSG00000130822.16",
               "ENSG00000159337.7",
               "ENSG00000214510.10",
               "ENSG00000282301.3",
               "ENSG00000164093.17",
               "ENSG00000181577.16",
               "ENSG00000107159.13",
               "ENSG00000205358.4",
               "ENSG00000125144.14",
               "ENSG00000019169.10",
               "ENSG00000160339.16",
               "ENSG00000198417.7",
               "ENSG00000175063.17",
               "ENSG00000164283.13",
               "ENSG00000147257.15",
               "ENSG00000204291.11",
               "ENSG00000117399.14",
               "ENSG00000164611.13",
               "ENSG00000101412.13",
               "ENSG00000089685.15",
               "ENSG00000145708.11"
)

top_genes_data <- volcano %>%
  dplyr::select(Gene, padj, diffexpressed) %>%
  filter(Gene %in% top_genes) %>%
  arrange(factor(Gene, levels = top_genes)) # Order genes based on vector order

# Save data
write.table(top_genes_data,
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/clean_script/s1s3/all_top5_dgea_genes.txt",
            sep = ",",
            row.names = FALSE,
            col.names = TRUE
)
