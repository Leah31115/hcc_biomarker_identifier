library("DESeq2")
library(ggplot2)
library(ggrepel)
library("vsn") # transform data according to variance
library("RColorBrewer") # For colouring heatmap
library("pheatmap") # for heatmap
library(tidyverse)

# Clear environment
rm(list=ls())

# Load in data
count_data <- read.csv("C:/Users/leahe/Documents/LIFE4137_independent_project/xena/clean_script/s1s4/stage1_stage4_rawcounts.csv", row.names = 1)
condition <- read.csv("C:/Users/leahe/Documents/LIFE4137_independent_project/xena/clean_script/s1s4/s1s4_raw_col_data.csv", header = TRUE, row.names = 1)

# Set the row names in condition as the column names in count_data
colnames(count_data) <- row.names(condition)
# Make values integers
count_data <- round(count_data, 0)
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = condition,
                              design = ~ condition)
dds <- DESeq(dds)
resultsNames(dds)

# Differential expression analysis
res <- results(dds, name="condition_Stage4_HCC_vs_Stage1_HCC") # y experiment (stage4) against x control (Stage1)

# P-values and adj-pvalues
resOrdered <- res[order(res$padj),]
summary(res) # Summary of adjusted pvalues FOR ALL SAMPLES

# How many adjusted p-values were less than 0.01 FOR ALL SAMPLES?
sum(res$padj < 0.01, na.rm=TRUE) # 0
sum(res$padj < 0.05, na.rm=TRUE) # 0

# MA plot 
plotMA(res, alpha = 0.05, ylim=c(-5, 5))

# Plot to compare log counts of the most significant gene (smallest padj value)
clean_d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                      returnData=TRUE)

ggplot(clean_d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10() +
  ylab("Log count") +
  xlab("Condition") +
  theme_bw()

# Top differentially expressed gene
which.min(res$padj) #1 - ENSG00000000003.15 gene (NOT SIGNIFICANT)
# log count statistics of the top differentially expressed gene
s1_gene_count <- clean_d %>%
  filter(condition=="Stage1_HCC")
mean(s1_gene_count[["count"]]) # 12.50067
median(s1_gene_count[["count"]])  # 12.52906

s4_gene_count <- clean_d %>%
  filter(condition=="Stage4_HCC")
mean(stage4_gene_count[["count"]]) # 12.34442
median(stage4_gene_count[["count"]]) # 12.05003


# Heatmap of the sample-to-sample distances
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

# PCA with customisation
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
volcano <- na.omit(volcano)
volcano$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
volcano$diffexpressed[volcano$log2FoldChange > 0.6 & volcano$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
volcano$diffexpressed[volcano$log2FoldChange < -0.6 & volcano$padj < 0.05] <- "DOWN"

# Label the top 5 up- and down-DEGs
top_genes <- c()
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


# Label the top_genes
volcano$label <- ifelse(volcano$Gene %in% top_genes, volcano$Gene, NA)

# Volcano plot
ggplot(data=volcano, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=label)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("grey", "firebrick2"),
                     labels=c("Not significant", "Upregulated")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  ylim(-5, 35) +
  xlim(-10, 10) +
  geom_text_repel(max.overlaps = 20, force = 60, size = 3, col = "black")

# Density plot
ggplot(volcano, aes(x=log2FoldChange, fill=diffexpressed)) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values = c("grey", "firebrick2"), name="Differential Gene Expression") +
  theme_classic() +
  labs(title="Density plot of Log2 Fold Changes", x="Log2 Fold Change", y="Density", legend="He") +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black",
                                fill = NA, linewidth = 0.8)
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
  select(Gene, padj, diffexpressed) %>%
  filter(Gene %in% top_genes) %>%
  arrange(factor(Gene, levels = top_genes)) # Order genes based on vector order

# Save data
write.table(top_genes_data,
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/clean_script/s1s4/all_top5_dgea_genes.txt",
            sep = ",",
            row.names = FALSE,
            col.names = TRUE
)

