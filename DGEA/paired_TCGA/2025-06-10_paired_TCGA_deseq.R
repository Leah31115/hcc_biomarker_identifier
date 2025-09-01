library("DESeq2")
library(ggplot2)
library(ggrepel)
library("vsn")
library("RColorBrewer")
library("pheatmap")
library(tidyverse)
library("magrittr") # for releveling

# Clear environment
rm(list=ls())

count_data <- read.csv("C:/Users/leahe/Documents/LIFE4137_independent_project/xena/paired_TCGA/SNT_cancer_pairs_rawcounts.csv", row.names = 1, check.names = FALSE)
condition <- read.csv("C:/Users/leahe/Documents/LIFE4137_independent_project/xena/paired_TCGA/paired_TCGA_col_data.csv", header = TRUE, row.names = 1)

# Set the row names in condition as the column names in count_data
colnames(count_data) <- row.names(condition)
# Make values integers
count_data <- round(count_data, 0)
# Make dds object for DESeq
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = condition,
                              design = ~sample+condition)
resultsNames(dds) # lists the coefficients
levels(dds$sample)
levels(dds$condition)

# Relevel so SNT is the control i.e Solid_tissue_normal_vs_HCC -> HCC_vs_Solid_tissue_normal)
# note - switching the order these were in the colData.csv file did not fix this problem.
# This is because R orders based on the alphabet!
dds$condition %<>% relevel("Solid_tissue_normal")
levels(dds$condition)

dds <- DESeq(dds)

# Differential expression analysis
res <- results(dds, name="condition_HCC_vs_Solid_tissue_normal") # y (HCC) vs x control (SNT)


# P-values and adj-pvalues
resOrdered <- res[order(res$padj),]
summary(res) # Summary of adjusted pvalues FOR ALL SAMPLES


# PCA with points coloured by sample and shaped by condition
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "sample"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=sample, shape = condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw(base_size = 14) +
  geom_text_repel(aes(label=colnames(assay(vsd))), size=3, force=30, color="black", max.overlaps = 20)

# PCA with points coloured by condition
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()
