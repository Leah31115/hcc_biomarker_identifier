library("DESeq2")
library(ggplot2)
library("vsn")
library("pheatmap")
library(tidyverse)

# Clear environment
rm(list=ls())

count_data <- read.csv("C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/healthy_vs_cancer/healthy_stages_rawcounts.csv", row.names = 1, check.names = FALSE)
condition <- read.csv("C:/Users/leahe/Documents/LIFE4137_independent_project/xena/clean_script/hstages/healthy_stages_raw_col_data.csv", header = TRUE, row.names = 1)
# Set the row names in condition as the column names in count_data
colnames(count_data) <- row.names(condition)
# Make values integers
count_data <- round(count_data, 0)
# Make dds object for DESeq
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = condition,
                              design = ~ condition)

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

# PCA where cancer samples are coloured by stage
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw(base_size = 16)
