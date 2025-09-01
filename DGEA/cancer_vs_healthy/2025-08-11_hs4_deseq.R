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

# Load in sample matrix and coldata 
count_data <- read.csv("C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/hs4/healthy_stage4_rawcounts.csv", row.names = 1)
condition <- read.csv("C:/Users/leahe/Documents/LIFE4137_independent_project/xena/clean_script/hs4/hs4_raw_col_data.csv", header = TRUE, row.names = 1)


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
res <- results(dds, name="condition_Stage4_HCC_vs_Healthy") # y experiment (stage4) against x control (healthy)

# P-values and adj-pvalues
resOrdered <- res[order(res$padj),]
summary(res) # Summary of adjusted pvalues FOR ALL SAMPLES

# Number of genes with padj < 0.01
sum(res$padj < 0.01, na.rm=TRUE) # 1985
# Number of genes with padj < 0.05
sum(res$padj < 0.05, na.rm=TRUE) # 3147

# MA plot 
plotMA(res, alpha = 0.05, ylim=c(-10, 10))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="limegreen", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=4, col="limegreen")
})

# Counts of most significant gene
clean_d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                      returnData=TRUE)

ggplot(clean_d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10() +
  ylab("Log count") +
  xlab("Condition") +
  theme_bw()

# Top differentially expressed gene
which.min(res$padj) # 33895 -  gene
# log count statistics of the top differentially expressed gene
healthy_gene_count <- clean_d %>%
  filter(condition=="Healthy")
mean(healthy_gene_count[["count"]]) # 0.9324913
median(healthy_gene_count[["count"]])  # 0.5

stage4_gene_count <- clean_d %>%
  filter(condition=="Stage4_HCC")
mean(stage4_gene_count[["count"]]) # 8.453062
median(stage4_gene_count[["count"]]) # 8.972071


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

# PCA
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()

# Volcano plot made manually
# Categorise differentially expressed genes
# Subset data
volcano <- data.frame(res@rownames, res$padj, res$log2FoldChange)
# Rename columns
colnames(volcano) <- c("Gene", "padj", "log2FoldChange")
volcano <- na.omit(volcano) # remove rows with NAs
volcano$diffexpressed <- "NO"
volcano

# Up-regulated genes if pvalue < 0.05 and 0.6 < log2Foldchange
volcano$diffexpressed[volcano$log2FoldChange > 0.6 & volcano$padj < 0.05] <- "UP"
# Down-regulated genes if pvalue < 0.05 and log2Foldchange < 0.6
volcano$diffexpressed[volcano$log2FoldChange < -0.6 & volcano$padj < 0.05] <- "DOWN"


# Label the top 5 differentially upregulated and top 5 differentially downregulated genes
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


# Label the top 10 differentially expressed genes
volcano$label <- ifelse(volcano$Gene %in% top_genes, volcano$Gene, NA)

# Plot volcano plot
ggplot(data=volcano, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=label)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("cornflowerblue", "grey", "firebrick2"),
                     labels=c("Downregulated", "Not significant", "Upregulated")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="grey25") +
  geom_hline(yintercept=-log10(0.05), col="grey25") +
  ylim(-5, 20) +
  xlim(-10, 10) +
  geom_text_repel(max.overlaps = 15, size = 3, force = 15, col = "black")

# Density plot
ggplot(volcano, aes(x=log2FoldChange, fill=diffexpressed)) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values = c("cornflowerblue", "grey", "firebrick2"), name="Differential Gene Expression") +
  theme_classic() +
  labs(title="Density plot of Log2 Fold Changes", x="Log2 Fold Change", y="Density", legend="He") +
  xlim(-8,6) +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black",
                                fill = NA, linewidth = 0.8)
  )

# GO ontology
# Upregulated pathways
upreg_genes <- as.vector(upreg$Gene)
upreg_GO <- enrichGO(gene = upreg_genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
as.data.frame(upreg_GO)
fit <- plot(barplot(upreg_GO, showCategory = 10))
fit
write.csv(as.data.frame(upreg_GO),file="C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/GO_analysis/upreg/hs4_upreg.csv")

# Downregulated pathways
downreg_genes <- as.vector(downreg$Gene)
downreg_GO <- enrichGO(gene = downreg_genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
as.data.frame(downreg_GO)
fit <- plot(barplot(downreg_GO, showCategory = 10))
fit
write.csv(as.data.frame(downreg_GO),file="C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/GO_analysis/downreg/hs4_downreg.csv")


# Top 20 upregulated significant genes (padj < 0.05)
upreg_20 <- upreg %>%
  arrange(desc(-padj))
upreg_20 <- upreg_20[1:20, ]
write.table(upreg_20["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/hs4/hs4_upreg_20_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)

# Top 50 upregulated significant genes (padj < 0.05)
upreg_50 <- upreg %>%
  arrange(desc(-padj))
upreg_50 <- upreg_50[1:50, ]
write.table(upreg_50["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/hs4/hs4_upreg_50_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)


# Top 10 upregulated significant genes (padj < 0.05)
upreg_10 <- upreg %>%
  arrange(desc(-padj))
upreg_10 <- upreg_10[1:10, ]
write.table(upreg_10["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/hs4/hs4_upreg_10_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)


# summary of the top 10 upregulated significant genes (padj < 0.05)
upreg_10_sum <- upreg_10 %>%
  dplyr::select(Gene, padj, log2FoldChange)
write.table(upreg_10_sum,
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/hs4/hs4_upreg_summary.txt",
            sep = ",",
            row.names = FALSE,
            col.names = FALSE
)

# Top 20 down-regulated significant genes (padj < 0.05)
downreg_20 <- downreg %>%
  arrange(desc(-padj))
downreg_20 <- downreg_20[1:20, ]
write.table(downreg_20["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/hs4/hs4_downreg_20_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)


# Top 50 down-regulated significant genes (padj < 0.05)
downreg_50 <- downreg %>%
  arrange(desc(-padj))
downreg_50 <- downreg_50[1:50, ]
write.table(downreg_50["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/hs4/hs4_downreg_50_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)


# Top 10 down-regulated significant genes (padj < 0.05)
downreg_10 <- downreg %>%
  arrange(desc(-padj))
downreg_10 <- downreg_10[1:10, ]
write.table(downreg_10["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/hs4/hs4_downreg_10_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)

# 10 down-regulated significant genes summaries (padj < 0.05)
downreg_10_sum <- downreg_10 %>%
  dplyr::select(Gene, padj, log2FoldChange)
write.table(downreg_10_sum,
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/hs4/hs4_downreg_summary.txt",
            sep = ",",
            row.names = FALSE,
            col.names = FALSE
)

# Upregulated, genes where padj < 0.05
upreg_05 <- upreg %>%
  arrange(desc(log2FoldChange))
write.table(upreg_05["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/hs4/hs4_upreg_05_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)


# Downregulated, genes where padj < 0.05
downreg_05 <- downreg %>%
  arrange(desc(-log2FoldChange))
write.table(downreg_05["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/hs4/hs4_downreg_05_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)


# Pvalues for the top significant genes across all healthy vs cancer DGEA comparisons
# Not including stage vs stages since they have specific gene versions
top_genes <- c("ENSG00000254772",
               "ENSG00000280614",
               "ENSG00000281181",
               "ENSG00000136305",
               "ENSG00000178922",
               "ENSG00000240801",
               "ENSG00000262160",
               "ENSG00000214736",
               "ENSG00000259948",
               "ENSG00000282034",
               "ENSG00000210196",
               "ENSG00000199753",
               "ENSG00000256393",
               "ENSG00000231131",
               "ENSG00000250982",
               "ENSG00000265452",
               "ENSG00000238249",
               "ENSG00000237638",
               "ENSG00000215006",
               "ENSG00000213793",
               "ENSG00000259040",
               "ENSG00000227063",
               "ENSG00000225953",
               "ENSG00000099977",
               "ENSG00000105583",
               "ENSG00000183054",
               "ENSG00000237039",
               "ENSG00000243176",
               "ENSG00000227077",
               "ENSG00000136149",
               "ENSG00000206557",
               "ENSG00000232368",
               "ENSG00000186940",
               "ENSG00000228526",
               "ENSG00000106031",
               "ENSG00000011052",
               "ENSG00000058673",
               "ENSG00000129965",
               "ENSG00000132207",
               "ENSG00000140488",
               "ENSG00000142684",
               "ENSG00000142694"
)

top_genes_data <- volcano %>%
  dplyr::select(Gene, padj, diffexpressed) %>%
  filter(Gene %in% top_genes) %>%
  arrange(factor(Gene, levels = top_genes)) # Order genes based on vector order

# Save data
write.table(top_genes_data,
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/hs4/top_dgea_genes.txt",
            sep = ",",
            row.names = FALSE,
            col.names = TRUE
)



# Pvalues for the top 5 up and down significant genes across all DGEA comparisons (including stage vs stage and paired)
top_genes <- c("ENSG00000254772",
               "ENSG00000280614",
               "ENSG00000281181",
               "ENSG00000136305",
               "ENSG00000178922",
               "ENSG00000210196",
               "ENSG00000199753",
               "ENSG00000256393",
               "ENSG00000231131",
               "ENSG00000250982",
               "ENSG00000240801",
               "ENSG00000213793",
               "ENSG00000099977",
               "ENSG00000237039",
               "ENSG00000227063",
               "ENSG00000243176",
               "ENSG00000227077",
               "ENSG00000136149",
               "ENSG00000011052",
               "ENSG00000058673",
               "ENSG00000105583",
               "ENSG00000129965",
               "ENSG00000287928",
               "ENSG00000288048",
               "ENSG00000234389",
               "ENSG00000273209",
               "ENSG00000176040",
               "ENSG00000213867",
               "ENSG00000236824",
               "ENSG00000130822",
               "ENSG00000159337",
               "ENSG00000214510",
               "ENSG00000282301",
               "ENSG00000164093",
               "ENSG00000181577",
               "ENSG00000107159",
               "ENSG00000205358",
               "ENSG00000125144",
               "ENSG00000019169",
               "ENSG00000160339",
               "ENSG00000198417",
               "ENSG00000175063",
               "ENSG00000164283",
               "ENSG00000147257",
               "ENSG00000204291",
               "ENSG00000117399",
               "ENSG00000164611",
               "ENSG00000101412",
               "ENSG00000089685",
               "ENSG00000145708"
)

top_genes_data <- volcano %>%
  dplyr::select(Gene, padj, diffexpressed) %>%
  filter(Gene %in% top_genes) %>%
  arrange(factor(Gene, levels = top_genes)) # Order genes based on vector order

# Save data
write.table(top_genes_data,
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/hs4/all_top5_dgea_genes.txt",
            sep = ",",
            row.names = FALSE,
            col.names = TRUE
)

