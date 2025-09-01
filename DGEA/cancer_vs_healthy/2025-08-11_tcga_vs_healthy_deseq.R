# All stages combined vs healthy DGEA

library("DESeq2")
library(ggplot2)
library(ggrepel)
library("magrittr")
library("vsn")
library("RColorBrewer")
library("pheatmap")
library(tidyverse)
library(reshape2)

# P values in ggplot figs
library(rstatix)
library(ggpubr) # for shapiro wilk test & box plots, pvalues

# GO ontology
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Clear environment
rm(list=ls())

# Load in sample raw data and sample condition
setwd("C:/Users/leahe/Documents/LIFE4137_independent_project/xena")
count_data <- read.csv(file = "removed_versions/healthy_vs_cancer/healthy_stages_rawcounts.csv", row.names = 1, check.names = FALSE)
condition <- read.csv(file = "clean_script/healthy_vs_tcga/healthy_tcga_raw_col_data.csv", header = TRUE, row.names = 1)
# Load in TCGA meta to make a SPINK13 gene counts figure to compare expression amongst sample type
cancer_meta <- read.table(file = "used_samples_metadata/TCGA_metadata.txt", sep = '\t', header = T)

# Pair row names
colnames(count_data) <- row.names(condition)
count_data <- round(count_data, 0) # Make values integers
# Make dds object for DESeq
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = condition,
                              design = ~ condition)

resultsNames(dds)
levels(dds$condition)

# Re-level so healthy is the control i.e healthy_vs_HCC -> HCC_vs_healthy
dds$condition %<>% relevel("Healthy")
levels(dds$condition)

dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name="condition_HCC_vs_Healthy") # y (HCC) vs x control (healthy)

# P-values and adj-pvalues
summary(res)
# Number of genes with padj < 0.01
sum(res$padj < 0.01, na.rm=TRUE) # 18079
# Number of genes with padj < 0.05
sum(res$padj < 0.05, na.rm=TRUE) # 20942


# Figure for SPINK13
cancer_meta <- cancer_meta %>%
  dplyr::select(sample, ajcc_pathologic_stage.diagnoses)

# Group subgroups in just one group e.g stage 4A = stage4
cancer_meta$stage <- cancer_meta$ajcc_pathologic_stage.diagnoses
cancer_meta$stage[cancer_meta$stage == "Stage IIIA" | cancer_meta$stage == "Stage IIIB" | cancer_meta$stage == "Stage IIIC"] <- "Stage III"
cancer_meta$stage[cancer_meta$stage == "Stage IVA" | cancer_meta$stage == "Stage IVB"] <- "Stage IV"
# Make sample names the rownames
rownames(cancer_meta) <- cancer_meta$sample
cancer_meta <- cancer_meta %>%
  dplyr::select(stage)

SPINK13 <- plotCounts(dds, gene="ENSG00000214510", intgroup="condition", 
                      returnData=TRUE)
# merge cancer meta with SPINK13 gene counts
SPINK13 <- merge(SPINK13, cancer_meta,
                    by = 'row.names', all = TRUE)
# Assign Healthy category to non-cancer samples
SPINK13$stage[is.na(SPINK13$stage)] <- "Healthy"
names(SPINK13)[names(SPINK13) == "Row.names"] <- "sample"

# Log counts used (outliers retained)
SPINK13 <- SPINK13 %>%
  select(count, condition, stage)

# Reorder following a precise order
SPINK13 <- SPINK13 %>%
  mutate(condition = fct_relevel(stage, 
                                 "Stage II", "Healthy", "Stage I", "Stage III", "Stage IV")
  )

plot <- ggplot(SPINK13, aes(x = condition, y = count, col=stage)) +
  geom_boxplot() +
  geom_point(alpha = 0.2, position = "jitter") +
  xlab("Condition") +
  ylab("Log gene counts") +
  theme_bw(base_size = 16) +
  theme(legend.box.margin=margin(,,,-13),
        legend.text = element_text(size=11), 
        legend.title=element_text(size=13)
  ) +
  ylim(0,12)
plot

my_comparisons <- list( c("Stage II", "Healthy"), c("Stage II", "Stage I"), c("Stage II", "Stage III"), c("Stage II", "Stage IV") )

plot + stat_compare_means(method = "wilcox", comparisons = my_comparisons, label.y = c(7.5, 8.5, 9.5, 10.5))



# MA plot 
plotMA(res, alpha = 0.05, ylim=c(-10, 10))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="limegreen", cex=2, lwd=2)
  text(baseMean, (log2FoldChange-0.5), topGene, pos=1, col="limegreen")
})

# Most significant DGE gene counts
clean_d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                      returnData=TRUE)

ggplot(clean_d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10() +
  ylab("Log count") +
  xlab("Condition") +
  theme_bw()

# Most significant DGE
which.min(res$padj) # 279 - ENSG00000011052 gene
healthy_gene_count <- clean_d %>%
  filter(condition=="Healthy")
mean(healthy_gene_count[["count"]]) # 7.481639
median(healthy_gene_count[["count"]])  # 7.685022

HCC_gene_count <- clean_d %>%
  filter(condition=="HCC")
mean(HCC_gene_count[["count"]]) # 0.5
median(HCC_gene_count[["count"]]) # 0.5


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
  theme_bw(base_size = 16)

# Volcano plot
# Subset data
volcano <- data.frame(res@rownames, res$padj, res$log2FoldChange)
# Rename columns
colnames(volcano) <- c("Gene", "padj", "log2FoldChange")
volcano <- na.omit(volcano)
volcano$diffexpressed <- "NO"

# Up-regulated genes if pvalue < 0.05 and 0.6 < log2Foldchange
volcano$diffexpressed[volcano$log2FoldChange > 0.6 & volcano$padj < 0.05] <- "UP"
# Down-regulated genes if pvalue < 0.05 and log2Foldchange < 0.6
volcano$diffexpressed[volcano$log2FoldChange < -0.6 & volcano$padj < 0.05] <- "DOWN"

# Label the top 5 DGEs
top_genes <- c()
# Get the top 5 most significant up-regulated genes
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

# Get the top 5 most significant down-regulated genes
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
replacement_values

# Label the top 10 DGEs
volcano$label <- ifelse(volcano$Gene %in% top_genes, volcano$Gene, NA)

# Convert geneID label to gene symbol
volcano$label[volcano$label=="ENSG00000210196"]<-"MT.TP"
volcano$label[volcano$label=="ENSG00000199753"]<-"SNORD104"
volcano$label[volcano$label=="ENSG00000256393"]<-"RPL41P5"
volcano$label[volcano$label=="ENSG00000250982"]<-"GAPDHP35"
volcano$label[volcano$label=="ENSG00000231131"]<-"LNCAROD"
volcano$label[volcano$label=="ENSG00000011052"]<-"NME1.NME2"
volcano$label[volcano$label=="ENSG00000058673"]<-"ZC3H11A"
volcano$label[volcano$label=="ENSG00000099977"]<-"DDT"
volcano$label[volcano$label=="ENSG00000105583"]<-"WDR83OS"
volcano$label[volcano$label=="ENSG00000129965"]<-"INS.IGF2"


# Plot volcano plot
ggplot(data=volcano, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=label)) +
  geom_point(alpha=0.5) + 
  theme_minimal(base_size = 16) +
  scale_color_manual(values=c("cornflowerblue", "grey", "firebrick2"),
                     labels=c("Downregulated", "No", "Upregulated")) +
  geom_vline(xintercept=c(-0.6, 0.6), col = "grey25") +
  geom_hline(yintercept=-log10(0.05), col = "grey25") +
  ylim(-5, 350) +
  xlim(-10, 10) +
  geom_text_repel(max.overlaps = Inf, size = 5, force = 10, col = "black") +
  theme( legend.box.margin=margin(,,,-13),
         legend.text = element_text(size=11), 
         legend.title=element_text(size=13)
  )

# Density plot
ggplot(volcano, aes(x=log2FoldChange, fill=diffexpressed)) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values = c("cornflowerblue", "grey", "firebrick2"), name = "Differential Gene Expression") +
  theme_classic() +
  labs(title = "Density plot of Log2 Fold Changes", x = "Log2 Fold Change", y = "Density", legend = "He") +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black",
                                fill = NA, linewidth = 0.8)
  ) +
  xlim(-8,6)


# GO ontology
# Upregulated pathways
upreg_genes <- as.vector(upreg$Gene)
upreg_GO <- enrichGO(gene = upreg_genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
as.data.frame(upreg_GO)
fit <- plot(barplot(upreg_GO, showCategory = 10))
fit
write.csv(as.data.frame(upreg_GO),file="C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/GO_analysis/upreg/healthy_vs_cancer_upreg.csv")


# Downregulated pathways
downreg_genes <- as.vector(downreg$Gene)
downreg_GO <- enrichGO(gene = downreg_genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
as.data.frame(downreg_GO)
fit <- plot(barplot(downreg_GO, showCategory = 10))
fit
write.csv(as.data.frame(downreg_GO),file="C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/GO_analysis/downreg/healthy_vs_cancer_downreg.csv")


# Top 20 upregulated significant genes (padj < 0.05)
upreg_20 <- upreg %>%
  arrange(desc(-padj))
upreg_20 <- upreg_20[1:20, ]
write.table(upreg_20["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/healthy_vs_cancer/h_vs_tcga_upreg_20_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)

# Top 50 upregulated significant genes (padj < 0.05)
upreg_50 <- upreg %>%
  arrange(desc(-padj))
upreg_50 <- upreg_50[1:50, ]
write.table(upreg_50["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/healthy_vs_cancer/h_vs_tcga_upreg_50_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)

# Top 10 upregulated significant genes (padj < 0.05)
upreg_10 <- upreg %>%
  arrange(desc(-padj))
upreg_10 <- upreg_10[1:10, ]
write.table(upreg_10["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/healthy_vs_cancer/h_vs_tcga_upreg_10_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)



# Top 20 down-regulated significant genes (padj < 0.05)
downreg_20 <- downreg %>%
  arrange(desc(-padj))
downreg_20 <- downreg_20[1:20, ]
write.table(downreg_20["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/healthy_vs_cancer/h_vs_tcga_downreg_20_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)

# Top 50 down-regulated significant genes (padj < 0.05)
downreg_50 <- downreg %>%
  arrange(desc(-padj))
downreg_50 <- downreg_50[1:50, ]
write.table(downreg_50["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/healthy_vs_cancer/h_vs_tcga_downreg_50_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)

# Top 10 down-regulated significant genes (padj < 0.05)
downreg_10 <- downreg %>%
  arrange(desc(-padj))
downreg_10 <- downreg_10[1:10, ]
write.table(downreg_10["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/healthy_vs_cancer/h_vs_tcga_downreg_10_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)

# Summary of the Top 10 down-regulated significant genes (padj < 0.05)
downreg_10_sum <- downreg_10 %>%
  dplyr::select(Gene, padj, log2FoldChange)
write.table(downreg_10_sum,
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/healthy_vs_cancer/h_vs_downreg_10_summary.txt",
            sep = ",",
            row.names = FALSE,
            col.names = FALSE
)

# Summary of the Top 10 down-regulated significant genes (padj < 0.05)
upreg_10_sum <- upreg_10 %>%
  dplyr::select(Gene, padj, log2FoldChange)
write.table(upreg_10_sum,
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/healthy_vs_cancer/h_vs_tcga_upreg_summary.txt",
            sep = ",",
            row.names = FALSE,
            col.names = FALSE
)


# Upregulated, genes where padj < 0.05
upreg_05 <- upreg %>%
  arrange(desc(log2FoldChange))
write.table(upreg_05["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/healthy_vs_cancer/h_vs_tcga_upreg_05_genes.txt",
            sep = "",
            row.names = FALSE,
            col.names = FALSE
)

# Downregulated, genes where padj < 0.05
downreg_05 <- downreg %>%
  arrange(desc(-log2FoldChange))
write.table(downreg_05["Gene"],
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/healthy_vs_cancer/h_vs_tcga_downreg_05_genes.txt",
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
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/healthy_vs_cancer/top_dgea_genes.txt",
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
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/removed_versions/healthy_vs_cancer/all_top5_dgea_genes.txt",
            sep = ",",
            row.names = FALSE,
            col.names = TRUE
)
