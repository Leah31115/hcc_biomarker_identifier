# Clear environment
rm(list=ls())

library(tidyverse)
library(survival)
library(survminer)

# Open TPM TCGA paired survival data
setwd("C:/Users/leahe/Documents/LIFE4137_independent_project/xena")
TCGA_samples <- read.table(file="cancer_samples/normalised_tcga_star_counts/TCGA-LIHC.star_tpm.tsv.gz", sep='\t', header=T, check.names = FALSE)
tcga_survival <- read.table(file="all_survival_LIHC_data.txt", sep='\t', header=T)
tcga_meta <- read.table(file="used_samples_metadata/TCGA_metadata.txt", sep='', header=T)

# Filter TPM TCGA for used samples
cancer_sample_names = tcga_meta[["samples"]]
cancer_sample_names

# Filter for cancer_samples
cancer_samples <- TCGA_samples %>%
  select("Ensembl_ID", all_of(cancer_sample_names))
# These genes are protein coding and upregulated from the paired TCGA analysis
downreg_genes <- c("ENSG00000263761.3",
                   "ENSG00000163217.2",
                   "ENSG00000104938.18",
                   "ENSG00000165682.14",
                   "ENSG00000205358.4")


# Filter rows for paired tcga dgea downregulated genes
filtered_counts <- cancer_samples[cancer_samples$Ensembl_ID %in% downreg_genes, ]

# Format paired counts dataframe
# Make ensemblID column the row names
rownames(filtered_counts) <- filtered_counts$Ensembl_ID
# Remove ensemblID column
filtered_counts <- filtered_counts[,!(names(filtered_counts) %in% "Ensembl_ID")]

# Make the columns (sample names) the rows
filtered_counts <- as.data.frame(t(filtered_counts))

# Add "A" to the end of all samples in the all samples
tcga_survival$sample <- paste0(tcga_survival$sample, "A")

# Filter survival data for used samples
filtered_tcga_survival <- tcga_survival[tcga_survival$sample %in% cancer_sample_names, ]

redacted_samples <- filtered_tcga_survival %>%
  filter(filtered_tcga_survival$Redaction == "Redacted")

# Filter out redacted samples
filtered_tcga_survival <- filtered_tcga_survival %>%
  filter(! filtered_tcga_survival$Redaction == "Redacted")

# Make sample column the row names
rownames(filtered_tcga_survival) <- filtered_tcga_survival$sample

# Remove unneeded columns
filtered_tcga_survival <- filtered_tcga_survival[, !names(filtered_tcga_survival) %in% c("sample", "Redaction", "DFI!", "DFI.time")]

# Merge survival data with gene counts 
survival_df <- merge(filtered_counts, filtered_tcga_survival, by = "row.names", 
                     all.x = TRUE) # all.x will give NA for missing values

# Rename gene ID to gene symbol name
names(survival_df)[names(survival_df) == "ENSG00000263761.3"] <- "GDF2"
names(survival_df)[names(survival_df) == "ENSG00000163217.2"] <- "BMP10"
names(survival_df)[names(survival_df) == "ENSG00000104938.18"] <- "CLEC4M"
names(survival_df)[names(survival_df) == "ENSG00000165682.14"] <- "CLEC1B"
names(survival_df)[names(survival_df) == "ENSG00000205358.4"] <- "MT1H"


# Apply function to new transformed gene count column
transform_tpm_counts <- function(count){
  transformed_count <- (2^count) - 0.001
  return(transformed_count)
}

survival_df["transformed_GDF2"] <- lapply(survival_df["GDF2"], transform_tpm_counts)
survival_df["transformed_BMP10"] <- lapply(survival_df["BMP10"], transform_tpm_counts)
survival_df["transformed_CLEC4M"] <- lapply(survival_df["CLEC4M"], transform_tpm_counts)
survival_df["transformed_CLEC1B"] <- lapply(survival_df["CLEC1B"], transform_tpm_counts)
survival_df["transformed_MT1H"] <- lapply(survival_df["MT1H"], transform_tpm_counts)


# 1 day = 0.0328767 days 
days_to_months_converter <- function(no_days){
  no_months <- round(no_days*0.0328767, 1)
  return(no_months)
}

survival_df["OS_months"] <- lapply(survival_df["OS.time"], days_to_months_converter)
survival_df["DSS_months"] <- lapply(survival_df["DSS.time"], days_to_months_converter)
survival_df["PFI_months"] <- lapply(survival_df["PFI.time"], days_to_months_converter)

# GDF2
# Make survival plot
survival_df$GDF2_median <- NA
survival_df$GDF2_median[survival_df$transformed_GDF2 <= quantile(survival_df$transformed_GDF2, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$GDF2_median[survival_df$transformed_GDF2 > quantile(survival_df$transformed_GDF2, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$GDF2_25p <- NA
survival_df$GDF2_25p[survival_df$transformed_GDF2 <= quantile(survival_df$transformed_GDF2, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$GDF2_25p[survival_df$transformed_GDF2 > quantile(survival_df$transformed_GDF2, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$GDF2_75p <- NA
survival_df$GDF2_75p[survival_df$transformed_GDF2 <= quantile(survival_df$transformed_GDF2, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$GDF2_75p[survival_df$transformed_GDF2 > quantile(survival_df$transformed_GDF2, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$GDF2_66p <- NA
survival_df$GDF2_66p[survival_df$transformed_GDF2 <= quantile(survival_df$transformed_GDF2, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$GDF2_66p[survival_df$transformed_GDF2 > quantile(survival_df$transformed_GDF2, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$GDF2_33p <- NA
survival_df$GDF2_33p[survival_df$transformed_GDF2 <= quantile(survival_df$transformed_GDF2, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$GDF2_33p[survival_df$transformed_GDF2 > quantile(survival_df$transformed_GDF2, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$GDF2_33_66 <- NA
survival_df$GDF2_33_66[survival_df$transformed_GDF2 <= quantile(survival_df$transformed_GDF2, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$GDF2_33_66[survival_df$transformed_GDF2 >= quantile(survival_df$transformed_GDF2, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - GDF2
for(strata in c("GDF2_median", "GDF2_25p", "GDF2_75p", "GDF2_66p", "GDF2_33p", "GDF2_33_66")){
  
  H_gTxt <- paste("HIGH GDF2", sep="")
  L_gTxt <- paste("LOW GDF2", sep="")
  suffix <- paste(datasetName, strata, sep="_")
  
  sc = "OS"
  tmp.sub <- subset(survival_df, !is.na(OS) & !is.na(OS_months))
  tmp.sub$usedSubset <- as.factor(tmp.sub[[strata]])
  tmp.sub <- subset(tmp.sub,!is.na(usedSubset))
  
  fit <- survfit(Surv(OS_months, OS) ~ usedSubset, data= tmp.sub)
  p <-  ggsurvplot(fit, palette = c("#ED0000FF", "black"), title=paste(strata,sep=" - "), xlab="Months", ylab="Overall Survival (OS)", 
                   break.time.by = 30,
                   xlim = c(0, 125),
                   font.x = c(18, "black"),
                   font.y = c(18, "black"),
                   font.tickslab = c(18, "plain", "black"),
                   legend = c(0.25, 0.15),
                   legend.title = "",
                   legend.labs =c(paste0 (H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                  paste0 (L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                   font.legend = c(18, "black"),
                   fun = function(y) y*100,
                   pval = TRUE, pval.coord=c(70, 80), pval.size = 5,
                   size = 1, censor.size = 2)
  
  tiff(filename = paste("survival_plots/downreg/GDF2/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
  
  sc = "DSS"
  tmp.sub <- subset(survival_df, !is.na(DSS) & !is.na(DSS_months))
  tmp.sub$usedSubset <- as.factor(tmp.sub[[strata]])
  tmp.sub <- subset(tmp.sub,!is.na(usedSubset))
  
  fit <- survfit(Surv(DSS_months, DSS) ~ usedSubset, data = tmp.sub)
  p <- ggsurvplot(fit, palette = c("#ED0000FF", "black"),title=paste(strata,sep=" - "), xlab="Months", ylab="Disease Specific Survival (DSS)", 
                  break.time.by = 30, 
                  xlim = c(0, 125),
                  font.x = c(18, "black"),
                  font.y = c(18, "black"),
                  font.tickslab = c(18, "plain", "black"),
                  legend = c(0.25, 0.15),
                  legend.title = "",
                  legend.labs =c(paste0 (H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                 paste0 (L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                  font.legend = c(18, "black"),
                  fun = function(y) y*100,
                  pval = TRUE, pval.coord = c(70, 80), pval.size = 5,
                  size = 1, censor.size = 2)
  
  tiff(filename = paste ("survival_plots/downreg/GDF2/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12,compression = c("none"),bg = "white")
  print(p)
  dev.off()
  
  
  sc = "PFI"
  tmp.sub <- subset(survival_df,!is.na(PFI) & !is.na(PFI_months))
  tmp.sub$usedSubset <- as.factor(tmp.sub[[strata]])
  tmp.sub <- subset(tmp.sub,!is.na(usedSubset))
  
  fit<- survfit(Surv(PFI_months, PFI) ~ usedSubset, data= tmp.sub)
  p <-  ggsurvplot(fit, palette = c("#ED0000FF","black"),title=paste(strata,sep=" - "), xlab="Months", ylab="Progression Free Interval (PFI)",
                   break.time.by = 30, 
                   xlim = c(0, 125),
                   font.x = c(18, "black"),
                   font.y = c(18, "black"),
                   font.tickslab = c(18, "plain", "black"),
                   legend = c(0.25, 0.15),
                   legend.title = "",
                   legend.labs = c(paste0(H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                   paste0(L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                   font.legend = c(18, "black"),
                   fun = function(y) y*100,
                   pval = TRUE, pval.coord = c(70, 80), pval.size = 5,
                   size = 1, censor.size = 2)
  
  tiff(filename =  paste("survival_plots/downreg/GDF2/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}