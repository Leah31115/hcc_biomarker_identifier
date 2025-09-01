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
# IF I WANT TO INVESTIGATE MORE GENES, ADD THEM TO THIS VECTOR!!!!
# These genes are protein coding and upregulated from the paired TCGA analysis
downreg_genes <- c("ENSG00000263761.3",
                   "ENSG00000163217.2",
                   "ENSG00000104938.18",
                   "ENSG00000165682.14",
                   "ENSG00000205358.4",
                   "ENSG00000019169.10",
                   "ENSG00000160339.16",
                   "ENSG00000145824.13",
                   "ENSG00000120057.5",
                   "ENSG00000136011.15")


# Filter for desired genes HERE!!!
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
names(survival_df)[names(survival_df) == "ENSG00000019169.10"] <- "MARCO"
names(survival_df)[names(survival_df) == "ENSG00000160339.16"] <- "FCN2"
names(survival_df)[names(survival_df) == "ENSG00000145824.13"] <- "CXCL14"
names(survival_df)[names(survival_df) == "ENSG00000120057.5"] <- "SFRP5"
names(survival_df)[names(survival_df) == "ENSG00000136011.15"] <- "STAB2"


# Apply function to new transformed gene count column
transform_tpm_counts <- function(count){
  transformed_count <- (2^count) + 0.001
  return(transformed_count)
}

survival_df["transformed_GDF2"] <- lapply(survival_df["GDF2"], transform_tpm_counts)
survival_df["transformed_BMP10"] <- lapply(survival_df["BMP10"], transform_tpm_counts)
survival_df["transformed_CLEC4M"] <- lapply(survival_df["CLEC4M"], transform_tpm_counts)
survival_df["transformed_CLEC1B"] <- lapply(survival_df["CLEC1B"], transform_tpm_counts)
survival_df["transformed_MT1H"] <- lapply(survival_df["MT1H"], transform_tpm_counts)
survival_df["transformed_MARCO"] <- lapply(survival_df["MARCO"], transform_tpm_counts)
survival_df["transformed_FCN2"] <- lapply(survival_df["FCN2"], transform_tpm_counts)
survival_df["transformed_CXCL14"] <- lapply(survival_df["CXCL14"], transform_tpm_counts)
survival_df["transformed_SFRP5"] <- lapply(survival_df["SFRP5"], transform_tpm_counts)
survival_df["transformed_STAB2"] <- lapply(survival_df["STAB2"], transform_tpm_counts)



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

# BMP10
# Make survival plot
survival_df$BMP10_median <- NA
survival_df$BMP10_median[survival_df$transformed_BMP10 <= quantile(survival_df$transformed_BMP10, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$BMP10_median[survival_df$transformed_BMP10 > quantile(survival_df$transformed_BMP10, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$BMP10_25p <- NA
survival_df$BMP10_25p[survival_df$transformed_BMP10 <= quantile(survival_df$transformed_BMP10, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$BMP10_25p[survival_df$transformed_BMP10 > quantile(survival_df$transformed_BMP10, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$BMP10_75p <- NA
survival_df$BMP10_75p[survival_df$transformed_BMP10 <= quantile(survival_df$transformed_BMP10, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$BMP10_75p[survival_df$transformed_BMP10 > quantile(survival_df$transformed_BMP10, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$BMP10_66p <- NA
survival_df$BMP10_66p[survival_df$transformed_BMP10 <= quantile(survival_df$transformed_BMP10, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$BMP10_66p[survival_df$transformed_BMP10 > quantile(survival_df$transformed_BMP10, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$BMP10_33p <- NA
survival_df$BMP10_33p[survival_df$transformed_BMP10 <= quantile(survival_df$transformed_BMP10, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$BMP10_33p[survival_df$transformed_BMP10 > quantile(survival_df$transformed_BMP10, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$BMP10_33_66 <- NA
survival_df$BMP10_33_66[survival_df$transformed_BMP10 <= quantile(survival_df$transformed_BMP10, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$BMP10_33_66[survival_df$transformed_BMP10 >= quantile(survival_df$transformed_BMP10, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - BMP10
for(strata in c("BMP10_median", "BMP10_25p", "BMP10_75p", "BMP10_66p", "BMP10_33p", "BMP10_33_66")){
  
  H_gTxt <- paste("HIGH BMP10", sep="")
  L_gTxt <- paste("LOW BMP10", sep="")
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
  
  tiff(filename = paste("survival_plots/downreg/BMP10/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/downreg/BMP10/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/downreg/BMP10/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}



# CLEC4M
# Make survival plot
survival_df$CLEC4M_median <- NA
survival_df$CLEC4M_median[survival_df$transformed_CLEC4M <= quantile(survival_df$transformed_CLEC4M, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$CLEC4M_median[survival_df$transformed_CLEC4M > quantile(survival_df$transformed_CLEC4M, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$CLEC4M_25p <- NA
survival_df$CLEC4M_25p[survival_df$transformed_CLEC4M <= quantile(survival_df$transformed_CLEC4M, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$CLEC4M_25p[survival_df$transformed_CLEC4M > quantile(survival_df$transformed_CLEC4M, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$CLEC4M_75p <- NA
survival_df$CLEC4M_75p[survival_df$transformed_CLEC4M <= quantile(survival_df$transformed_CLEC4M, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$CLEC4M_75p[survival_df$transformed_CLEC4M > quantile(survival_df$transformed_CLEC4M, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$CLEC4M_66p <- NA
survival_df$CLEC4M_66p[survival_df$transformed_CLEC4M <= quantile(survival_df$transformed_CLEC4M, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$CLEC4M_66p[survival_df$transformed_CLEC4M > quantile(survival_df$transformed_CLEC4M, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$CLEC4M_33p <- NA
survival_df$CLEC4M_33p[survival_df$transformed_CLEC4M <= quantile(survival_df$transformed_CLEC4M, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CLEC4M_33p[survival_df$transformed_CLEC4M > quantile(survival_df$transformed_CLEC4M, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$CLEC4M_33_66 <- NA
survival_df$CLEC4M_33_66[survival_df$transformed_CLEC4M <= quantile(survival_df$transformed_CLEC4M, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CLEC4M_33_66[survival_df$transformed_CLEC4M >= quantile(survival_df$transformed_CLEC4M, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - CLEC4M
for(strata in c("CLEC4M_median", "CLEC4M_25p", "CLEC4M_75p", "CLEC4M_66p", "CLEC4M_33p", "CLEC4M_33_66")){
  
  H_gTxt <- paste("HIGH CLEC4M", sep="")
  L_gTxt <- paste("LOW CLEC4M", sep="")
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
  
  tiff(filename = paste("survival_plots/downreg/CLEC4M/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/downreg/CLEC4M/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/downreg/CLEC4M/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

# CLEC1B
# Make survival plot
survival_df$CLEC1B_median <- NA
survival_df$CLEC1B_median[survival_df$transformed_CLEC1B <= quantile(survival_df$transformed_CLEC1B, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$CLEC1B_median[survival_df$transformed_CLEC1B > quantile(survival_df$transformed_CLEC1B, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$CLEC1B_25p <- NA
survival_df$CLEC1B_25p[survival_df$transformed_CLEC1B <= quantile(survival_df$transformed_CLEC1B, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$CLEC1B_25p[survival_df$transformed_CLEC1B > quantile(survival_df$transformed_CLEC1B, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$CLEC1B_75p <- NA
survival_df$CLEC1B_75p[survival_df$transformed_CLEC1B <= quantile(survival_df$transformed_CLEC1B, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$CLEC1B_75p[survival_df$transformed_CLEC1B > quantile(survival_df$transformed_CLEC1B, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$CLEC1B_66p <- NA
survival_df$CLEC1B_66p[survival_df$transformed_CLEC1B <= quantile(survival_df$transformed_CLEC1B, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$CLEC1B_66p[survival_df$transformed_CLEC1B > quantile(survival_df$transformed_CLEC1B, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$CLEC1B_33p <- NA
survival_df$CLEC1B_33p[survival_df$transformed_CLEC1B <= quantile(survival_df$transformed_CLEC1B, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CLEC1B_33p[survival_df$transformed_CLEC1B > quantile(survival_df$transformed_CLEC1B, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$CLEC1B_33_66 <- NA
survival_df$CLEC1B_33_66[survival_df$transformed_CLEC1B <= quantile(survival_df$transformed_CLEC1B, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CLEC1B_33_66[survival_df$transformed_CLEC1B >= quantile(survival_df$transformed_CLEC1B, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - CLEC1B
for(strata in c("CLEC1B_median", "CLEC1B_25p", "CLEC1B_75p", "CLEC1B_66p", "CLEC1B_33p", "CLEC1B_33_66")){
  
  H_gTxt <- paste("HIGH CLEC1B", sep="")
  L_gTxt <- paste("LOW CLEC1B", sep="")
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
  
  tiff(filename = paste("survival_plots/downreg/CLEC1B/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/downreg/CLEC1B/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/downreg/CLEC1B/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

# MT1H
# Make survival plot
survival_df$MT1H_median <- NA
survival_df$MT1H_median[survival_df$transformed_MT1H <= quantile(survival_df$transformed_MT1H, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$MT1H_median[survival_df$transformed_MT1H > quantile(survival_df$transformed_MT1H, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$MT1H_25p <- NA
survival_df$MT1H_25p[survival_df$transformed_MT1H <= quantile(survival_df$transformed_MT1H, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$MT1H_25p[survival_df$transformed_MT1H > quantile(survival_df$transformed_MT1H, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$MT1H_75p <- NA
survival_df$MT1H_75p[survival_df$transformed_MT1H <= quantile(survival_df$transformed_MT1H, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$MT1H_75p[survival_df$transformed_MT1H > quantile(survival_df$transformed_MT1H, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$MT1H_66p <- NA
survival_df$MT1H_66p[survival_df$transformed_MT1H <= quantile(survival_df$transformed_MT1H, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$MT1H_66p[survival_df$transformed_MT1H > quantile(survival_df$transformed_MT1H, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$MT1H_33p <- NA
survival_df$MT1H_33p[survival_df$transformed_MT1H <= quantile(survival_df$transformed_MT1H, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$MT1H_33p[survival_df$transformed_MT1H > quantile(survival_df$transformed_MT1H, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$MT1H_33_66 <- NA
survival_df$MT1H_33_66[survival_df$transformed_MT1H <= quantile(survival_df$transformed_MT1H, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$MT1H_33_66[survival_df$transformed_MT1H >= quantile(survival_df$transformed_MT1H, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - MT1H
for(strata in c("MT1H_median", "MT1H_25p", "MT1H_75p", "MT1H_66p", "MT1H_33p", "MT1H_33_66")){
  
  H_gTxt <- paste("HIGH MT1H", sep="")
  L_gTxt <- paste("LOW MT1H", sep="")
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
  
  tiff(filename = paste("survival_plots/downreg/MT1H/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/downreg/MT1H/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/downreg/MT1H/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}


# MARCO
# Make survival plot
survival_df$MARCO_median <- NA
survival_df$MARCO_median[survival_df$transformed_MARCO <= quantile(survival_df$transformed_MARCO, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$MARCO_median[survival_df$transformed_MARCO > quantile(survival_df$transformed_MARCO, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$MARCO_25p <- NA
survival_df$MARCO_25p[survival_df$transformed_MARCO <= quantile(survival_df$transformed_MARCO, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$MARCO_25p[survival_df$transformed_MARCO > quantile(survival_df$transformed_MARCO, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$MARCO_75p <- NA
survival_df$MARCO_75p[survival_df$transformed_MARCO <= quantile(survival_df$transformed_MARCO, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$MARCO_75p[survival_df$transformed_MARCO > quantile(survival_df$transformed_MARCO, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$MARCO_66p <- NA
survival_df$MARCO_66p[survival_df$transformed_MARCO <= quantile(survival_df$transformed_MARCO, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$MARCO_66p[survival_df$transformed_MARCO > quantile(survival_df$transformed_MARCO, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$MARCO_33p <- NA
survival_df$MARCO_33p[survival_df$transformed_MARCO <= quantile(survival_df$transformed_MARCO, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$MARCO_33p[survival_df$transformed_MARCO > quantile(survival_df$transformed_MARCO, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$MARCO_33_66 <- NA
survival_df$MARCO_33_66[survival_df$transformed_MARCO <= quantile(survival_df$transformed_MARCO, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$MARCO_33_66[survival_df$transformed_MARCO >= quantile(survival_df$transformed_MARCO, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - MARCO
for(strata in c("MARCO_median", "MARCO_25p", "MARCO_75p", "MARCO_66p", "MARCO_33p", "MARCO_33_66")){
  
  H_gTxt <- paste("HIGH MARCO", sep="")
  L_gTxt <- paste("LOW MARCO", sep="")
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
  
  tiff(filename = paste("survival_plots/downreg/MARCO/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/downreg/MARCO/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/downreg/MARCO/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}


# FCN2
# Make survival plot
survival_df$FCN2_median <- NA
survival_df$FCN2_median[survival_df$transformed_FCN2 <= quantile(survival_df$transformed_FCN2, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$FCN2_median[survival_df$transformed_FCN2 > quantile(survival_df$transformed_FCN2, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$FCN2_25p <- NA
survival_df$FCN2_25p[survival_df$transformed_FCN2 <= quantile(survival_df$transformed_FCN2, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$FCN2_25p[survival_df$transformed_FCN2 > quantile(survival_df$transformed_FCN2, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$FCN2_75p <- NA
survival_df$FCN2_75p[survival_df$transformed_FCN2 <= quantile(survival_df$transformed_FCN2, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$FCN2_75p[survival_df$transformed_FCN2 > quantile(survival_df$transformed_FCN2, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$FCN2_66p <- NA
survival_df$FCN2_66p[survival_df$transformed_FCN2 <= quantile(survival_df$transformed_FCN2, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$FCN2_66p[survival_df$transformed_FCN2 > quantile(survival_df$transformed_FCN2, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$FCN2_33p <- NA
survival_df$FCN2_33p[survival_df$transformed_FCN2 <= quantile(survival_df$transformed_FCN2, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$FCN2_33p[survival_df$transformed_FCN2 > quantile(survival_df$transformed_FCN2, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$FCN2_33_66 <- NA
survival_df$FCN2_33_66[survival_df$transformed_FCN2 <= quantile(survival_df$transformed_FCN2, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$FCN2_33_66[survival_df$transformed_FCN2 >= quantile(survival_df$transformed_FCN2, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - FCN2
for(strata in c("FCN2_median", "FCN2_25p", "FCN2_75p", "FCN2_66p", "FCN2_33p", "FCN2_33_66")){
  
  H_gTxt <- paste("HIGH FCN2", sep="")
  L_gTxt <- paste("LOW FCN2", sep="")
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
  
  tiff(filename = paste("survival_plots/downreg/FCN2/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/downreg/FCN2/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/downreg/FCN2/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}


# CXCL14
# Make survival plot
survival_df$CXCL14_median <- NA
survival_df$CXCL14_median[survival_df$transformed_CXCL14 <= quantile(survival_df$transformed_CXCL14, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$CXCL14_median[survival_df$transformed_CXCL14 > quantile(survival_df$transformed_CXCL14, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$CXCL14_25p <- NA
survival_df$CXCL14_25p[survival_df$transformed_CXCL14 <= quantile(survival_df$transformed_CXCL14, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$CXCL14_25p[survival_df$transformed_CXCL14 > quantile(survival_df$transformed_CXCL14, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$CXCL14_75p <- NA
survival_df$CXCL14_75p[survival_df$transformed_CXCL14 <= quantile(survival_df$transformed_CXCL14, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$CXCL14_75p[survival_df$transformed_CXCL14 > quantile(survival_df$transformed_CXCL14, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$CXCL14_66p <- NA
survival_df$CXCL14_66p[survival_df$transformed_CXCL14 <= quantile(survival_df$transformed_CXCL14, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$CXCL14_66p[survival_df$transformed_CXCL14 > quantile(survival_df$transformed_CXCL14, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$CXCL14_33p <- NA
survival_df$CXCL14_33p[survival_df$transformed_CXCL14 <= quantile(survival_df$transformed_CXCL14, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CXCL14_33p[survival_df$transformed_CXCL14 > quantile(survival_df$transformed_CXCL14, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$CXCL14_33_66 <- NA
survival_df$CXCL14_33_66[survival_df$transformed_CXCL14 <= quantile(survival_df$transformed_CXCL14, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CXCL14_33_66[survival_df$transformed_CXCL14 >= quantile(survival_df$transformed_CXCL14, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - CXCL14
for(strata in c("CXCL14_median", "CXCL14_25p", "CXCL14_75p", "CXCL14_66p", "CXCL14_33p", "CXCL14_33_66")){
  
  H_gTxt <- paste("HIGH CXCL14", sep="")
  L_gTxt <- paste("LOW CXCL14", sep="")
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
  
  tiff(filename = paste("survival_plots/downreg/CXCL14/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/downreg/CXCL14/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/downreg/CXCL14/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}



# SFRP5
# Make survival plot
survival_df$SFRP5_median <- NA
survival_df$SFRP5_median[survival_df$transformed_SFRP5 <= quantile(survival_df$transformed_SFRP5, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$SFRP5_median[survival_df$transformed_SFRP5 > quantile(survival_df$transformed_SFRP5, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$SFRP5_25p <- NA
survival_df$SFRP5_25p[survival_df$transformed_SFRP5 <= quantile(survival_df$transformed_SFRP5, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$SFRP5_25p[survival_df$transformed_SFRP5 > quantile(survival_df$transformed_SFRP5, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$SFRP5_75p <- NA
survival_df$SFRP5_75p[survival_df$transformed_SFRP5 <= quantile(survival_df$transformed_SFRP5, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$SFRP5_75p[survival_df$transformed_SFRP5 > quantile(survival_df$transformed_SFRP5, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$SFRP5_66p <- NA
survival_df$SFRP5_66p[survival_df$transformed_SFRP5 <= quantile(survival_df$transformed_SFRP5, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$SFRP5_66p[survival_df$transformed_SFRP5 > quantile(survival_df$transformed_SFRP5, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$SFRP5_33p <- NA
survival_df$SFRP5_33p[survival_df$transformed_SFRP5 <= quantile(survival_df$transformed_SFRP5, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$SFRP5_33p[survival_df$transformed_SFRP5 > quantile(survival_df$transformed_SFRP5, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$SFRP5_33_66 <- NA
survival_df$SFRP5_33_66[survival_df$transformed_SFRP5 <= quantile(survival_df$transformed_SFRP5, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$SFRP5_33_66[survival_df$transformed_SFRP5 >= quantile(survival_df$transformed_SFRP5, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - SFRP5
for(strata in c("SFRP5_median", "SFRP5_25p", "SFRP5_75p", "SFRP5_66p", "SFRP5_33p", "SFRP5_33_66")){
  
  H_gTxt <- paste("HIGH SFRP5", sep="")
  L_gTxt <- paste("LOW SFRP5", sep="")
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
  
  tiff(filename = paste("survival_plots/downreg/SFRP5/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/downreg/SFRP5/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/downreg/SFRP5/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}



# STAB2
# Make survival plot
survival_df$STAB2_median <- NA
survival_df$STAB2_median[survival_df$transformed_STAB2 <= quantile(survival_df$transformed_STAB2, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$STAB2_median[survival_df$transformed_STAB2 > quantile(survival_df$transformed_STAB2, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$STAB2_25p <- NA
survival_df$STAB2_25p[survival_df$transformed_STAB2 <= quantile(survival_df$transformed_STAB2, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$STAB2_25p[survival_df$transformed_STAB2 > quantile(survival_df$transformed_STAB2, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$STAB2_75p <- NA
survival_df$STAB2_75p[survival_df$transformed_STAB2 <= quantile(survival_df$transformed_STAB2, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$STAB2_75p[survival_df$transformed_STAB2 > quantile(survival_df$transformed_STAB2, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$STAB2_66p <- NA
survival_df$STAB2_66p[survival_df$transformed_STAB2 <= quantile(survival_df$transformed_STAB2, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$STAB2_66p[survival_df$transformed_STAB2 > quantile(survival_df$transformed_STAB2, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$STAB2_33p <- NA
survival_df$STAB2_33p[survival_df$transformed_STAB2 <= quantile(survival_df$transformed_STAB2, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$STAB2_33p[survival_df$transformed_STAB2 > quantile(survival_df$transformed_STAB2, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$STAB2_33_66 <- NA
survival_df$STAB2_33_66[survival_df$transformed_STAB2 <= quantile(survival_df$transformed_STAB2, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$STAB2_33_66[survival_df$transformed_STAB2 >= quantile(survival_df$transformed_STAB2, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - STAB2
for(strata in c("STAB2_median", "STAB2_25p", "STAB2_75p", "STAB2_66p", "STAB2_33p", "STAB2_33_66")){
  
  H_gTxt <- paste("HIGH STAB2", sep="")
  L_gTxt <- paste("LOW STAB2", sep="")
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
  
  tiff(filename = paste("survival_plots/downreg/STAB2/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/downreg/STAB2/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/downreg/STAB2/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}




# Cox regression
meta <- c("Row.names", "X_PATIENT")

# Median
strata <- c("transformed_GDF2", "GDF2_median", "transformed_BMP10", "BMP10_median", "transformed_CLEC4M", "CLEC4M_median", "transformed_CLEC1B", "CLEC1B_median", "transformed_MT1H", "MT1H_median",
            "transformed_MARCO", "MARCO_median", "transformed_FCN2", "FCN2_median", "transformed_CXCL14", "CXCL14_median", "transformed_SFRP5", "SFRP5_median", "transformed_STAB2", "STAB2_median"
            )
# OS
os <- c("OS", "OS_months")
os_median_df <- survival_df %>% 
  select(all_of(c(meta, os, strata)))
os_median_df <- subset(os_median_df,  !is.na(OS) & !is.na(OS_months))
cox <- coxph(Surv(OS_months, OS) ~ GDF2_median + BMP10_median + CLEC4M_median + CLEC1B_median + MT1H_median + MARCO_median + FCN2_median + CXCL14_median + SFRP5_median + STAB2_median, data= os_median_df)
summary(cox)
# Remove STAB2
cox <- coxph(Surv(OS_months, OS) ~ GDF2_median + BMP10_median + CLEC4M_median + CLEC1B_median + MT1H_median + MARCO_median + FCN2_median + CXCL14_median + SFRP5_median, data= os_median_df)
summary(cox)
# Remove CXCL14
cox <- coxph(Surv(OS_months, OS) ~ GDF2_median + BMP10_median + CLEC4M_median + CLEC1B_median + MT1H_median + MARCO_median + FCN2_median + SFRP5_median, data= os_median_df)
summary(cox)
# Remove MARCO
cox <- coxph(Surv(OS_months, OS) ~ GDF2_median + BMP10_median + CLEC4M_median + CLEC1B_median + MT1H_median + FCN2_median + SFRP5_median, data= os_median_df)
summary(cox)
# Remove GDF2
cox <- coxph(Surv(OS_months, OS) ~ BMP10_median + CLEC4M_median + CLEC1B_median + MT1H_median + FCN2_median, data= os_median_df)
summary(cox)
# Remove BMP10
cox <- coxph(Surv(OS_months, OS) ~ CLEC4M_median + CLEC1B_median + MT1H_median + FCN2_median, data= os_median_df)
summary(cox)
# Remove MT1H
cox <- coxph(Surv(OS_months, OS) ~ CLEC4M_median + CLEC1B_median + FCN2_median, data= os_median_df)
summary(cox)
# Remove FCN2
cox <- coxph(Surv(OS_months, OS) ~ CLEC4M_median + CLEC1B_median, data= os_median_df)
summary(cox)
# Remove CLEC4m
cox <- coxph(Surv(OS_months, OS) ~ CLEC1B_median, data= os_median_df)
summary(cox)

# DSS
dss <- c("DSS", "DSS_months")
dss_median_df <- survival_df %>% 
  select(all_of(c(meta, dss, strata)))
dss_median_df <- subset(dss_median_df,  !is.na(DSS) & !is.na(DSS_months))
cox <- coxph(Surv(DSS_months, DSS) ~ GDF2_median + BMP10_median + CLEC4M_median + CLEC1B_median + MT1H_median + MARCO_median + FCN2_median + CXCL14_median + SFRP5_median + STAB2_median, data= dss_median_df)
summary(cox)
# Remove MARCO
cox <- coxph(Surv(DSS_months, DSS) ~ GDF2_median + BMP10_median + CLEC4M_median + CLEC1B_median + MT1H_median + FCN2_median + CXCL14_median + SFRP5_median + STAB2_median, data= dss_median_df)
summary(cox)
# Remove FCN2
cox <- coxph(Surv(DSS_months, DSS) ~ GDF2_median + BMP10_median + CLEC4M_median + CLEC1B_median + MT1H_median + CXCL14_median + SFRP5_median + STAB2_median, data= dss_median_df)
summary(cox)
# Remove MT1H
cox <- coxph(Surv(DSS_months, DSS) ~ GDF2_median + BMP10_median + CLEC4M_median + CLEC1B_median + CXCL14_median + SFRP5_median + STAB2_median, data= dss_median_df)
summary(cox)
# Remove SFRP
cox <- coxph(Surv(DSS_months, DSS) ~ GDF2_median + BMP10_median + CLEC4M_median + CLEC1B_median + CXCL14_median + STAB2_median, data= dss_median_df)
summary(cox)
# Remove Stab
cox <- coxph(Surv(DSS_months, DSS) ~ GDF2_median + BMP10_median + CLEC4M_median + CLEC1B_median + CXCL14_median, data= dss_median_df)
summary(cox)
# Remove CLEC4M
cox <- coxph(Surv(DSS_months, DSS) ~ GDF2_median + BMP10_median + CLEC1B_median + CXCL14_median, data= dss_median_df)
summary(cox)
# Remove CXLC14
cox <- coxph(Surv(DSS_months, DSS) ~ GDF2_median + BMP10_median + CLEC1B_median, data= dss_median_df)
summary(cox)
# Remove BMP
cox <- coxph(Surv(DSS_months, DSS) ~ GDF2_median + CLEC1B_median, data= dss_median_df)
summary(cox)
# Remove GDF2
cox <- coxph(Surv(DSS_months, DSS) ~ CLEC1B_median, data= dss_median_df)
summary(cox)


# PFI
pfi <- c("PFI", "PFI_months")
pfi_median_df <- survival_df %>% 
  select(all_of(c(meta, pfi, strata)))
pfi_median_df <- subset(pfi_median_df,  !is.na(PFI) & !is.na(PFI_months))
cox <- coxph(Surv(PFI_months, PFI) ~ GDF2_median + BMP10_median + CLEC4M_median + CLEC1B_median + MT1H_median + MARCO_median + FCN2_median + CXCL14_median + SFRP5_median + STAB2_median, data= pfi_median_df)
summary(cox)
# Remove FCN2
cox <- coxph(Surv(PFI_months, PFI) ~ GDF2_median + BMP10_median + CLEC4M_median + CLEC1B_median + MT1H_median + MARCO_median + CXCL14_median + SFRP5_median + STAB2_median, data= pfi_median_df)
summary(cox)
# Remove MT1H
cox <- coxph(Surv(PFI_months, PFI) ~ GDF2_median + BMP10_median + CLEC4M_median + CLEC1B_median + MARCO_median + CXCL14_median + SFRP5_median + STAB2_median, data= pfi_median_df)
summary(cox)
# Remove STAB
cox <- coxph(Surv(PFI_months, PFI) ~ GDF2_median + BMP10_median + CLEC4M_median + CLEC1B_median + MARCO_median + CXCL14_median + SFRP5_median, data= pfi_median_df)
summary(cox)
# Remove BMP10
cox <- coxph(Surv(PFI_months, PFI) ~ GDF2_median + CLEC4M_median + CLEC1B_median + MARCO_median + CXCL14_median + SFRP5_median, data= pfi_median_df)
summary(cox)
# Remove GDF2
cox <- coxph(Surv(PFI_months, PFI) ~ CLEC4M_median + CLEC1B_median + MARCO_median + CXCL14_median + SFRP5_median, data= pfi_median_df)
summary(cox)
# Remove CXCL14
cox <- coxph(Surv(PFI_months, PFI) ~ CLEC4M_median + CLEC1B_median + MARCO_median + SFRP5_median, data= pfi_median_df)
summary(cox)
# Remove SFRP5
cox <- coxph(Surv(PFI_months, PFI) ~ CLEC4M_median + CLEC1B_median + MARCO_median, data= pfi_median_df)
summary(cox)
# Remove MARCO
cox <- coxph(Surv(PFI_months, PFI) ~ CLEC4M_median + CLEC1B_median, data= pfi_median_df)
summary(cox)
# Remove CLECL4m
cox <- coxph(Surv(PFI_months, PFI) ~ CLEC1B_median, data= pfi_median_df)
summary(cox)

