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
# The last 3 are for CHCHD2 (identified as sig upreg in paired test) and pseudogenes CHCHD2P2/9
# The next last 3 are for RPL41 and its pseudogenes (pseudogene upreg in cancer vs healthy)
upreg_genes <- c("ENSG00000175063.17",
                 "ENSG00000164283.13",
                 "ENSG00000147257.15",
                 "ENSG00000204291.11",
                 "ENSG00000117399.14",
                 "ENSG00000089685.15",
                 "ENSG00000101412.13",
                 "ENSG00000164611.13",
                 "ENSG00000110492.15",
                 "ENSG00000131747.15",
                 "ENSG00000106153.13",
                 "ENSG00000215006.4",
                 "ENSG00000186940.6",
                 "ENSG00000229117.9",
                 "ENSG00000227063.5",
                 "ENSG00000256393.1",
                 "ENSG00000107159.13"
                 )


# Filter for desired genes HERE!!!
# Filter rows for paired tcga dgea downregulated genes
filtered_counts <- cancer_samples[cancer_samples$Ensembl_ID %in% upreg_genes, ]

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
names(survival_df)[names(survival_df) == "ENSG00000175063.17"] <- "UBE2C"
names(survival_df)[names(survival_df) == "ENSG00000164283.13"] <- "ESM1"
names(survival_df)[names(survival_df) == "ENSG00000147257.15"] <- "GPC3"
names(survival_df)[names(survival_df) == "ENSG00000204291.11"] <- "COL15A1"
names(survival_df)[names(survival_df) == "ENSG00000117399.14"] <- "CDC20"
names(survival_df)[names(survival_df) == "ENSG00000089685.15"] <- "BIRC5"
names(survival_df)[names(survival_df) == "ENSG00000101412.13"] <- "E2F1"
names(survival_df)[names(survival_df) == "ENSG00000164611.13"] <- "PTTG1"
names(survival_df)[names(survival_df) == "ENSG00000110492.15"] <- "MDK"
names(survival_df)[names(survival_df) == "ENSG00000131747.15"] <- "TOP2A"
names(survival_df)[names(survival_df) == "ENSG00000107159.13"] <- "CA9"
# CHCHD2 and its pseudogenes
names(survival_df)[names(survival_df) == "ENSG00000106153.13"] <- "CHCHD2"
names(survival_df)[names(survival_df) == "ENSG00000215006.4"] <- "CHCHD2P2"
names(survival_df)[names(survival_df) == "ENSG00000186940.6"] <- "CHCHD2P9"
# RPL41 and its pseudogenes
names(survival_df)[names(survival_df) == "ENSG00000229117.9"] <- "RPL41"
names(survival_df)[names(survival_df) == "ENSG00000227063.5"] <- "RPL41P1"
names(survival_df)[names(survival_df) == "ENSG00000256393.1"] <- "RPL41P5"


# Apply function to new transformed gene count column
transform_tpm_counts <- function(count){
  transformed_count <- (2^count) - 0.001
  return(transformed_count)
}

survival_df["transformed_UBE2C"] <- lapply(survival_df["UBE2C"], transform_tpm_counts)
survival_df["transformed_ESM1"] <- lapply(survival_df["ESM1"], transform_tpm_counts)
survival_df["transformed_GPC3"] <- lapply(survival_df["GPC3"], transform_tpm_counts)
survival_df["transformed_COL15A1"] <- lapply(survival_df["COL15A1"], transform_tpm_counts)
survival_df["transformed_CDC20"] <- lapply(survival_df["CDC20"], transform_tpm_counts)
survival_df["transformed_BIRC5"] <- lapply(survival_df["BIRC5"], transform_tpm_counts)
survival_df["transformed_E2F1"] <- lapply(survival_df["E2F1"], transform_tpm_counts)
survival_df["transformed_PTTG1"] <- lapply(survival_df["PTTG1"], transform_tpm_counts)
survival_df["transformed_MDK"] <- lapply(survival_df["MDK"], transform_tpm_counts)
survival_df["transformed_TOP2A"] <- lapply(survival_df["TOP2A"], transform_tpm_counts)
survival_df["transformed_CHCHD2"] <- lapply(survival_df["CHCHD2"], transform_tpm_counts)
survival_df["transformed_CHCHD2P2"] <- lapply(survival_df["CHCHD2P2"], transform_tpm_counts)
survival_df["transformed_CHCHD2P9"] <- lapply(survival_df["CHCHD2P9"], transform_tpm_counts)
survival_df["transformed_RPL41"] <- lapply(survival_df["RPL41"], transform_tpm_counts)
survival_df["transformed_RPL41P1"] <- lapply(survival_df["RPL41P1"], transform_tpm_counts)
survival_df["transformed_RPL41P5"] <- lapply(survival_df["RPL41P5"], transform_tpm_counts)
survival_df["transformed_CA9"] <- lapply(survival_df["CA9"], transform_tpm_counts)

# Convert days to months
# 1 day = 0.0328767 months
days_to_months_converter <- function(no_days){
  no_months <- round(no_days*0.0328767, 1)
  return(no_months)
}

survival_df["OS_months"] <- lapply(survival_df["OS.time"], days_to_months_converter)
survival_df["DSS_months"] <- lapply(survival_df["DSS.time"], days_to_months_converter)
survival_df["PFI_months"] <- lapply(survival_df["PFI.time"], days_to_months_converter)


# UBE2C
# Make survival plot
survival_df$UBE2C_median <- NA
survival_df$UBE2C_median[survival_df$transformed_UBE2C <= quantile(survival_df$transformed_UBE2C, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$UBE2C_median[survival_df$transformed_UBE2C > quantile(survival_df$transformed_UBE2C, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$UBE2C_25p <- NA
survival_df$UBE2C_25p[survival_df$transformed_UBE2C <= quantile(survival_df$transformed_UBE2C, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$UBE2C_25p[survival_df$transformed_UBE2C > quantile(survival_df$transformed_UBE2C, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$UBE2C_75p <- NA
survival_df$UBE2C_75p[survival_df$transformed_UBE2C <= quantile(survival_df$transformed_UBE2C, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$UBE2C_75p[survival_df$transformed_UBE2C > quantile(survival_df$transformed_UBE2C, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$UBE2C_66p <- NA
survival_df$UBE2C_66p[survival_df$transformed_UBE2C <= quantile(survival_df$transformed_UBE2C, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$UBE2C_66p[survival_df$transformed_UBE2C > quantile(survival_df$transformed_UBE2C, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$UBE2C_33p <- NA
survival_df$UBE2C_33p[survival_df$transformed_UBE2C <= quantile(survival_df$transformed_UBE2C, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$UBE2C_33p[survival_df$transformed_UBE2C > quantile(survival_df$transformed_UBE2C, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$UBE2C_33_66 <- NA
survival_df$UBE2C_33_66[survival_df$transformed_UBE2C <= quantile(survival_df$transformed_UBE2C, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$UBE2C_33_66[survival_df$transformed_UBE2C >= quantile(survival_df$transformed_UBE2C, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - UBE2C
for(strata in c("UBE2C_median", "UBE2C_25p", "UBE2C_75p", "UBE2C_66p", "UBE2C_33p", "UBE2C_33_66")){
  
  H_gTxt <- paste("HIGH UBE2C", sep="")
  L_gTxt <- paste("LOW UBE2C", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/upreg/UBE2C/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/upreg/UBE2C/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/upreg/UBE2C/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

# GPC3
# Make survival plot
survival_df$GPC3_median <- NA
survival_df$GPC3_median[survival_df$transformed_GPC3 <= quantile(survival_df$transformed_GPC3, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$GPC3_median[survival_df$transformed_GPC3 > quantile(survival_df$transformed_GPC3, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$GPC3_25p <- NA
survival_df$GPC3_25p[survival_df$transformed_GPC3 <= quantile(survival_df$transformed_GPC3, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$GPC3_25p[survival_df$transformed_GPC3 > quantile(survival_df$transformed_GPC3, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$GPC3_75p <- NA
survival_df$GPC3_75p[survival_df$transformed_GPC3 <= quantile(survival_df$transformed_GPC3, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$GPC3_75p[survival_df$transformed_GPC3 > quantile(survival_df$transformed_GPC3, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$GPC3_66p <- NA
survival_df$GPC3_66p[survival_df$transformed_GPC3 <= quantile(survival_df$transformed_GPC3, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$GPC3_66p[survival_df$transformed_GPC3 > quantile(survival_df$transformed_GPC3, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$GPC3_33p <- NA
survival_df$GPC3_33p[survival_df$transformed_GPC3 <= quantile(survival_df$transformed_GPC3, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$GPC3_33p[survival_df$transformed_GPC3 > quantile(survival_df$transformed_GPC3, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$GPC3_33_66 <- NA
survival_df$GPC3_33_66[survival_df$transformed_GPC3 <= quantile(survival_df$transformed_GPC3, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$GPC3_33_66[survival_df$transformed_GPC3 >= quantile(survival_df$transformed_GPC3, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - GPC3
for(strata in c("GPC3_median", "GPC3_25p", "GPC3_75p", "GPC3_66p", "GPC3_33p", "GPC3_33_66")){
  
  H_gTxt <- paste("HIGH GPC3", sep="")
  L_gTxt <- paste("LOW GPC3", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/upreg/GPC3/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/upreg/GPC3/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/upreg/GPC3/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

# COL15A1
# Make survival plot
survival_df$COL15A1_median <- NA
survival_df$COL15A1_median[survival_df$transformed_COL15A1 <= quantile(survival_df$transformed_COL15A1, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$COL15A1_median[survival_df$transformed_COL15A1 > quantile(survival_df$transformed_COL15A1, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$COL15A1_25p <- NA
survival_df$COL15A1_25p[survival_df$transformed_COL15A1 <= quantile(survival_df$transformed_COL15A1, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$COL15A1_25p[survival_df$transformed_COL15A1 > quantile(survival_df$transformed_COL15A1, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$COL15A1_75p <- NA
survival_df$COL15A1_75p[survival_df$transformed_COL15A1 <= quantile(survival_df$transformed_COL15A1, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$COL15A1_75p[survival_df$transformed_COL15A1 > quantile(survival_df$transformed_COL15A1, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$COL15A1_66p <- NA
survival_df$COL15A1_66p[survival_df$transformed_COL15A1 <= quantile(survival_df$transformed_COL15A1, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$COL15A1_66p[survival_df$transformed_COL15A1 > quantile(survival_df$transformed_COL15A1, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$COL15A1_33p <- NA
survival_df$COL15A1_33p[survival_df$transformed_COL15A1 <= quantile(survival_df$transformed_COL15A1, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$COL15A1_33p[survival_df$transformed_COL15A1 > quantile(survival_df$transformed_COL15A1, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$COL15A1_33_66 <- NA
survival_df$COL15A1_33_66[survival_df$transformed_COL15A1 <= quantile(survival_df$transformed_COL15A1, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$COL15A1_33_66[survival_df$transformed_COL15A1 >= quantile(survival_df$transformed_COL15A1, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - COL15A1
for(strata in c("COL15A1_median", "COL15A1_25p", "COL15A1_75p", "COL15A1_66p", "COL15A1_33p", "COL15A1_33_66")){
  
  H_gTxt <- paste("HIGH COL15A1", sep="")
  L_gTxt <- paste("LOW COL15A1", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/upreg/COL15A1/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/upreg/COL15A1/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/upreg/COL15A1/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

# CDC20
# Make survival plot
survival_df$CDC20_median <- NA
survival_df$CDC20_median[survival_df$transformed_CDC20 <= quantile(survival_df$transformed_CDC20, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$CDC20_median[survival_df$transformed_CDC20 > quantile(survival_df$transformed_CDC20, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$CDC20_25p <- NA
survival_df$CDC20_25p[survival_df$transformed_CDC20 <= quantile(survival_df$transformed_CDC20, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$CDC20_25p[survival_df$transformed_CDC20 > quantile(survival_df$transformed_CDC20, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$CDC20_75p <- NA
survival_df$CDC20_75p[survival_df$transformed_CDC20 <= quantile(survival_df$transformed_CDC20, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$CDC20_75p[survival_df$transformed_CDC20 > quantile(survival_df$transformed_CDC20, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$CDC20_66p <- NA
survival_df$CDC20_66p[survival_df$transformed_CDC20 <= quantile(survival_df$transformed_CDC20, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$CDC20_66p[survival_df$transformed_CDC20 > quantile(survival_df$transformed_CDC20, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$CDC20_33p <- NA
survival_df$CDC20_33p[survival_df$transformed_CDC20 <= quantile(survival_df$transformed_CDC20, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CDC20_33p[survival_df$transformed_CDC20 > quantile(survival_df$transformed_CDC20, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$CDC20_33_66 <- NA
survival_df$CDC20_33_66[survival_df$transformed_CDC20 <= quantile(survival_df$transformed_CDC20, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CDC20_33_66[survival_df$transformed_CDC20 >= quantile(survival_df$transformed_CDC20, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - CDC20
for(strata in c("CDC20_median", "CDC20_25p", "CDC20_75p", "CDC20_66p", "CDC20_33p", "CDC20_33_66")){
  
  H_gTxt <- paste("HIGH CDC20", sep="")
  L_gTxt <- paste("LOW CDC20", sep="")
  suffix <- paste(datasetName, strata, sep="_")
  
  sc = "OS"
  tmp.sub <- subset(survival_df, !is.na(OS) & !is.na(OS_months))
  tmp.sub$usedSubset <- as.factor(tmp.sub[[strata]])
  tmp.sub <- subset(tmp.sub,!is.na(usedSubset))
  
  fit <- survfit(Surv(OS_months, OS) ~ usedSubset, data= tmp.sub)
  p <-  ggsurvplot(fit, palette = c("#ED0000FF", "black"), title=paste(strata,sep=" - "), xlab="Months", ylab="Overall Survival (OS)", 
                   break.time.by = 30,
                   xlim = c(0, 125),
                   font.x = c(20, "black"),
                   font.y = c(20, "black"),
                   font.tickslab = c(20, "plain", "black"),
                   legend = c(0.3, 0.15),
                   legend.title = "",
                   legend.labs =c(paste0 (H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                  paste0 (L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                   font.legend = c(20, "black"),
                   fun = function(y) y*100,
                   pval = TRUE, pval.coord=c(70, 80), pval.size = 7,
                   size = 1, censor.size = 2)
  
  tiff(filename = paste("survival_plots/removed_sample/upreg/CDC20/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
                  font.x = c(20, "black"),
                  font.y = c(20, "black"),
                  font.tickslab = c(18, "plain", "black"),
                  legend = c(0.3, 0.15),
                  legend.title = "",
                  legend.labs =c(paste0 (H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                 paste0 (L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                  font.legend = c(20, "black"),
                  fun = function(y) y*100,
                  pval = TRUE, pval.coord = c(70, 80), pval.size = 7,
                  size = 1, censor.size = 2)
  
  tiff(filename = paste ("survival_plots/removed_sample/upreg/CDC20/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
                   font.x = c(20, "black"),
                   font.y = c(20, "black"),
                   font.tickslab = c(20, "plain", "black"),
                   legend = c(0.3, 0.15),
                   legend.title = "",
                   legend.labs = c(paste0(H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                   paste0(L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                   font.legend = c(18, "black"),
                   fun = function(y) y*100,
                   pval = TRUE, pval.coord = c(70, 80), pval.size = 7,
                   size = 1, censor.size = 2)
  
  tiff(filename =  paste("survival_plots/removed_sample/upreg/CDC20/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

# BIRC5
# Make survival plot
survival_df$BIRC5_median <- NA
survival_df$BIRC5_median[survival_df$transformed_BIRC5 <= quantile(survival_df$transformed_BIRC5, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$BIRC5_median[survival_df$transformed_BIRC5 > quantile(survival_df$transformed_BIRC5, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$BIRC5_25p <- NA
survival_df$BIRC5_25p[survival_df$transformed_BIRC5 <= quantile(survival_df$transformed_BIRC5, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$BIRC5_25p[survival_df$transformed_BIRC5 > quantile(survival_df$transformed_BIRC5, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$BIRC5_75p <- NA
survival_df$BIRC5_75p[survival_df$transformed_BIRC5 <= quantile(survival_df$transformed_BIRC5, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$BIRC5_75p[survival_df$transformed_BIRC5 > quantile(survival_df$transformed_BIRC5, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$BIRC5_66p <- NA
survival_df$BIRC5_66p[survival_df$transformed_BIRC5 <= quantile(survival_df$transformed_BIRC5, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$BIRC5_66p[survival_df$transformed_BIRC5 > quantile(survival_df$transformed_BIRC5, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$BIRC5_33p <- NA
survival_df$BIRC5_33p[survival_df$transformed_BIRC5 <= quantile(survival_df$transformed_BIRC5, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$BIRC5_33p[survival_df$transformed_BIRC5 > quantile(survival_df$transformed_BIRC5, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$BIRC5_33_66 <- NA
survival_df$BIRC5_33_66[survival_df$transformed_BIRC5 <= quantile(survival_df$transformed_BIRC5, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$BIRC5_33_66[survival_df$transformed_BIRC5 >= quantile(survival_df$transformed_BIRC5, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - BIRC5
for(strata in c("BIRC5_median", "BIRC5_25p", "BIRC5_75p", "BIRC5_66p", "BIRC5_33p", "BIRC5_33_66")){
  
  H_gTxt <- paste("HIGH BIRC5", sep="")
  L_gTxt <- paste("LOW BIRC5", sep="")
  suffix <- paste(datasetName, strata, sep="_")
  
  sc = "OS"
  tmp.sub <- subset(survival_df, !is.na(OS) & !is.na(OS_months))
  tmp.sub$usedSubset <- as.factor(tmp.sub[[strata]])
  tmp.sub <- subset(tmp.sub,!is.na(usedSubset))
  
  fit <- survfit(Surv(OS_months, OS) ~ usedSubset, data= tmp.sub)
  p <-  ggsurvplot(fit, palette = c("#ED0000FF", "black"), title=paste(strata,sep=" - "), xlab="Months", ylab="Overall Survival (OS)", 
                   break.time.by = 30,
                   xlim = c(0, 125),
                   font.x = c(20, "black"),
                   font.y = c(20, "black"),
                   font.tickslab = c(18, "plain", "black"),
                   legend = c(0.3, 0.15),
                   legend.title = "",
                   legend.labs =c(paste0 (H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                  paste0 (L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                   font.legend = c(20, "black"),
                   fun = function(y) y*100,
                   pval = TRUE, pval.coord=c(70, 80), pval.size = 7,
                   size = 1, censor.size = 2)
  
  tiff(filename = paste("survival_plots/removed_sample/upreg/BIRC5/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
                  font.x = c(20, "black"),
                  font.y = c(20, "black"),
                  font.tickslab = c(18, "plain", "black"),
                  legend = c(0.3, 0.15),
                  legend.title = "",
                  legend.labs =c(paste0 (H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                 paste0 (L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                  font.legend = c(20, "black"),
                  fun = function(y) y*100,
                  pval = TRUE, pval.coord = c(70, 80), pval.size = 7,
                  size = 1, censor.size = 2)
  
  tiff(filename = paste ("survival_plots/removed_sample/upreg/BIRC5/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
                   font.x = c(20, "black"),
                   font.y = c(20, "black"),
                   font.tickslab = c(18, "plain", "black"),
                   legend = c(0.3, 0.15),
                   legend.title = "",
                   legend.labs = c(paste0(H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                   paste0(L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                   font.legend = c(20, "black"),
                   fun = function(y) y*100,
                   pval = TRUE, pval.coord = c(70, 80), pval.size = 7,
                   size = 1, censor.size = 2)
  
  tiff(filename =  paste("survival_plots/removed_sample/upreg/BIRC5/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}


# E2F1
# Make survival plot
survival_df$E2F1_median <- NA
survival_df$E2F1_median[survival_df$transformed_E2F1 <= quantile(survival_df$transformed_E2F1, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$E2F1_median[survival_df$transformed_E2F1 > quantile(survival_df$transformed_E2F1, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$E2F1_25p <- NA
survival_df$E2F1_25p[survival_df$transformed_E2F1 <= quantile(survival_df$transformed_E2F1, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$E2F1_25p[survival_df$transformed_E2F1 > quantile(survival_df$transformed_E2F1, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$E2F1_75p <- NA
survival_df$E2F1_75p[survival_df$transformed_E2F1 <= quantile(survival_df$transformed_E2F1, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$E2F1_75p[survival_df$transformed_E2F1 > quantile(survival_df$transformed_E2F1, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$E2F1_66p <- NA
survival_df$E2F1_66p[survival_df$transformed_E2F1 <= quantile(survival_df$transformed_E2F1, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$E2F1_66p[survival_df$transformed_E2F1 > quantile(survival_df$transformed_E2F1, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$E2F1_33p <- NA
survival_df$E2F1_33p[survival_df$transformed_E2F1 <= quantile(survival_df$transformed_E2F1, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$E2F1_33p[survival_df$transformed_E2F1 > quantile(survival_df$transformed_E2F1, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$E2F1_33_66 <- NA
survival_df$E2F1_33_66[survival_df$transformed_E2F1 <= quantile(survival_df$transformed_E2F1, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$E2F1_33_66[survival_df$transformed_E2F1 >= quantile(survival_df$transformed_E2F1, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - E2F1
for(strata in c("E2F1_median", "E2F1_25p", "E2F1_75p", "E2F1_66p", "E2F1_33p", "E2F1_33_66")){
  
  H_gTxt <- paste("HIGH E2F1", sep="")
  L_gTxt <- paste("LOW E2F1", sep="")
  suffix <- paste(datasetName, strata, sep="_")
  
  sc = "OS"
  tmp.sub <- subset(survival_df, !is.na(OS) & !is.na(OS_months))
  tmp.sub$usedSubset <- as.factor(tmp.sub[[strata]])
  tmp.sub <- subset(tmp.sub,!is.na(usedSubset))
  
  fit <- survfit(Surv(OS_months, OS) ~ usedSubset, data= tmp.sub)
  p <-  ggsurvplot(fit, palette = c("#ED0000FF", "black"), title=paste(strata,sep=" - "), xlab="Months", ylab="Overall Survival (OS)", 
                   break.time.by = 30,
                   xlim = c(0, 125),
                   font.x = c(20, "black"),
                   font.y = c(20, "black"),
                   font.tickslab = c(20, "plain", "black"),
                   legend = c(0.3, 0.15),
                   legend.title = "",
                   legend.labs =c(paste0 (H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                  paste0 (L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                   font.legend = c(20, "black"),
                   fun = function(y) y*100,
                   pval = TRUE, pval.coord=c(70, 80), pval.size = 7,
                   size = 1, censor.size = 2)
  
  tiff(filename = paste("survival_plots/removed_sample/upreg/E2F1/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
                  font.x = c(20, "black"),
                  font.y = c(20, "black"),
                  font.tickslab = c(20, "plain", "black"),
                  legend = c(0.3, 0.15),
                  legend.title = "",
                  legend.labs =c(paste0 (H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                 paste0 (L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                  font.legend = c(20, "black"),
                  fun = function(y) y*100,
                  pval = TRUE, pval.coord = c(70, 80), pval.size = 7,
                  size = 1, censor.size = 2)
  
  tiff(filename = paste ("survival_plots/removed_sample/upreg/E2F1/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
                   font.x = c(20, "black"),
                   font.y = c(20, "black"),
                   font.tickslab = c(20, "plain", "black"),
                   legend = c(0.3, 0.15),
                   legend.title = "",
                   legend.labs = c(paste0(H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                   paste0(L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                   font.legend = c(18, "black"),
                   fun = function(y) y*100,
                   pval = TRUE, pval.coord = c(70, 80), pval.size = 7,
                   size = 1, censor.size = 2)
  
  tiff(filename =  paste("survival_plots/removed_sample/upreg/E2F1/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}


# PTTG1
# Make survival plot
survival_df$PTTG1_median <- NA
survival_df$PTTG1_median[survival_df$transformed_PTTG1 <= quantile(survival_df$transformed_PTTG1, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$PTTG1_median[survival_df$transformed_PTTG1 > quantile(survival_df$transformed_PTTG1, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$PTTG1_25p <- NA
survival_df$PTTG1_25p[survival_df$transformed_PTTG1 <= quantile(survival_df$transformed_PTTG1, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$PTTG1_25p[survival_df$transformed_PTTG1 > quantile(survival_df$transformed_PTTG1, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$PTTG1_75p <- NA
survival_df$PTTG1_75p[survival_df$transformed_PTTG1 <= quantile(survival_df$transformed_PTTG1, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$PTTG1_75p[survival_df$transformed_PTTG1 > quantile(survival_df$transformed_PTTG1, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$PTTG1_66p <- NA
survival_df$PTTG1_66p[survival_df$transformed_PTTG1 <= quantile(survival_df$transformed_PTTG1, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$PTTG1_66p[survival_df$transformed_PTTG1 > quantile(survival_df$transformed_PTTG1, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$PTTG1_33p <- NA
survival_df$PTTG1_33p[survival_df$transformed_PTTG1 <= quantile(survival_df$transformed_PTTG1, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$PTTG1_33p[survival_df$transformed_PTTG1 > quantile(survival_df$transformed_PTTG1, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$PTTG1_33_66 <- NA
survival_df$PTTG1_33_66[survival_df$transformed_PTTG1 <= quantile(survival_df$transformed_PTTG1, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$PTTG1_33_66[survival_df$transformed_PTTG1 >= quantile(survival_df$transformed_PTTG1, 0.66, na.rm=TRUE)] <- "HIGH"

# SURVIVAL CURVES - PTTG1
for(strata in c("PTTG1_median", "PTTG1_25p", "PTTG1_75p", "PTTG1_66p", "PTTG1_33p", "PTTG1_33_66")){
  
  H_gTxt <- paste("HIGH PTTG1", sep="")
  L_gTxt <- paste("LOW PTTG1", sep="")
  suffix <- paste(datasetName, strata, sep="_")
  
  sc = "OS"
  tmp.sub <- subset(survival_df, !is.na(OS) & !is.na(OS_months))
  tmp.sub$usedSubset <- as.factor(tmp.sub[[strata]])
  tmp.sub <- subset(tmp.sub,!is.na(usedSubset))
  
  fit <- survfit(Surv(OS_months, OS) ~ usedSubset, data= tmp.sub)
  p <-  ggsurvplot(fit, palette = c("#ED0000FF", "black"), title=paste(strata,sep=" - "), xlab="Months", ylab="Overall Survival (OS)", 
                   break.time.by = 30,
                   xlim = c(0, 125),
                   font.x = c(20, "black"),
                   font.y = c(20, "black"),
                   font.tickslab = c(20, "plain", "black"),
                   legend = c(0.3, 0.15),
                   legend.title = "",
                   legend.labs =c(paste0 (H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                  paste0 (L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                   font.legend = c(20, "black"),
                   fun = function(y) y*100,
                   pval = TRUE, pval.coord=c(70, 80), pval.size = 7,
                   size = 1, censor.size = 2)
  
  tiff(filename = paste("survival_plots/removed_sample/upreg/PTTG1/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
                  font.x = c(20, "black"),
                  font.y = c(20, "black"),
                  font.tickslab = c(20, "plain", "black"),
                  legend = c(0.3, 0.15),
                  legend.title = "",
                  legend.labs =c(paste0 (H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                 paste0 (L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                  font.legend = c(20, "black"),
                  fun = function(y) y*100,
                  pval = TRUE, pval.coord = c(70, 80), pval.size = 7,
                  size = 1, censor.size = 2)
  
  tiff(filename = paste ("survival_plots/removed_sample/upreg/PTTG1/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
                   font.x = c(20, "black"),
                   font.y = c(20, "black"),
                   font.tickslab = c(20, "plain", "black"),
                   legend = c(0.3, 0.15),
                   legend.title = "",
                   legend.labs = c(paste0(H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                   paste0(L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                   font.legend = c(20, "black"),
                   fun = function(y) y*100,
                   pval = TRUE, pval.coord = c(70, 80), pval.size = 7,
                   size = 1, censor.size = 2)
  
  tiff(filename =  paste("survival_plots/removed_sample/upreg/PTTG1/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}


# MDK
# Make survival plot
survival_df$MDK_median <- NA
survival_df$MDK_median[survival_df$transformed_MDK <= quantile(survival_df$transformed_MDK, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$MDK_median[survival_df$transformed_MDK > quantile(survival_df$transformed_MDK, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$MDK_25p <- NA
survival_df$MDK_25p[survival_df$transformed_MDK <= quantile(survival_df$transformed_MDK, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$MDK_25p[survival_df$transformed_MDK > quantile(survival_df$transformed_MDK, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$MDK_75p <- NA
survival_df$MDK_75p[survival_df$transformed_MDK <= quantile(survival_df$transformed_MDK, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$MDK_75p[survival_df$transformed_MDK > quantile(survival_df$transformed_MDK, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$MDK_66p <- NA
survival_df$MDK_66p[survival_df$transformed_MDK <= quantile(survival_df$transformed_MDK, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$MDK_66p[survival_df$transformed_MDK > quantile(survival_df$transformed_MDK, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$MDK_33p <- NA
survival_df$MDK_33p[survival_df$transformed_MDK <= quantile(survival_df$transformed_MDK, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$MDK_33p[survival_df$transformed_MDK > quantile(survival_df$transformed_MDK, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$MDK_33_66 <- NA
survival_df$MDK_33_66[survival_df$transformed_MDK <= quantile(survival_df$transformed_MDK, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$MDK_33_66[survival_df$transformed_MDK >= quantile(survival_df$transformed_MDK, 0.66, na.rm=TRUE)] <- "HIGH"

# SURVIVAL CURVES - MDK
for(strata in c("MDK_median", "MDK_25p", "MDK_75p", "MDK_66p", "MDK_33p", "MDK_33_66")){
  
  H_gTxt <- paste("HIGH MDK", sep="")
  L_gTxt <- paste("LOW MDK", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/upreg/MDK/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/upreg/MDK/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/upreg/MDK/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

# TOP2A
# Make survival plot
survival_df$TOP2A_median <- NA
survival_df$TOP2A_median[survival_df$transformed_TOP2A <= quantile(survival_df$transformed_TOP2A, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$TOP2A_median[survival_df$transformed_TOP2A > quantile(survival_df$transformed_TOP2A, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$TOP2A_25p <- NA
survival_df$TOP2A_25p[survival_df$transformed_TOP2A <= quantile(survival_df$transformed_TOP2A, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$TOP2A_25p[survival_df$transformed_TOP2A > quantile(survival_df$transformed_TOP2A, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$TOP2A_75p <- NA
survival_df$TOP2A_75p[survival_df$transformed_TOP2A <= quantile(survival_df$transformed_TOP2A, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$TOP2A_75p[survival_df$transformed_TOP2A > quantile(survival_df$transformed_TOP2A, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$TOP2A_66p <- NA
survival_df$TOP2A_66p[survival_df$transformed_TOP2A <= quantile(survival_df$transformed_TOP2A, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$TOP2A_66p[survival_df$transformed_TOP2A > quantile(survival_df$transformed_TOP2A, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$TOP2A_33p <- NA
survival_df$TOP2A_33p[survival_df$transformed_TOP2A <= quantile(survival_df$transformed_TOP2A, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$TOP2A_33p[survival_df$transformed_TOP2A > quantile(survival_df$transformed_TOP2A, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$TOP2A_33_66 <- NA
survival_df$TOP2A_33_66[survival_df$transformed_TOP2A <= quantile(survival_df$transformed_TOP2A, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$TOP2A_33_66[survival_df$transformed_TOP2A >= quantile(survival_df$transformed_TOP2A, 0.66, na.rm=TRUE)] <- "HIGH"

# SURVIVAL CURVES - TOP2A
for(strata in c("TOP2A_median", "TOP2A_25p", "TOP2A_75p", "TOP2A_66p", "TOP2A_33p", "TOP2A_33_66")){
  
  H_gTxt <- paste("HIGH TOP2A", sep="")
  L_gTxt <- paste("LOW TOP2A", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/upreg/TOP2A/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/upreg/TOP2A/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/upreg/TOP2A/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}



# Plot already made in previous survival script
# ESM1
# Make survival plot
survival_df$ESM1_median <- NA
survival_df$ESM1_median[survival_df$transformed_ESM1 <= quantile(survival_df$transformed_ESM1, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$ESM1_median[survival_df$transformed_ESM1 > quantile(survival_df$transformed_ESM1, 0.5, na.rm=TRUE)] <- "HIGH"


# CHCHD2
# Make survival plot
survival_df$CHCHD2_median <- NA
survival_df$CHCHD2_median[survival_df$transformed_CHCHD2 <= quantile(survival_df$transformed_CHCHD2, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2_median[survival_df$transformed_CHCHD2 > quantile(survival_df$transformed_CHCHD2, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$CHCHD2_25p <- NA
survival_df$CHCHD2_25p[survival_df$transformed_CHCHD2 <= quantile(survival_df$transformed_CHCHD2, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2_25p[survival_df$transformed_CHCHD2 > quantile(survival_df$transformed_CHCHD2, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$CHCHD2_75p <- NA
survival_df$CHCHD2_75p[survival_df$transformed_CHCHD2 <= quantile(survival_df$transformed_CHCHD2, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2_75p[survival_df$transformed_CHCHD2 > quantile(survival_df$transformed_CHCHD2, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$CHCHD2_66p <- NA
survival_df$CHCHD2_66p[survival_df$transformed_CHCHD2 <= quantile(survival_df$transformed_CHCHD2, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2_66p[survival_df$transformed_CHCHD2 > quantile(survival_df$transformed_CHCHD2, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$CHCHD2_33p <- NA
survival_df$CHCHD2_33p[survival_df$transformed_CHCHD2 <= quantile(survival_df$transformed_CHCHD2, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2_33p[survival_df$transformed_CHCHD2 > quantile(survival_df$transformed_CHCHD2, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$CHCHD2_33_66 <- NA
survival_df$CHCHD2_33_66[survival_df$transformed_CHCHD2 <= quantile(survival_df$transformed_CHCHD2, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2_33_66[survival_df$transformed_CHCHD2 >= quantile(survival_df$transformed_CHCHD2, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - CHCHD2
for(strata in c("CHCHD2_median", "CHCHD2_25p", "CHCHD2_75p", "CHCHD2_66p", "CHCHD2_33p", "CHCHD2_33_66")){
  
  H_gTxt <- paste("HIGH CHCHD2", sep="")
  L_gTxt <- paste("LOW CHCHD2", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/upreg/CHCHD2/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/upreg/CHCHD2/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/upreg/CHCHD2/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}


# CHCHD2P2
# Make survival plot
survival_df$CHCHD2P2_median <- NA
survival_df$CHCHD2P2_median[survival_df$transformed_CHCHD2P2 <= quantile(survival_df$transformed_CHCHD2P2, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2P2_median[survival_df$transformed_CHCHD2P2 > quantile(survival_df$transformed_CHCHD2P2, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$CHCHD2P2_25p <- NA
survival_df$CHCHD2P2_25p[survival_df$transformed_CHCHD2P2 <= quantile(survival_df$transformed_CHCHD2P2, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2P2_25p[survival_df$transformed_CHCHD2P2 > quantile(survival_df$transformed_CHCHD2P2, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$CHCHD2P2_75p <- NA
survival_df$CHCHD2P2_75p[survival_df$transformed_CHCHD2P2 <= quantile(survival_df$transformed_CHCHD2P2, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2P2_75p[survival_df$transformed_CHCHD2P2 > quantile(survival_df$transformed_CHCHD2P2, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$CHCHD2P2_66p <- NA
survival_df$CHCHD2P2_66p[survival_df$transformed_CHCHD2P2 <= quantile(survival_df$transformed_CHCHD2P2, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2P2_66p[survival_df$transformed_CHCHD2P2 > quantile(survival_df$transformed_CHCHD2P2, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$CHCHD2P2_33p <- NA
survival_df$CHCHD2P2_33p[survival_df$transformed_CHCHD2P2 <= quantile(survival_df$transformed_CHCHD2P2, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2P2_33p[survival_df$transformed_CHCHD2P2 > quantile(survival_df$transformed_CHCHD2P2, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$CHCHD2P2_33_66 <- NA
survival_df$CHCHD2P2_33_66[survival_df$transformed_CHCHD2P2 <= quantile(survival_df$transformed_CHCHD2P2, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2P2_33_66[survival_df$transformed_CHCHD2P2 >= quantile(survival_df$transformed_CHCHD2P2, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - CHCHD2P2
for(strata in c("CHCHD2P2_median", "CHCHD2P2_25p", "CHCHD2P2_75p", "CHCHD2P2_66p", "CHCHD2P2_33p", "CHCHD2P2_33_66")){
  
  H_gTxt <- paste("HIGH CHCHD2P2", sep="")
  L_gTxt <- paste("LOW CHCHD2P2", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/upreg/CHCHD2P2/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/upreg/CHCHD2P2/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/upreg/CHCHD2P2/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

# CHCHD2P9
# Make survival plot
survival_df$CHCHD2P9_median <- NA
survival_df$CHCHD2P9_median[survival_df$transformed_CHCHD2P9 <= quantile(survival_df$transformed_CHCHD2P9, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2P9_median[survival_df$transformed_CHCHD2P9 > quantile(survival_df$transformed_CHCHD2P9, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$CHCHD2P9_25p <- NA
survival_df$CHCHD2P9_25p[survival_df$transformed_CHCHD2P9 <= quantile(survival_df$transformed_CHCHD2P9, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2P9_25p[survival_df$transformed_CHCHD2P9 > quantile(survival_df$transformed_CHCHD2P9, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$CHCHD2P9_75p <- NA
survival_df$CHCHD2P9_75p[survival_df$transformed_CHCHD2P9 <= quantile(survival_df$transformed_CHCHD2P9, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2P9_75p[survival_df$transformed_CHCHD2P9 > quantile(survival_df$transformed_CHCHD2P9, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$CHCHD2P9_66p <- NA
survival_df$CHCHD2P9_66p[survival_df$transformed_CHCHD2P9 <= quantile(survival_df$transformed_CHCHD2P9, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2P9_66p[survival_df$transformed_CHCHD2P9 > quantile(survival_df$transformed_CHCHD2P9, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$CHCHD2P9_33p <- NA
survival_df$CHCHD2P9_33p[survival_df$transformed_CHCHD2P9 <= quantile(survival_df$transformed_CHCHD2P9, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2P9_33p[survival_df$transformed_CHCHD2P9 > quantile(survival_df$transformed_CHCHD2P9, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$CHCHD2P9_33_66 <- NA
survival_df$CHCHD2P9_33_66[survival_df$transformed_CHCHD2P9 <= quantile(survival_df$transformed_CHCHD2P9, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CHCHD2P9_33_66[survival_df$transformed_CHCHD2P9 >= quantile(survival_df$transformed_CHCHD2P9, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - CHCHD2P9
for(strata in c("CHCHD2P9_median", "CHCHD2P9_25p", "CHCHD2P9_75p", "CHCHD2P9_66p", "CHCHD2P9_33p", "CHCHD2P9_33_66")){
  
  H_gTxt <- paste("HIGH CHCHD2P9", sep="")
  L_gTxt <- paste("LOW CHCHD2P9", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/upreg/CHCHD2P9/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/upreg/CHCHD2P9/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/upreg/CHCHD2P9/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}




# RPL41 
# Make survival plot
survival_df$RPL41_median <- NA
survival_df$RPL41_median[survival_df$transformed_RPL41 <= quantile(survival_df$transformed_RPL41, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$RPL41_median[survival_df$transformed_RPL41 > quantile(survival_df$transformed_RPL41, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$RPL41_25p <- NA
survival_df$RPL41_25p[survival_df$transformed_RPL41 <= quantile(survival_df$transformed_RPL41, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$RPL41_25p[survival_df$transformed_RPL41 > quantile(survival_df$transformed_RPL41, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$RPL41_75p <- NA
survival_df$RPL41_75p[survival_df$transformed_RPL41 <= quantile(survival_df$transformed_RPL41, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$RPL41_75p[survival_df$transformed_RPL41 > quantile(survival_df$transformed_RPL41, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$RPL41_66p <- NA
survival_df$RPL41_66p[survival_df$transformed_RPL41 <= quantile(survival_df$transformed_RPL41, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$RPL41_66p[survival_df$transformed_RPL41 > quantile(survival_df$transformed_RPL41, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$RPL41_33p <- NA
survival_df$RPL41_33p[survival_df$transformed_RPL41 <= quantile(survival_df$transformed_RPL41, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$RPL41_33p[survival_df$transformed_RPL41 > quantile(survival_df$transformed_RPL41, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$RPL41_33_66 <- NA
survival_df$RPL41_33_66[survival_df$transformed_RPL41 <= quantile(survival_df$transformed_RPL41, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$RPL41_33_66[survival_df$transformed_RPL41 >= quantile(survival_df$transformed_RPL41, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - RPL41
for(strata in c("RPL41_median", "RPL41_25p", "RPL41_75p", "RPL41_66p", "RPL41_33p", "RPL41_33_66")){
  
  H_gTxt <- paste("HIGH RPL41", sep="")
  L_gTxt <- paste("LOW RPL41", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/upreg/RPL41/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/upreg/RPL41/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/upreg/RPL41/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}




# RPL41P1 
# Make survival plot
survival_df$RPL41P1_median <- NA
survival_df$RPL41P1_median[survival_df$transformed_RPL41P1 <= quantile(survival_df$transformed_RPL41P1, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$RPL41P1_median[survival_df$transformed_RPL41P1 > quantile(survival_df$transformed_RPL41P1, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$RPL41P1_25p <- NA
survival_df$RPL41P1_25p[survival_df$transformed_RPL41P1 <= quantile(survival_df$transformed_RPL41P1, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$RPL41P1_25p[survival_df$transformed_RPL41P1 > quantile(survival_df$transformed_RPL41P1, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$RPL41P1_75p <- NA
survival_df$RPL41P1_75p[survival_df$transformed_RPL41P1 <= quantile(survival_df$transformed_RPL41P1, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$RPL41P1_75p[survival_df$transformed_RPL41P1 > quantile(survival_df$transformed_RPL41P1, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$RPL41P1_66p <- NA
survival_df$RPL41P1_66p[survival_df$transformed_RPL41P1 <= quantile(survival_df$transformed_RPL41P1, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$RPL41P1_66p[survival_df$transformed_RPL41P1 > quantile(survival_df$transformed_RPL41P1, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$RPL41P1_33p <- NA
survival_df$RPL41P1_33p[survival_df$transformed_RPL41P1 <= quantile(survival_df$transformed_RPL41P1, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$RPL41P1_33p[survival_df$transformed_RPL41P1 > quantile(survival_df$transformed_RPL41P1, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$RPL41P1_33_66 <- NA
survival_df$RPL41P1_33_66[survival_df$transformed_RPL41P1 <= quantile(survival_df$transformed_RPL41P1, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$RPL41P1_33_66[survival_df$transformed_RPL41P1 >= quantile(survival_df$transformed_RPL41P1, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - RPL41P1
for(strata in c("RPL41P1_median", "RPL41P1_25p", "RPL41P1_75p", "RPL41P1_66p", "RPL41P1_33p", "RPL41P1_33_66")){
  
  H_gTxt <- paste("HIGH RPL41P1", sep="")
  L_gTxt <- paste("LOW RPL41P1", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/upreg/RPL41P1/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/upreg/RPL41P1/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/upreg/RPL41P1/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

# RPL41P5 
# Make survival plot
survival_df$RPL41P5_median <- NA
survival_df$RPL41P5_median[survival_df$transformed_RPL41P5 <= quantile(survival_df$transformed_RPL41P5, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$RPL41P5_median[survival_df$transformed_RPL41P5 > quantile(survival_df$transformed_RPL41P5, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$RPL41P5_25p <- NA
survival_df$RPL41P5_25p[survival_df$transformed_RPL41P5 <= quantile(survival_df$transformed_RPL41P5, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$RPL41P5_25p[survival_df$transformed_RPL41P5 > quantile(survival_df$transformed_RPL41P5, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$RPL41P5_75p <- NA
survival_df$RPL41P5_75p[survival_df$transformed_RPL41P5 <= quantile(survival_df$transformed_RPL41P5, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$RPL41P5_75p[survival_df$transformed_RPL41P5 > quantile(survival_df$transformed_RPL41P5, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$RPL41P5_66p <- NA
survival_df$RPL41P5_66p[survival_df$transformed_RPL41P5 <= quantile(survival_df$transformed_RPL41P5, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$RPL41P5_66p[survival_df$transformed_RPL41P5 > quantile(survival_df$transformed_RPL41P5, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$RPL41P5_33p <- NA
survival_df$RPL41P5_33p[survival_df$transformed_RPL41P5 <= quantile(survival_df$transformed_RPL41P5, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$RPL41P5_33p[survival_df$transformed_RPL41P5 > quantile(survival_df$transformed_RPL41P5, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$RPL41P5_33_66 <- NA
survival_df$RPL41P5_33_66[survival_df$transformed_RPL41P5 <= quantile(survival_df$transformed_RPL41P5, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$RPL41P5_33_66[survival_df$transformed_RPL41P5 >= quantile(survival_df$transformed_RPL41P5, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - RPL41P5
for(strata in c("RPL41P5_median", "RPL41P5_25p", "RPL41P5_75p", "RPL41P5_66p", "RPL41P5_33p", "RPL41P5_33_66")){
  
  H_gTxt <- paste("HIGH RPL41P5", sep="")
  L_gTxt <- paste("LOW RPL41P5", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/upreg/RPL41P5/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/upreg/RPL41P5/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/upreg/RPL41P5/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}



# CA9 
# Make survival plot
survival_df$CA9_median <- NA
survival_df$CA9_median[survival_df$transformed_CA9 <= quantile(survival_df$transformed_CA9, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$CA9_median[survival_df$transformed_CA9 > quantile(survival_df$transformed_CA9, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$CA9_25p <- NA
survival_df$CA9_25p[survival_df$transformed_CA9 <= quantile(survival_df$transformed_CA9, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$CA9_25p[survival_df$transformed_CA9 > quantile(survival_df$transformed_CA9, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$CA9_75p <- NA
survival_df$CA9_75p[survival_df$transformed_CA9 <= quantile(survival_df$transformed_CA9, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$CA9_75p[survival_df$transformed_CA9 > quantile(survival_df$transformed_CA9, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$CA9_66p <- NA
survival_df$CA9_66p[survival_df$transformed_CA9 <= quantile(survival_df$transformed_CA9, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$CA9_66p[survival_df$transformed_CA9 > quantile(survival_df$transformed_CA9, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$CA9_33p <- NA
survival_df$CA9_33p[survival_df$transformed_CA9 <= quantile(survival_df$transformed_CA9, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CA9_33p[survival_df$transformed_CA9 > quantile(survival_df$transformed_CA9, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$CA9_33_66 <- NA
survival_df$CA9_33_66[survival_df$transformed_CA9 <= quantile(survival_df$transformed_CA9, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CA9_33_66[survival_df$transformed_CA9 >= quantile(survival_df$transformed_CA9, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - CA9
for(strata in c("CA9_median", "CA9_25p", "CA9_75p", "CA9_66p", "CA9_33p", "CA9_33_66")){
  
  H_gTxt <- paste("HIGH CA9", sep="")
  L_gTxt <- paste("LOW CA9", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/upreg/CA9/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/upreg/CA9/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/upreg/CA9/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}



# Cox regression
meta <- c("Row.names", "X_PATIENT")

# Median
strata <- c("transformed_UBE2C", "UBE2C_median", "transformed_ESM1", "ESM1_median", "transformed_GPC3", "GPC3_median", "transformed_COL15A1", "COL15A1_median", "transformed_CDC20", "CDC20_median", "transformed_BIRC5", "BIRC5_median",
            "transformed_E2F1", "E2F1_median", "transformed_PTTG1", "PTTG1_median", "transformed_MDK", "MDK_median", "transformed_TOP2A", "TOP2A_median"
            )
# OS
os <- c("OS", "OS_months")
os_median_df <- survival_df %>% 
  select(all_of(c(meta, os, strata)))
os_median_df <- subset(os_median_df,  !is.na(OS) & !is.na(OS_months))
cox <- coxph(Surv(OS_months, OS) ~ UBE2C_median + ESM1_median + GPC3_median + COL15A1_median + CDC20_median + BIRC5_median + E2F1_median + PTTG1_median + MDK_median + TOP2A_median, data= os_median_df)
summary(cox)
# Remove MDK
cox <- coxph(Surv(OS_months, OS) ~ UBE2C_median + ESM1_median + GPC3_median + COL15A1_median + CDC20_median + BIRC5_median + E2F1_median + PTTG1_median + TOP2A_median, data= os_median_df)
summary(cox)
# Remove BIRC5
cox <- coxph(Surv(OS_months, OS) ~ UBE2C_median + ESM1_median + GPC3_median + COL15A1_median + CDC20_median + E2F1_median + PTTG1_median + TOP2A_median, data= os_median_df)
summary(cox)
# Remove E2F1
cox <- coxph(Surv(OS_months, OS) ~ UBE2C_median + ESM1_median + GPC3_median + COL15A1_median + CDC20_median + PTTG1_median + TOP2A_median, data= os_median_df)
summary(cox)
# Remove GPC3
cox <- coxph(Surv(OS_months, OS) ~ UBE2C_median + ESM1_median + COL15A1_median + CDC20_median + PTTG1_median + TOP2A_median, data= os_median_df)
summary(cox)
# Remove ESM1
cox <- coxph(Surv(OS_months, OS) ~ UBE2C_median + COL15A1_median + CDC20_median + PTTG1_median + TOP2A_median, data= os_median_df)
summary(cox)
# Remove CDC20
cox <- coxph(Surv(OS_months, OS) ~ UBE2C_median + COL15A1_median + PTTG1_median + TOP2A_median, data= os_median_df)
summary(cox)
# Remove UBE2C
cox <- coxph(Surv(OS_months, OS) ~ COL15A1_median + PTTG1_median + TOP2A_median, data= os_median_df)
summary(cox)
# Remove TOP2A
cox <- coxph(Surv(OS_months, OS) ~ COL15A1_median + PTTG1_median, data= os_median_df)
summary(cox)


# DSS
dss <- c("DSS", "DSS_months")
dss_median_df <- survival_df %>% 
  select(all_of(c(meta, dss, strata)))
dss_median_df <- subset(dss_median_df,  !is.na(DSS) & !is.na(DSS_months))
cox <- coxph(Surv(DSS_months, DSS) ~ UBE2C_median + ESM1_median + GPC3_median + COL15A1_median + CDC20_median + BIRC5_median + E2F1_median + PTTG1_median + MDK_median + TOP2A_median, data= dss_median_df)
summary(cox)
# Remove E2F1
cox <- coxph(Surv(DSS_months, DSS) ~ UBE2C_median + ESM1_median + GPC3_median + COL15A1_median + CDC20_median + BIRC5_median + PTTG1_median + MDK_median + TOP2A_median, data= dss_median_df)
summary(cox)
# Remove MDK
cox <- coxph(Surv(DSS_months, DSS) ~ UBE2C_median + ESM1_median + GPC3_median + COL15A1_median + CDC20_median + BIRC5_median + PTTG1_median + TOP2A_median, data= dss_median_df)
summary(cox)
# Remove PTTG1
cox <- coxph(Surv(DSS_months, DSS) ~ UBE2C_median + ESM1_median + GPC3_median + COL15A1_median + CDC20_median + BIRC5_median + TOP2A_median, data= dss_median_df)
summary(cox)
# Remove UBE2C
cox <- coxph(Surv(DSS_months, DSS) ~ ESM1_median + GPC3_median + COL15A1_median + CDC20_median + BIRC5_median + TOP2A_median, data= dss_median_df)
summary(cox)
# Remove BIRC5
cox <- coxph(Surv(DSS_months, DSS) ~ ESM1_median + GPC3_median + COL15A1_median + CDC20_median + TOP2A_median, data= dss_median_df)
summary(cox)
# Remove GPC3
cox <- coxph(Surv(DSS_months, DSS) ~ ESM1_median + COL15A1_median + CDC20_median + TOP2A_median, data= dss_median_df)
summary(cox)
# Remove ESM1
cox <- coxph(Surv(DSS_months, DSS) ~ COL15A1_median + CDC20_median + TOP2A_median, data= dss_median_df)
summary(cox)
# Remove CDC20
cox <- coxph(Surv(DSS_months, DSS) ~ COL15A1_median + TOP2A_median, data= dss_median_df)
summary(cox)



# PFI
pfi <- c("PFI", "PFI_months")
pfi_median_df <- survival_df %>% 
  select(all_of(c(meta, pfi, strata)))
pfi_median_df <- subset(pfi_median_df,  !is.na(PFI) & !is.na(PFI_months))
cox <- coxph(Surv(PFI_months, PFI) ~ UBE2C_median + ESM1_median + GPC3_median + COL15A1_median + CDC20_median + BIRC5_median + E2F1_median + PTTG1_median + MDK_median + TOP2A_median, data= pfi_median_df)
summary(cox)
# Remove UBE2C
cox <- coxph(Surv(PFI_months, PFI) ~ ESM1_median + GPC3_median + COL15A1_median + CDC20_median + BIRC5_median + E2F1_median + PTTG1_median + MDK_median + TOP2A_median, data= pfi_median_df)
summary(cox)
# Remove GPC3
cox <- coxph(Surv(PFI_months, PFI) ~ ESM1_median + COL15A1_median + CDC20_median + BIRC5_median + E2F1_median + PTTG1_median + MDK_median + TOP2A_median, data= pfi_median_df)
summary(cox)
# Remove MDK
cox <- coxph(Surv(PFI_months, PFI) ~ ESM1_median + COL15A1_median + CDC20_median + BIRC5_median + E2F1_median + PTTG1_median + TOP2A_median, data= pfi_median_df)
summary(cox)
# Remove CDC20
cox <- coxph(Surv(PFI_months, PFI) ~ ESM1_median + COL15A1_median + BIRC5_median + E2F1_median + PTTG1_median + TOP2A_median, data= pfi_median_df)
summary(cox)
# Remove BIRC5
cox <- coxph(Surv(PFI_months, PFI) ~ ESM1_median + COL15A1_median + E2F1_median + PTTG1_median + TOP2A_median, data= pfi_median_df)
summary(cox)
# Remove E2F1
cox <- coxph(Surv(PFI_months, PFI) ~ ESM1_median + COL15A1_median + PTTG1_median + TOP2A_median, data= pfi_median_df)
summary(cox)
# Remove ESM1
cox <- coxph(Surv(PFI_months, PFI) ~ COL15A1_median + PTTG1_median + TOP2A_median, data= pfi_median_df)
summary(cox)
# Remove COL15A1
cox <- coxph(Surv(PFI_months, PFI) ~ PTTG1_median + TOP2A_median, data= pfi_median_df)
summary(cox)
# Remove TOP2A
cox <- coxph(Surv(PFI_months, PFI) ~ PTTG1_median, data= pfi_median_df)
summary(cox)
