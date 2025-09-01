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
# These genes are protein coding and downregulated from the paired TCGA analysis (TCGA-BC-A10X pair removed) and healthy_vs_tcga
downreg_genes <- c("ENSG00000205358.4",
                   "ENSG00000125144.14",
                   "ENSG00000019169.10",
                   "ENSG00000160339.16",
                   "ENSG00000198417.7",
                   "ENSG00000205364.4",
                   "ENSG00000105697.9",
                   "ENSG00000140505.7",
                   "ENSG00000263761.3",
                   "ENSG00000145708.11",
                   "ENSG00000282301.3",
                   "ENSG00000183054.11",
                   "ENSG00000214967.5"
                   )


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
names(survival_df)[names(survival_df) == "ENSG00000205358.4"] <- "MT1H"
names(survival_df)[names(survival_df) == "ENSG00000125144.14"] <- "MT1G"
names(survival_df)[names(survival_df) == "ENSG00000019169.10"] <- "MARCO"
names(survival_df)[names(survival_df) == "ENSG00000160339.16"] <- "FCN2"
names(survival_df)[names(survival_df) == "ENSG00000198417.7"] <- "MT1F"
names(survival_df)[names(survival_df) == "ENSG00000205364.4"] <- "MT1M"
names(survival_df)[names(survival_df) == "ENSG00000105697.9"] <- "HAMP"
names(survival_df)[names(survival_df) == "ENSG00000140505.7"] <- "CYP1A2"
names(survival_df)[names(survival_df) == "ENSG00000263761.3"] <- "GDF2"
names(survival_df)[names(survival_df) == "ENSG00000145708.11"] <- "CRHBP"
names(survival_df)[names(survival_df) == "ENSG00000282301.3"] <- "CYP3A7.CYP3A51P"

# From TCGA vs GTEX
names(survival_df)[names(survival_df) == "ENSG00000183054.11"] <- "RGPD6"
names(survival_df)[names(survival_df) == "ENSG00000214967.5"] <- "NPIPA7"



# Apply function to new transformed gene count column
transform_tpm_counts <- function(count){
  transformed_count <- (2^count) - 0.001
  return(transformed_count)
}


survival_df["transformed_MT1H"] <- lapply(survival_df["MT1H"], transform_tpm_counts)
survival_df["transformed_MT1G"] <- lapply(survival_df["MT1G"], transform_tpm_counts)
survival_df["transformed_MARCO"] <- lapply(survival_df["MARCO"], transform_tpm_counts)
survival_df["transformed_FCN2"] <- lapply(survival_df["FCN2"], transform_tpm_counts)
survival_df["transformed_MT1F"] <- lapply(survival_df["MT1F"], transform_tpm_counts)
survival_df["transformed_MT1M"] <- lapply(survival_df["MT1M"], transform_tpm_counts)
survival_df["transformed_HAMP"] <- lapply(survival_df["HAMP"], transform_tpm_counts)
survival_df["transformed_CYP1A2"] <- lapply(survival_df["CYP1A2"], transform_tpm_counts)
survival_df["transformed_GDF2"] <- lapply(survival_df["GDF2"], transform_tpm_counts)
survival_df["transformed_CRHBP"] <- lapply(survival_df["CRHBP"], transform_tpm_counts)
survival_df["transformed_CYP3A7.CYP3A51P"] <- lapply(survival_df["CYP3A7.CYP3A51P"], transform_tpm_counts)

# From TCGA vs GTEX
survival_df["transformed_RGPD6"] <- lapply(survival_df["RGPD6"], transform_tpm_counts)
survival_df["transformed_NPIPA7"] <- lapply(survival_df["NPIPA7"], transform_tpm_counts)


# 1 day = 0.0328767 days 
days_to_months_converter <- function(no_days){
  no_months <- round(no_days*0.0328767, 1)
  return(no_months)
}

survival_df["OS_months"] <- lapply(survival_df["OS.time"], days_to_months_converter)
survival_df["DSS_months"] <- lapply(survival_df["DSS.time"], days_to_months_converter)
survival_df["PFI_months"] <- lapply(survival_df["PFI.time"], days_to_months_converter)


# MT1G
# Make survival plot
survival_df$MT1G_median <- NA
survival_df$MT1G_median[survival_df$transformed_MT1G <= quantile(survival_df$transformed_MT1G, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$MT1G_median[survival_df$transformed_MT1G > quantile(survival_df$transformed_MT1G, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$MT1G_25p <- NA
survival_df$MT1G_25p[survival_df$transformed_MT1G <= quantile(survival_df$transformed_MT1G, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$MT1G_25p[survival_df$transformed_MT1G > quantile(survival_df$transformed_MT1G, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$MT1G_75p <- NA
survival_df$MT1G_75p[survival_df$transformed_MT1G <= quantile(survival_df$transformed_MT1G, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$MT1G_75p[survival_df$transformed_MT1G > quantile(survival_df$transformed_MT1G, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$MT1G_66p <- NA
survival_df$MT1G_66p[survival_df$transformed_MT1G <= quantile(survival_df$transformed_MT1G, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$MT1G_66p[survival_df$transformed_MT1G > quantile(survival_df$transformed_MT1G, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$MT1G_33p <- NA
survival_df$MT1G_33p[survival_df$transformed_MT1G <= quantile(survival_df$transformed_MT1G, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$MT1G_33p[survival_df$transformed_MT1G > quantile(survival_df$transformed_MT1G, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$MT1G_33_66 <- NA
survival_df$MT1G_33_66[survival_df$transformed_MT1G <= quantile(survival_df$transformed_MT1G, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$MT1G_33_66[survival_df$transformed_MT1G >= quantile(survival_df$transformed_MT1G, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - MT1G
for(strata in c("MT1G_median", "MT1G_25p", "MT1G_75p", "MT1G_66p", "MT1G_33p", "MT1G_33_66")){
  
  H_gTxt <- paste("HIGH MT1G", sep="")
  L_gTxt <- paste("LOW MT1G", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/downreg/MT1G/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/downreg/MT1G/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/downreg/MT1G/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

# MT1F
# Make survival plot
survival_df$MT1F_median <- NA
survival_df$MT1F_median[survival_df$transformed_MT1F <= quantile(survival_df$transformed_MT1F, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$MT1F_median[survival_df$transformed_MT1F > quantile(survival_df$transformed_MT1F, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$MT1F_25p <- NA
survival_df$MT1F_25p[survival_df$transformed_MT1F <= quantile(survival_df$transformed_MT1F, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$MT1F_25p[survival_df$transformed_MT1F > quantile(survival_df$transformed_MT1F, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$MT1F_75p <- NA
survival_df$MT1F_75p[survival_df$transformed_MT1F <= quantile(survival_df$transformed_MT1F, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$MT1F_75p[survival_df$transformed_MT1F > quantile(survival_df$transformed_MT1F, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$MT1F_66p <- NA
survival_df$MT1F_66p[survival_df$transformed_MT1F <= quantile(survival_df$transformed_MT1F, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$MT1F_66p[survival_df$transformed_MT1F > quantile(survival_df$transformed_MT1F, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$MT1F_33p <- NA
survival_df$MT1F_33p[survival_df$transformed_MT1F <= quantile(survival_df$transformed_MT1F, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$MT1F_33p[survival_df$transformed_MT1F > quantile(survival_df$transformed_MT1F, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$MT1F_33_66 <- NA
survival_df$MT1F_33_66[survival_df$transformed_MT1F <= quantile(survival_df$transformed_MT1F, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$MT1F_33_66[survival_df$transformed_MT1F >= quantile(survival_df$transformed_MT1F, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - MT1F
for(strata in c("MT1F_median", "MT1F_25p", "MT1F_75p", "MT1F_66p", "MT1F_33p", "MT1F_33_66")){
  
  H_gTxt <- paste("HIGH MT1F", sep="")
  L_gTxt <- paste("LOW MT1F", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/downreg/MT1F/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/downreg/MT1F/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/downreg/MT1F/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

# MT1M
# Make survival plot
survival_df$MT1M_median <- NA
survival_df$MT1M_median[survival_df$transformed_MT1M <= quantile(survival_df$transformed_MT1M, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$MT1M_median[survival_df$transformed_MT1M > quantile(survival_df$transformed_MT1M, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$MT1M_25p <- NA
survival_df$MT1M_25p[survival_df$transformed_MT1M <= quantile(survival_df$transformed_MT1M, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$MT1M_25p[survival_df$transformed_MT1M > quantile(survival_df$transformed_MT1M, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$MT1M_75p <- NA
survival_df$MT1M_75p[survival_df$transformed_MT1M <= quantile(survival_df$transformed_MT1M, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$MT1M_75p[survival_df$transformed_MT1M > quantile(survival_df$transformed_MT1M, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$MT1M_66p <- NA
survival_df$MT1M_66p[survival_df$transformed_MT1M <= quantile(survival_df$transformed_MT1M, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$MT1M_66p[survival_df$transformed_MT1M > quantile(survival_df$transformed_MT1M, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$MT1M_33p <- NA
survival_df$MT1M_33p[survival_df$transformed_MT1M <= quantile(survival_df$transformed_MT1M, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$MT1M_33p[survival_df$transformed_MT1M > quantile(survival_df$transformed_MT1M, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$MT1M_33_66 <- NA
survival_df$MT1M_33_66[survival_df$transformed_MT1M <= quantile(survival_df$transformed_MT1M, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$MT1M_33_66[survival_df$transformed_MT1M >= quantile(survival_df$transformed_MT1M, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - MT1M
for(strata in c("MT1M_median", "MT1M_25p", "MT1M_75p", "MT1M_66p", "MT1M_33p", "MT1M_33_66")){
  
  H_gTxt <- paste("HIGH MT1M", sep="")
  L_gTxt <- paste("LOW MT1M", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/downreg/MT1M/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/downreg/MT1M/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/downreg/MT1M/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

# HAMP
# Make survival plot
survival_df$HAMP_median <- NA
survival_df$HAMP_median[survival_df$transformed_HAMP <= quantile(survival_df$transformed_HAMP, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$HAMP_median[survival_df$transformed_HAMP > quantile(survival_df$transformed_HAMP, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$HAMP_25p <- NA
survival_df$HAMP_25p[survival_df$transformed_HAMP <= quantile(survival_df$transformed_HAMP, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$HAMP_25p[survival_df$transformed_HAMP > quantile(survival_df$transformed_HAMP, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$HAMP_75p <- NA
survival_df$HAMP_75p[survival_df$transformed_HAMP <= quantile(survival_df$transformed_HAMP, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$HAMP_75p[survival_df$transformed_HAMP > quantile(survival_df$transformed_HAMP, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$HAMP_66p <- NA
survival_df$HAMP_66p[survival_df$transformed_HAMP <= quantile(survival_df$transformed_HAMP, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$HAMP_66p[survival_df$transformed_HAMP > quantile(survival_df$transformed_HAMP, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$HAMP_33p <- NA
survival_df$HAMP_33p[survival_df$transformed_HAMP <= quantile(survival_df$transformed_HAMP, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$HAMP_33p[survival_df$transformed_HAMP > quantile(survival_df$transformed_HAMP, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$HAMP_33_66 <- NA
survival_df$HAMP_33_66[survival_df$transformed_HAMP <= quantile(survival_df$transformed_HAMP, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$HAMP_33_66[survival_df$transformed_HAMP >= quantile(survival_df$transformed_HAMP, 0.66, na.rm=TRUE)] <- "HIGH"

# SURVIVAL CURVES - HAMP
for(strata in c("HAMP_median", "HAMP_25p", "HAMP_75p", "HAMP_66p", "HAMP_33p", "HAMP_33_66")){
  
  H_gTxt <- paste("HIGH HAMP", sep="")
  L_gTxt <- paste("LOW HAMP", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/downreg/HAMP/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/downreg/HAMP/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/downreg/HAMP/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}


# CYP1A2
# Make survival plot
survival_df$CYP1A2_median <- NA
survival_df$CYP1A2_median[survival_df$transformed_CYP1A2 <= quantile(survival_df$transformed_CYP1A2, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$CYP1A2_median[survival_df$transformed_CYP1A2 > quantile(survival_df$transformed_CYP1A2, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$CYP1A2_25p <- NA
survival_df$CYP1A2_25p[survival_df$transformed_CYP1A2 <= quantile(survival_df$transformed_CYP1A2, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$CYP1A2_25p[survival_df$transformed_CYP1A2 > quantile(survival_df$transformed_CYP1A2, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$CYP1A2_75p <- NA
survival_df$CYP1A2_75p[survival_df$transformed_CYP1A2 <= quantile(survival_df$transformed_CYP1A2, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$CYP1A2_75p[survival_df$transformed_CYP1A2 > quantile(survival_df$transformed_CYP1A2, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$CYP1A2_66p <- NA
survival_df$CYP1A2_66p[survival_df$transformed_CYP1A2 <= quantile(survival_df$transformed_CYP1A2, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$CYP1A2_66p[survival_df$transformed_CYP1A2 > quantile(survival_df$transformed_CYP1A2, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$CYP1A2_33p <- NA
survival_df$CYP1A2_33p[survival_df$transformed_CYP1A2 <= quantile(survival_df$transformed_CYP1A2, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CYP1A2_33p[survival_df$transformed_CYP1A2 > quantile(survival_df$transformed_CYP1A2, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$CYP1A2_33_66 <- NA
survival_df$CYP1A2_33_66[survival_df$transformed_CYP1A2 <= quantile(survival_df$transformed_CYP1A2, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CYP1A2_33_66[survival_df$transformed_CYP1A2 >= quantile(survival_df$transformed_CYP1A2, 0.66, na.rm=TRUE)] <- "HIGH"

# SURVIVAL CURVES - CYP1A2
for(strata in c("CYP1A2_median", "CYP1A2_25p", "CYP1A2_75p", "CYP1A2_66p", "CYP1A2_33p", "CYP1A2_33_66")){
  
  H_gTxt <- paste("HIGH CYP1A2", sep="")
  L_gTxt <- paste("LOW CYP1A2", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/downreg/CYP1A2/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/downreg/CYP1A2/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/downreg/CYP1A2/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}


# CRHBP
# Make survival plot
survival_df$CRHBP_median <- NA
survival_df$CRHBP_median[survival_df$transformed_CRHBP <= quantile(survival_df$transformed_CRHBP, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$CRHBP_median[survival_df$transformed_CRHBP > quantile(survival_df$transformed_CRHBP, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$CRHBP_25p <- NA
survival_df$CRHBP_25p[survival_df$transformed_CRHBP <= quantile(survival_df$transformed_CRHBP, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$CRHBP_25p[survival_df$transformed_CRHBP > quantile(survival_df$transformed_CRHBP, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$CRHBP_75p <- NA
survival_df$CRHBP_75p[survival_df$transformed_CRHBP <= quantile(survival_df$transformed_CRHBP, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$CRHBP_75p[survival_df$transformed_CRHBP > quantile(survival_df$transformed_CRHBP, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$CRHBP_66p <- NA
survival_df$CRHBP_66p[survival_df$transformed_CRHBP <= quantile(survival_df$transformed_CRHBP, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$CRHBP_66p[survival_df$transformed_CRHBP > quantile(survival_df$transformed_CRHBP, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$CRHBP_33p <- NA
survival_df$CRHBP_33p[survival_df$transformed_CRHBP <= quantile(survival_df$transformed_CRHBP, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CRHBP_33p[survival_df$transformed_CRHBP > quantile(survival_df$transformed_CRHBP, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$CRHBP_33_66 <- NA
survival_df$CRHBP_33_66[survival_df$transformed_CRHBP <= quantile(survival_df$transformed_CRHBP, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CRHBP_33_66[survival_df$transformed_CRHBP >= quantile(survival_df$transformed_CRHBP, 0.66, na.rm=TRUE)] <- "HIGH"

# SURVIVAL CURVES - CRHBP
for(strata in c("CRHBP_median", "CRHBP_25p", "CRHBP_75p", "CRHBP_66p", "CRHBP_33p", "CRHBP_33_66")){
  
  H_gTxt <- paste("HIGH CRHBP", sep="")
  L_gTxt <- paste("LOW CRHBP", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/downreg/CRHBP/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
                  pval = TRUE, pval.coord = c(70, 85), pval.size = 7,
                  size = 1, censor.size = 2)
  
  tiff(filename = paste ("survival_plots/removed_sample/downreg/CRHBP/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/downreg/CRHBP/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}


# CYP3A7.CYP3A51P
# Make survival plot
survival_df$CYP3A7.CYP3A51P_median <- NA
survival_df$CYP3A7.CYP3A51P_median[survival_df$transformed_CYP3A7.CYP3A51P <= quantile(survival_df$transformed_CYP3A7.CYP3A51P, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$CYP3A7.CYP3A51P_median[survival_df$transformed_CYP3A7.CYP3A51P > quantile(survival_df$transformed_CYP3A7.CYP3A51P, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$CYP3A7.CYP3A51P_25p <- NA
survival_df$CYP3A7.CYP3A51P_25p[survival_df$transformed_CYP3A7.CYP3A51P <= quantile(survival_df$transformed_CYP3A7.CYP3A51P, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$CYP3A7.CYP3A51P_25p[survival_df$transformed_CYP3A7.CYP3A51P > quantile(survival_df$transformed_CYP3A7.CYP3A51P, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$CYP3A7.CYP3A51P_75p <- NA
survival_df$CYP3A7.CYP3A51P_75p[survival_df$transformed_CYP3A7.CYP3A51P <= quantile(survival_df$transformed_CYP3A7.CYP3A51P, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$CYP3A7.CYP3A51P_75p[survival_df$transformed_CYP3A7.CYP3A51P > quantile(survival_df$transformed_CYP3A7.CYP3A51P, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$CYP3A7.CYP3A51P_66p <- NA
survival_df$CYP3A7.CYP3A51P_66p[survival_df$transformed_CYP3A7.CYP3A51P <= quantile(survival_df$transformed_CYP3A7.CYP3A51P, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$CYP3A7.CYP3A51P_66p[survival_df$transformed_CYP3A7.CYP3A51P > quantile(survival_df$transformed_CYP3A7.CYP3A51P, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$CYP3A7.CYP3A51P_33p <- NA
survival_df$CYP3A7.CYP3A51P_33p[survival_df$transformed_CYP3A7.CYP3A51P <= quantile(survival_df$transformed_CYP3A7.CYP3A51P, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CYP3A7.CYP3A51P_33p[survival_df$transformed_CYP3A7.CYP3A51P > quantile(survival_df$transformed_CYP3A7.CYP3A51P, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$CYP3A7.CYP3A51P_33_66 <- NA
survival_df$CYP3A7.CYP3A51P_33_66[survival_df$transformed_CYP3A7.CYP3A51P <= quantile(survival_df$transformed_CYP3A7.CYP3A51P, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CYP3A7.CYP3A51P_33_66[survival_df$transformed_CYP3A7.CYP3A51P >= quantile(survival_df$transformed_CYP3A7.CYP3A51P, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - CYP3A7.CYP3A51P
for(strata in c("CYP3A7.CYP3A51P_median", "CYP3A7.CYP3A51P_25p", "CYP3A7.CYP3A51P_75p", "CYP3A7.CYP3A51P_66p", "CYP3A7.CYP3A51P_33p", "CYP3A7.CYP3A51P_33_66")){
  
  H_gTxt <- paste("HIGH CYP3A7.CYP3A51P", sep="")
  L_gTxt <- paste("LOW CYP3A7.CYP3A51P", sep="")
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
  
  tiff(filename = paste("survival_plots/removed_sample/downreg/CYP3A7.CYP3A51P/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/removed_sample/downreg/CYP3A7.CYP3A51P/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/removed_sample/downreg/CYP3A7.CYP3A51P/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}



# Plot already made before
# GDF2
# Make survival plot
survival_df$GDF2_median <- NA
survival_df$GDF2_median[survival_df$transformed_GDF2 <= quantile(survival_df$transformed_GDF2, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$GDF2_median[survival_df$transformed_GDF2 > quantile(survival_df$transformed_GDF2, 0.5, na.rm=TRUE)] <- "HIGH"


# MT1H
# Make survival plot
survival_df$MT1H_median <- NA
survival_df$MT1H_median[survival_df$transformed_MT1H <= quantile(survival_df$transformed_MT1H, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$MT1H_median[survival_df$transformed_MT1H > quantile(survival_df$transformed_MT1H, 0.5, na.rm=TRUE)] <- "HIGH"


# MARCO
# Make survival plot
survival_df$MARCO_median <- NA
survival_df$MARCO_median[survival_df$transformed_MARCO <= quantile(survival_df$transformed_MARCO, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$MARCO_median[survival_df$transformed_MARCO > quantile(survival_df$transformed_MARCO, 0.5, na.rm=TRUE)] <- "HIGH"


# FCN2
# Make survival plot
survival_df$FCN2_median <- NA
survival_df$FCN2_median[survival_df$transformed_FCN2 <= quantile(survival_df$transformed_FCN2, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$FCN2_median[survival_df$transformed_FCN2 > quantile(survival_df$transformed_FCN2, 0.5, na.rm=TRUE)] <- "HIGH"


# From TCGA vs GTEX (survival plots already made in other survival analysis script)
# RGPD6
# Make survival plot
survival_df$RGPD6_median <- NA
survival_df$RGPD6_median[survival_df$transformed_RGPD6 <= quantile(survival_df$transformed_RGPD6, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$RGPD6_median[survival_df$transformed_RGPD6 > quantile(survival_df$transformed_RGPD6, 0.5, na.rm=TRUE)] <- "HIGH"

# NPIPA7
# Make survival plot
survival_df$NPIPA7_median <- NA
survival_df$NPIPA7_median[survival_df$transformed_NPIPA7 <= quantile(survival_df$transformed_NPIPA7, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$NPIPA7_median[survival_df$transformed_NPIPA7 > quantile(survival_df$transformed_NPIPA7, 0.5, na.rm=TRUE)] <- "HIGH"


# Cox regression
meta <- c("Row.names", "X_PATIENT")

# Median
strata <- c("transformed_MT1H", "MT1H_median", "transformed_MT1G", "MT1G_median", "transformed_MARCO", "MARCO_median", "transformed_FCN2", "FCN2_median", "transformed_MT1F", "MT1F_median",
            "transformed_MT1M", "MT1M_median", "transformed_HAMP", "HAMP_median", "transformed_CYP1A2", "CYP1A2_median", "transformed_CRHBP", "CRHBP_median", "transformed_GDF2", "GDF2_median",
            "transformed_RGPD6", "RGPD6_median", "transformed_NPIPA7", "NPIPA7_median", "transformed_CYP3A7.CYP3A51P", "CYP3A7.CYP3A51P_median"
)

# OS
os <- c("OS", "OS_months")
os_median_df <- survival_df %>% 
  select(all_of(c(meta, os, strata)))
os_median_df <- subset(os_median_df,  !is.na(OS) & !is.na(OS_months))
cox <- coxph(Surv(OS_months, OS) ~ MT1H_median + MT1G_median + MARCO_median + FCN2_median + MT1F_median + MT1M_median + HAMP_median + CYP1A2_median + GDF2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= os_median_df)
summary(cox)
# Remove GDF2
cox <- coxph(Surv(OS_months, OS) ~ MT1H_median + MT1G_median + MARCO_median + FCN2_median + MT1F_median + MT1M_median + HAMP_median + CYP1A2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= os_median_df)
summary(cox)
# Remove MT1M
cox <- coxph(Surv(OS_months, OS) ~ MT1H_median + MT1G_median + MARCO_median + FCN2_median + MT1F_median + HAMP_median + CYP1A2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= os_median_df)
summary(cox)
# Remove MARCO
cox <- coxph(Surv(OS_months, OS) ~ MT1H_median + MT1G_median + FCN2_median + MT1F_median + HAMP_median + CYP1A2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= os_median_df)
summary(cox)
# Remove HAMP
cox <- coxph(Surv(OS_months, OS) ~ MT1H_median + MT1G_median + FCN2_median + MT1F_median + CYP1A2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= os_median_df)
summary(cox)
# Remove MT1G
cox <- coxph(Surv(OS_months, OS) ~ MT1H_median + FCN2_median + MT1F_median + CYP1A2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= os_median_df)
summary(cox)
# MT1F
cox <- coxph(Surv(OS_months, OS) ~ MT1H_median + FCN2_median + CYP1A2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= os_median_df)
summary(cox)
# remove RGPD6
cox <- coxph(Surv(OS_months, OS) ~ MT1H_median + FCN2_median + CYP1A2_median + CRHBP_median + NPIPA7_median, data= os_median_df)
summary(cox)
# Remove MT1H
cox <- coxph(Surv(OS_months, OS) ~ FCN2_median + CYP1A2_median + CRHBP_median + NPIPA7_median, data= os_median_df)
summary(cox)
# Remove FCN2
cox <- coxph(Surv(OS_months, OS) ~ CYP1A2_median + CRHBP_median + NPIPA7_median, data= os_median_df)
summary(cox)
# Remove NPIPA7
cox <- coxph(Surv(OS_months, OS) ~ CYP1A2_median + CRHBP_median, data= os_median_df)
summary(cox)
# Remove CYP1A2
cox <- coxph(Surv(OS_months, OS) ~ CRHBP_median, data= os_median_df)
summary(cox)


# DSS
dss <- c("DSS", "DSS_months")
dss_median_df <- survival_df %>% 
  select(all_of(c(meta, dss, strata)))
dss_median_df <- subset(dss_median_df,  !is.na(DSS) & !is.na(DSS_months))
cox <- coxph(Surv(DSS_months, DSS) ~ MT1H_median + MT1G_median + MARCO_median + FCN2_median + MT1F_median + MT1M_median + HAMP_median + CYP1A2_median + GDF2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= dss_median_df)
summary(cox)
# Remove MT1H
cox <- coxph(Surv(DSS_months, DSS) ~ MT1G_median + MARCO_median + FCN2_median + MT1F_median + MT1M_median + HAMP_median + CYP1A2_median + GDF2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= dss_median_df)
summary(cox)
# Remove FCN2
cox <- coxph(Surv(DSS_months, DSS) ~ MT1G_median + MARCO_median + MT1F_median + MT1M_median + HAMP_median + CYP1A2_median + GDF2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= dss_median_df)
summary(cox)
# Remove MT1G
cox <- coxph(Surv(DSS_months, DSS) ~ MARCO_median + MT1F_median + MT1M_median + HAMP_median + CYP1A2_median + GDF2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= dss_median_df)
summary(cox)
# Remove RGPD6
cox <- coxph(Surv(DSS_months, DSS) ~ MARCO_median + MT1F_median + MT1M_median + HAMP_median + CYP1A2_median + GDF2_median + CRHBP_median + NPIPA7_median, data= dss_median_df)
summary(cox)
# Remove NPIPA7
cox <- coxph(Surv(DSS_months, DSS) ~ MARCO_median + MT1F_median + MT1M_median + HAMP_median + CYP1A2_median + GDF2_median + CRHBP_median, data= dss_median_df)
summary(cox)
# Remove MARCO
cox <- coxph(Surv(DSS_months, DSS) ~ MT1F_median + MT1M_median + HAMP_median + CYP1A2_median + GDF2_median + CRHBP_median, data= dss_median_df)
summary(cox)
# Remove MT1F
cox <- coxph(Surv(DSS_months, DSS) ~ MT1M_median + HAMP_median + CYP1A2_median + GDF2_median + CRHBP_median, data= dss_median_df)
summary(cox)
# Remove MT1M
cox <- coxph(Surv(DSS_months, DSS) ~ HAMP_median + CYP1A2_median + GDF2_median + CRHBP_median, data= dss_median_df)
summary(cox)
# Remove GDF2
cox <- coxph(Surv(DSS_months, DSS) ~ HAMP_median + CYP1A2_median + CRHBP_median, data= dss_median_df)
summary(cox)
# Remove HAMP
cox <- coxph(Surv(DSS_months, DSS) ~ CYP1A2_median + CRHBP_median, data= dss_median_df)
summary(cox)
# Remove CYP1A2
cox <- coxph(Surv(DSS_months, DSS) ~ CRHBP_median, data= dss_median_df)
summary(cox)



# PFI
pfi <- c("PFI", "PFI_months")
pfi_median_df <- survival_df %>% 
  select(all_of(c(meta, pfi, strata)))
pfi_median_df <- subset(pfi_median_df,  !is.na(PFI) & !is.na(PFI_months))
cox <- coxph(Surv(PFI_months, PFI) ~ MT1H_median + MT1G_median + MARCO_median + FCN2_median + MT1F_median + MT1M_median + HAMP_median + CYP1A2_median + GDF2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= pfi_median_df)
summary(cox)
# Remove GDF2
cox <- coxph(Surv(PFI_months, PFI) ~ MT1H_median + MT1G_median + MARCO_median + FCN2_median + MT1F_median + MT1M_median + HAMP_median + CYP1A2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= pfi_median_df)
summary(cox)
# Remove MT1M
cox <- coxph(Surv(PFI_months, PFI) ~ MT1H_median + MT1G_median + MARCO_median + FCN2_median + MT1F_median + HAMP_median + CYP1A2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= pfi_median_df)
summary(cox)
# Remove FCN2
cox <- coxph(Surv(PFI_months, PFI) ~ MT1H_median + MT1G_median + MARCO_median + MT1F_median + HAMP_median + CYP1A2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= pfi_median_df)
summary(cox)
# Remove MARCO
cox <- coxph(Surv(PFI_months, PFI) ~ MT1H_median + MT1G_median + MT1F_median + HAMP_median + CYP1A2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= pfi_median_df)
summary(cox)
# Remove MT1G
cox <- coxph(Surv(PFI_months, PFI) ~ MT1H_median + MT1F_median + HAMP_median + CYP1A2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= pfi_median_df)
summary(cox)
# Remove HAMP
cox <- coxph(Surv(PFI_months, PFI) ~ MT1H_median + MT1F_median + CYP1A2_median + CRHBP_median + RGPD6_median + NPIPA7_median, data= pfi_median_df)
summary(cox)
# Remove RGPD6
cox <- coxph(Surv(PFI_months, PFI) ~ MT1H_median + MT1F_median + CYP1A2_median + CRHBP_median + NPIPA7_median, data= pfi_median_df)
summary(cox)
# Remove MT1H
cox <- coxph(Surv(PFI_months, PFI) ~ MT1F_median + CYP1A2_median + CRHBP_median + NPIPA7_median, data= pfi_median_df)
summary(cox)
# Remove MT1F
cox <- coxph(Surv(PFI_months, PFI) ~ CYP1A2_median + CRHBP_median + NPIPA7_median, data= pfi_median_df)
summary(cox)
# Remove CYP1A2
cox <- coxph(Surv(PFI_months, PFI) ~ CRHBP_median + NPIPA7_median, data= pfi_median_df)
summary(cox)
