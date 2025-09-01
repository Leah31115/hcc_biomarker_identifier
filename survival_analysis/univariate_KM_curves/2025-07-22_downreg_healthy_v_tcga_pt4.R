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
downreg_genes <- c("ENSG00000011052.21",
                   "ENSG00000136305.11",
                   "ENSG00000142694.6",
                   "ENSG00000183054.11",
                   "ENSG00000188223.9",
                   "ENSG00000205923.3",
                   "ENSG00000211675.2",
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
names(survival_df)[names(survival_df) == "ENSG00000011052.21"] <- "NME1_NME2"
names(survival_df)[names(survival_df) == "ENSG00000136305.11"] <- "CIDEB"
names(survival_df)[names(survival_df) == "ENSG00000142694.6"] <- "EVA1B"
names(survival_df)[names(survival_df) == "ENSG00000183054.11"] <- "RGPD6"
names(survival_df)[names(survival_df) == "ENSG00000188223.9"] <- "AD000671.1"
names(survival_df)[names(survival_df) == "ENSG00000205923.3"] <- "CEMP1"
names(survival_df)[names(survival_df) == "ENSG00000211675.2"] <- "IGLC1"
names(survival_df)[names(survival_df) == "ENSG00000214967.5"] <- "NPIPA7"


# Antilog gene count function
antilog_counts <- function(tpm){
  antilog_count <- (2^tpm) - 0.001
  return(antilog_count)
}

survival_df["antilog_NME1_NME2"] <- lapply(survival_df["NME1_NME2"], antilog_counts)
survival_df["antilog_CIDEB"] <- lapply(survival_df["CIDEB"], antilog_counts)
survival_df["antilog_EVA1B"] <- lapply(survival_df["EVA1B"], antilog_counts)
survival_df["antilog_RGPD6"] <- lapply(survival_df["RGPD6"], antilog_counts)
survival_df["antilog_AD000671.1"] <- lapply(survival_df["AD000671.1"], antilog_counts)
survival_df["antilog_CEMP1"] <- lapply(survival_df["CEMP1"], antilog_counts)
survival_df["antilog_IGLC1"] <- lapply(survival_df["IGLC1"], antilog_counts)
survival_df["antilog_NPIPA7"] <- lapply(survival_df["NPIPA7"], antilog_counts)

# 1 day = 0.0328767 days 
days_to_months_converter <- function(no_days){
  no_months <- round(no_days*0.0328767, 1)
  return(no_months)
}

survival_df["OS_months"] <- lapply(survival_df["OS.time"], days_to_months_converter)
survival_df["DSS_months"] <- lapply(survival_df["DSS.time"], days_to_months_converter)
survival_df["PFI_months"] <- lapply(survival_df["PFI.time"], days_to_months_converter)

# NO COMPARISON - ALL 0 counts
# NME1_NME2
# Make survival plot
survival_df$NME1_NME2_median <- NA
survival_df$NME1_NME2_median[survival_df$antilog_NME1_NME2 <= quantile(survival_df$antilog_NME1_NME2, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$NME1_NME2_median[survival_df$antilog_NME1_NME2 > quantile(survival_df$antilog_NME1_NME2, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$NME1_NME2_25p <- NA
survival_df$NME1_NME2_25p[survival_df$antilog_NME1_NME2 <= quantile(survival_df$antilog_NME1_NME2, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$NME1_NME2_25p[survival_df$antilog_NME1_NME2 > quantile(survival_df$antilog_NME1_NME2, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$NME1_NME2_75p <- NA
survival_df$NME1_NME2_75p[survival_df$antilog_NME1_NME2 <= quantile(survival_df$antilog_NME1_NME2, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$NME1_NME2_75p[survival_df$antilog_NME1_NME2 > quantile(survival_df$antilog_NME1_NME2, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$NME1_NME2_66p <- NA
survival_df$NME1_NME2_66p[survival_df$antilog_NME1_NME2 <= quantile(survival_df$antilog_NME1_NME2, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$NME1_NME2_66p[survival_df$antilog_NME1_NME2 > quantile(survival_df$antilog_NME1_NME2, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$NME1_NME2_33p <- NA
survival_df$NME1_NME2_33p[survival_df$antilog_NME1_NME2 <= quantile(survival_df$antilog_NME1_NME2, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$NME1_NME2_33p[survival_df$antilog_NME1_NME2 > quantile(survival_df$antilog_NME1_NME2, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$NME1_NME2_33_66 <- NA
survival_df$NME1_NME2_33_66[survival_df$antilog_NME1_NME2 <= quantile(survival_df$antilog_NME1_NME2, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$NME1_NME2_33_66[survival_df$antilog_NME1_NME2 >= quantile(survival_df$antilog_NME1_NME2, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - NME1_NME2
for(strata in c("NME1_NME2_median", "NME1_NME2_25p", "NME1_NME2_75p", "NME1_NME2_66p", "NME1_NME2_33p", "NME1_NME2_33_66")){
  
  H_gTxt <- paste("HIGH NME1_NME2", sep="")
  L_gTxt <- paste("LOW NME1_NME2", sep="")
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
  
  tiff(filename = paste("survival_plots/downreg/healthy_vs_tcga/NME1-NME2/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/downreg/healthy_vs_tcga/NME1-NME2/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/downreg/healthy_vs_tcga/NME1-NME2/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

# NO COMPARISON - ALL 0 counts
# CIDEB
# Make survival plot
survival_df$CIDEB_median <- NA
survival_df$CIDEB_median[survival_df$antilog_CIDEB <= quantile(survival_df$antilog_CIDEB, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$CIDEB_median[survival_df$antilog_CIDEB > quantile(survival_df$antilog_CIDEB, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$CIDEB_25p <- NA
survival_df$CIDEB_25p[survival_df$antilog_CIDEB <= quantile(survival_df$antilog_CIDEB, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$CIDEB_25p[survival_df$antilog_CIDEB > quantile(survival_df$antilog_CIDEB, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$CIDEB_75p <- NA
survival_df$CIDEB_75p[survival_df$antilog_CIDEB <= quantile(survival_df$antilog_CIDEB, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$CIDEB_75p[survival_df$antilog_CIDEB > quantile(survival_df$antilog_CIDEB, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$CIDEB_66p <- NA
survival_df$CIDEB_66p[survival_df$antilog_CIDEB <= quantile(survival_df$antilog_CIDEB, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$CIDEB_66p[survival_df$antilog_CIDEB > quantile(survival_df$antilog_CIDEB, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$CIDEB_33p <- NA
survival_df$CIDEB_33p[survival_df$antilog_CIDEB <= quantile(survival_df$antilog_CIDEB, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CIDEB_33p[survival_df$antilog_CIDEB > quantile(survival_df$antilog_CIDEB, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$CIDEB_33_66 <- NA
survival_df$CIDEB_33_66[survival_df$antilog_CIDEB <= quantile(survival_df$antilog_CIDEB, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CIDEB_33_66[survival_df$antilog_CIDEB >= quantile(survival_df$antilog_CIDEB, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - CIDEB
for(strata in c("CIDEB_median", "CIDEB_25p", "CIDEB_75p", "CIDEB_66p", "CIDEB_33p", "CIDEB_33_66")){
  
  H_gTxt <- paste("HIGH CIDEB", sep="")
  L_gTxt <- paste("LOW CIDEB", sep="")
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
  
  tiff(filename = paste("survival_plots/downreg/healthy_vs_tcga/CIDEB/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/downreg/healthy_vs_tcga/CIDEB/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/downreg/healthy_vs_tcga/CIDEB/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}


# NO COMPARISON - ALL 0 counts
# EVA1B
# Make survival plot
survival_df$EVA1B_median <- NA
survival_df$EVA1B_median[survival_df$antilog_EVA1B <= quantile(survival_df$antilog_EVA1B, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$EVA1B_median[survival_df$antilog_EVA1B > quantile(survival_df$antilog_EVA1B, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$EVA1B_25p <- NA
survival_df$EVA1B_25p[survival_df$antilog_EVA1B <= quantile(survival_df$antilog_EVA1B, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$EVA1B_25p[survival_df$antilog_EVA1B > quantile(survival_df$antilog_EVA1B, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$EVA1B_75p <- NA
survival_df$EVA1B_75p[survival_df$antilog_EVA1B <= quantile(survival_df$antilog_EVA1B, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$EVA1B_75p[survival_df$antilog_EVA1B > quantile(survival_df$antilog_EVA1B, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$EVA1B_66p <- NA
survival_df$EVA1B_66p[survival_df$antilog_EVA1B <= quantile(survival_df$antilog_EVA1B, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$EVA1B_66p[survival_df$antilog_EVA1B > quantile(survival_df$antilog_EVA1B, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$EVA1B_33p <- NA
survival_df$EVA1B_33p[survival_df$antilog_EVA1B <= quantile(survival_df$antilog_EVA1B, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$EVA1B_33p[survival_df$antilog_EVA1B > quantile(survival_df$antilog_EVA1B, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$EVA1B_33_66 <- NA
survival_df$EVA1B_33_66[survival_df$antilog_EVA1B <= quantile(survival_df$antilog_EVA1B, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$EVA1B_33_66[survival_df$antilog_EVA1B >= quantile(survival_df$antilog_EVA1B, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - EVA1B
for(strata in c("EVA1B_median", "EVA1B_25p", "EVA1B_75p", "EVA1B_66p", "EVA1B_33p", "EVA1B_33_66")){
  
  H_gTxt <- paste("HIGH EVA1B", sep="")
  L_gTxt <- paste("LOW EVA1B", sep="")
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
  
  tiff(filename = paste("survival_plots/downreg/healthy_vs_tcga/EVA1B/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/downreg/healthy_vs_tcga/EVA1B/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/downreg/healthy_vs_tcga/EVA1B/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

# RGPD6
# Make survival plot
survival_df$RGPD6_median <- NA
survival_df$RGPD6_median[survival_df$antilog_RGPD6 <= quantile(survival_df$antilog_RGPD6, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$RGPD6_median[survival_df$antilog_RGPD6 > quantile(survival_df$antilog_RGPD6, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$RGPD6_25p <- NA
survival_df$RGPD6_25p[survival_df$antilog_RGPD6 <= quantile(survival_df$antilog_RGPD6, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$RGPD6_25p[survival_df$antilog_RGPD6 > quantile(survival_df$antilog_RGPD6, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$RGPD6_75p <- NA
survival_df$RGPD6_75p[survival_df$antilog_RGPD6 <= quantile(survival_df$antilog_RGPD6, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$RGPD6_75p[survival_df$antilog_RGPD6 > quantile(survival_df$antilog_RGPD6, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$RGPD6_66p <- NA
survival_df$RGPD6_66p[survival_df$antilog_RGPD6 <= quantile(survival_df$antilog_RGPD6, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$RGPD6_66p[survival_df$antilog_RGPD6 > quantile(survival_df$antilog_RGPD6, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$RGPD6_33p <- NA
survival_df$RGPD6_33p[survival_df$antilog_RGPD6 <= quantile(survival_df$antilog_RGPD6, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$RGPD6_33p[survival_df$antilog_RGPD6 > quantile(survival_df$antilog_RGPD6, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$RGPD6_33_66 <- NA
survival_df$RGPD6_33_66[survival_df$antilog_RGPD6 <= quantile(survival_df$antilog_RGPD6, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$RGPD6_33_66[survival_df$antilog_RGPD6 >= quantile(survival_df$antilog_RGPD6, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - RGPD6
for(strata in c("RGPD6_median", "RGPD6_25p", "RGPD6_75p", "RGPD6_66p", "RGPD6_33p", "RGPD6_33_66")){
  
  H_gTxt <- paste("HIGH RGPD6", sep="")
  L_gTxt <- paste("LOW RGPD6", sep="")
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
  
  tiff(filename = paste("survival_plots/downreg/healthy_vs_tcga/RGPD6/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/downreg/healthy_vs_tcga/RGPD6/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/downreg/healthy_vs_tcga/RGPD6/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

# NO COMPARISON - ALL 0 counts
# AD000671.1
# Make survival plot
survival_df$AD000671.1_median <- NA
survival_df$AD000671.1_median[survival_df$antilog_AD000671.1 <= quantile(survival_df$antilog_AD000671.1, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$AD000671.1_median[survival_df$antilog_AD000671.1 > quantile(survival_df$antilog_AD000671.1, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$AD000671.1_25p <- NA
survival_df$AD000671.1_25p[survival_df$antilog_AD000671.1 <= quantile(survival_df$antilog_AD000671.1, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$AD000671.1_25p[survival_df$antilog_AD000671.1 > quantile(survival_df$antilog_AD000671.1, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$AD000671.1_75p <- NA
survival_df$AD000671.1_75p[survival_df$antilog_AD000671.1 <= quantile(survival_df$antilog_AD000671.1, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$AD000671.1_75p[survival_df$antilog_AD000671.1 > quantile(survival_df$antilog_AD000671.1, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$AD000671.1_66p <- NA
survival_df$AD000671.1_66p[survival_df$antilog_AD000671.1 <= quantile(survival_df$antilog_AD000671.1, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$AD000671.1_66p[survival_df$antilog_AD000671.1 > quantile(survival_df$antilog_AD000671.1, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$AD000671.1_33p <- NA
survival_df$AD000671.1_33p[survival_df$antilog_AD000671.1 <= quantile(survival_df$antilog_AD000671.1, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$AD000671.1_33p[survival_df$antilog_AD000671.1 > quantile(survival_df$antilog_AD000671.1, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$AD000671.1_33_66 <- NA
survival_df$AD000671.1_33_66[survival_df$antilog_AD000671.1 <= quantile(survival_df$antilog_AD000671.1, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$AD000671.1_33_66[survival_df$antilog_AD000671.1 >= quantile(survival_df$antilog_AD000671.1, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - AD000671.1
for(strata in c("AD000671.1_median", "AD000671.1_25p", "AD000671.1_75p", "AD000671.1_66p", "AD000671.1_33p", "AD000671.1_33_66")){
  
  H_gTxt <- paste("HIGH AD000671.1", sep="")
  L_gTxt <- paste("LOW AD000671.1", sep="")
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
  
  tiff(filename = paste("survival_plots/downreg/healthy_vs_tcga/AD000671.1/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/downreg/healthy_vs_tcga/AD000671.1/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/downreg/healthy_vs_tcga/AD000671.1/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

# NO COMPARISON - ALL 0 counts
# CEMP1
# Make survival plot
survival_df$CEMP1_median <- NA
survival_df$CEMP1_median[survival_df$antilog_CEMP1 <= quantile(survival_df$antilog_CEMP1, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$CEMP1_median[survival_df$antilog_CEMP1 > quantile(survival_df$antilog_CEMP1, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$CEMP1_25p <- NA
survival_df$CEMP1_25p[survival_df$antilog_CEMP1 <= quantile(survival_df$antilog_CEMP1, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$CEMP1_25p[survival_df$antilog_CEMP1 > quantile(survival_df$antilog_CEMP1, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$CEMP1_75p <- NA
survival_df$CEMP1_75p[survival_df$antilog_CEMP1 <= quantile(survival_df$antilog_CEMP1, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$CEMP1_75p[survival_df$antilog_CEMP1 > quantile(survival_df$antilog_CEMP1, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$CEMP1_66p <- NA
survival_df$CEMP1_66p[survival_df$antilog_CEMP1 <= quantile(survival_df$antilog_CEMP1, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$CEMP1_66p[survival_df$antilog_CEMP1 > quantile(survival_df$antilog_CEMP1, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$CEMP1_33p <- NA
survival_df$CEMP1_33p[survival_df$antilog_CEMP1 <= quantile(survival_df$antilog_CEMP1, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CEMP1_33p[survival_df$antilog_CEMP1 > quantile(survival_df$antilog_CEMP1, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$CEMP1_33_66 <- NA
survival_df$CEMP1_33_66[survival_df$antilog_CEMP1 <= quantile(survival_df$antilog_CEMP1, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$CEMP1_33_66[survival_df$antilog_CEMP1 >= quantile(survival_df$antilog_CEMP1, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - CEMP1
for(strata in c("CEMP1_median", "CEMP1_25p", "CEMP1_75p", "CEMP1_66p", "CEMP1_33p", "CEMP1_33_66")){
  
  H_gTxt <- paste("HIGH CEMP1", sep="")
  L_gTxt <- paste("LOW CEMP1", sep="")
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
  
  tiff(filename = paste("survival_plots/downreg/healthy_vs_tcga/CEMP1/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/downreg/healthy_vs_tcga/CEMP1/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/downreg/healthy_vs_tcga/CEMP1/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

# NO COMPARISON - ALL 0 counts
# IGLC1
# Make survival plot
survival_df$IGLC1_median <- NA
survival_df$IGLC1_median[survival_df$antilog_IGLC1 <= quantile(survival_df$antilog_IGLC1, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$IGLC1_median[survival_df$antilog_IGLC1 > quantile(survival_df$antilog_IGLC1, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$IGLC1_25p <- NA
survival_df$IGLC1_25p[survival_df$antilog_IGLC1 <= quantile(survival_df$antilog_IGLC1, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$IGLC1_25p[survival_df$antilog_IGLC1 > quantile(survival_df$antilog_IGLC1, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$IGLC1_75p <- NA
survival_df$IGLC1_75p[survival_df$antilog_IGLC1 <= quantile(survival_df$antilog_IGLC1, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$IGLC1_75p[survival_df$antilog_IGLC1 > quantile(survival_df$antilog_IGLC1, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$IGLC1_66p <- NA
survival_df$IGLC1_66p[survival_df$antilog_IGLC1 <= quantile(survival_df$antilog_IGLC1, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$IGLC1_66p[survival_df$antilog_IGLC1 > quantile(survival_df$antilog_IGLC1, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$IGLC1_33p <- NA
survival_df$IGLC1_33p[survival_df$antilog_IGLC1 <= quantile(survival_df$antilog_IGLC1, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$IGLC1_33p[survival_df$antilog_IGLC1 > quantile(survival_df$antilog_IGLC1, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$IGLC1_33_66 <- NA
survival_df$IGLC1_33_66[survival_df$antilog_IGLC1 <= quantile(survival_df$antilog_IGLC1, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$IGLC1_33_66[survival_df$antilog_IGLC1 >= quantile(survival_df$antilog_IGLC1, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - IGLC1
for(strata in c("IGLC1_median", "IGLC1_25p", "IGLC1_75p", "IGLC1_66p", "IGLC1_33p", "IGLC1_33_66")){
  
  H_gTxt <- paste("HIGH IGLC1", sep="")
  L_gTxt <- paste("LOW IGLC1", sep="")
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
  
  tiff(filename = paste("survival_plots/downreg/healthy_vs_tcga/IGLC1/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/downreg/healthy_vs_tcga/IGLC1/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/downreg/healthy_vs_tcga/IGLC1/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}


# NPIPA7
# Make survival plot
survival_df$NPIPA7_median <- NA
survival_df$NPIPA7_median[survival_df$antilog_NPIPA7 <= quantile(survival_df$antilog_NPIPA7, 0.5, na.rm=TRUE)] <- "LOW"
survival_df$NPIPA7_median[survival_df$antilog_NPIPA7 > quantile(survival_df$antilog_NPIPA7, 0.5, na.rm=TRUE)] <- "HIGH"

survival_df$NPIPA7_25p <- NA
survival_df$NPIPA7_25p[survival_df$antilog_NPIPA7 <= quantile(survival_df$antilog_NPIPA7, 0.25, na.rm=TRUE)] <- "LOW"
survival_df$NPIPA7_25p[survival_df$antilog_NPIPA7 > quantile(survival_df$antilog_NPIPA7, 0.25, na.rm=TRUE)] <- "HIGH"

survival_df$NPIPA7_75p <- NA
survival_df$NPIPA7_75p[survival_df$antilog_NPIPA7 <= quantile(survival_df$antilog_NPIPA7, 0.75, na.rm=TRUE)] <- "LOW"
survival_df$NPIPA7_75p[survival_df$antilog_NPIPA7 > quantile(survival_df$antilog_NPIPA7, 0.75, na.rm=TRUE)] <- "HIGH"

survival_df$NPIPA7_66p <- NA
survival_df$NPIPA7_66p[survival_df$antilog_NPIPA7 <= quantile(survival_df$antilog_NPIPA7, 0.66, na.rm=TRUE)] <- "LOW"
survival_df$NPIPA7_66p[survival_df$antilog_NPIPA7 > quantile(survival_df$antilog_NPIPA7, 0.66, na.rm=TRUE)] <- "HIGH"

survival_df$NPIPA7_33p <- NA
survival_df$NPIPA7_33p[survival_df$antilog_NPIPA7 <= quantile(survival_df$antilog_NPIPA7, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$NPIPA7_33p[survival_df$antilog_NPIPA7 > quantile(survival_df$antilog_NPIPA7, 0.33, na.rm=TRUE)] <- "HIGH"

survival_df$NPIPA7_33_66 <- NA
survival_df$NPIPA7_33_66[survival_df$antilog_NPIPA7 <= quantile(survival_df$antilog_NPIPA7, 0.33, na.rm=TRUE)] <- "LOW"
survival_df$NPIPA7_33_66[survival_df$antilog_NPIPA7 >= quantile(survival_df$antilog_NPIPA7, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "TCGA"

# SURVIVAL CURVES - NPIPA7
for(strata in c("NPIPA7_median", "NPIPA7_25p", "NPIPA7_75p", "NPIPA7_66p", "NPIPA7_33p", "NPIPA7_33_66")){
  
  H_gTxt <- paste("HIGH NPIPA7", sep="")
  L_gTxt <- paste("LOW NPIPA7", sep="")
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
  
  tiff(filename = paste("survival_plots/downreg/healthy_vs_tcga/NPIPA7/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
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
  
  tiff(filename = paste ("survival_plots/downreg/healthy_vs_tcga/NPIPA7/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
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
  
  tiff(filename =  paste("survival_plots/downreg/healthy_vs_tcga/NPIPA7/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}



# Cox regression
meta <- c("Row.names", "X_PATIENT")

# Median
strata <- c("antilog_RGPD6", "RGPD6_median", "antilog_NPIPA7", "NPIPA7_median"
            )

# OS
os <- c("OS", "OS_months")
os_median_df <- survival_df %>% 
  select(all_of(c(meta, os, strata)))
os_median_df <- subset(os_median_df,  !is.na(OS) & !is.na(OS_months))
cox <- coxph(Surv(OS_months, OS) ~ RGPD6_median + NPIPA7_median, data= os_median_df)
summary(cox)
# Remove RGPD6
cox <- coxph(Surv(OS_months, OS) ~ NPIPA7_median, data= os_median_df)
summary(cox)

# DSS
dss <- c("DSS", "DSS_months")
dss_median_df <- survival_df %>% 
  select(all_of(c(meta, dss, strata)))
dss_median_df <- subset(dss_median_df,  !is.na(DSS) & !is.na(DSS_months))
cox <- coxph(Surv(DSS_months, DSS) ~ RGPD6_median + NPIPA7_median, data= dss_median_df)
summary(cox)
# Remove RGPD6
cox <- coxph(Surv(DSS_months, DSS) ~ NPIPA7_median, data= dss_median_df)
summary(cox)

# PFI
pfi <- c("PFI", "PFI_months")
pfi_median_df <- survival_df %>% 
  select(all_of(c(meta, pfi, strata)))
pfi_median_df <- subset(pfi_median_df,  !is.na(PFI) & !is.na(PFI_months))
cox <- coxph(Surv(PFI_months, PFI) ~ RGPD6_median + NPIPA7_median, data = pfi_median_df)
summary(cox)
# remove rgpd6
cox <- coxph(Surv(PFI_months, PFI) ~ NPIPA7_median, data = pfi_median_df)
summary(cox)
