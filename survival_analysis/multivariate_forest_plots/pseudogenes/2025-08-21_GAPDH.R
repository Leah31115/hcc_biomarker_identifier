library(tidyverse)
library(ggpubr) # for shapiro wilk test & box plots, pvalues
library(reshape2) # paired t tests, melt

# P values in ggplot figs
library(rstatix)

# survival
library(survival)
library(survminer)
library(tidytidbits)
library(survivalAnalysis)

# Clear environment
rm(list=ls())

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

GAPDH_geneid <- c("ENSG00000232718.1",
           "ENSG00000251013.1",
           "ENSG00000225502.2",
           "ENSG00000250982.2",
           "ENSG00000218582.2",
           "ENSG00000216624.3",
           "ENSG00000213376.4",
           "ENSG00000213352.3",
           "ENSG00000226540.2",
           "ENSG00000223460.1",
           "ENSG00000233876.2",
           "ENSG00000249489.1",
           "ENSG00000249210.1",
           "ENSG00000231739.2",
           "ENSG00000236056.1",
           "ENSG00000231072.1",
           "ENSG00000235587.2",
           "ENSG00000249840.2",
           "ENSG00000213594.4",
           "ENSG00000236993.2",
           "ENSG00000248180.1",
           "ENSG00000229570.2",
           "ENSG00000224333.1",
           "ENSG00000228232.1",
           "ENSG00000248626.1",
           "ENSG00000111640.15",
           "ENSG00000219930.4",
           "ENSG00000234005.3",
           "ENSG00000239873.2",
           "ENSG00000248415.1",
           "ENSG00000226443.3",
           "ENSG00000224055.1",
           "ENSG00000231907.2"
)


# GAPDH
GAPDH_df <- TCGA_samples %>%
  filter(Ensembl_ID %in% c(GAPDH_geneid))

# Rename row.names column and make the rows the gene names
rownames(GAPDH_df) <- GAPDH_df$Ensembl_ID
GAPDH_df <- GAPDH_df[,!(names(GAPDH_df) %in% "Ensembl_ID")]
# Filter for used cancer samples
GAPDH_df <- as.data.frame(t(GAPDH_df))
GAPDH_df <- tibble::rownames_to_column(GAPDH_df, "sample")
GAPDH_df <- GAPDH_df %>%
  filter(sample %in% c(cancer_sample_names))
# Make sample the row names
rownames(GAPDH_df) <- GAPDH_df$sample
GAPDH_df <- GAPDH_df[,!(names(GAPDH_df) %in% "sample")]

# Add meta
rownames(tcga_meta) <- tcga_meta$sample
tcga_meta <- tcga_meta[,!(names(tcga_meta) %in% c("sample", "samples"))]
GAPDH_df <- merge(GAPDH_df, tcga_meta,
                  by = 'row.names', all = TRUE)
# Remove A suffix to pair samples with survival data
names(GAPDH_df)[names(GAPDH_df) == "Row.names"] <- "sample"
GAPDH_df$sample <- substring(GAPDH_df$sample, 1, nchar(GAPDH_df$sample)-1)
rownames(GAPDH_df) <- GAPDH_df$sample
GAPDH_df <- GAPDH_df[,!(names(GAPDH_df) %in% "sample")]

colnames(GAPDH_df)


names(GAPDH_df)[names(GAPDH_df) == "primary_diagnosis.diagnoses"] <- "primary_diagnosis"
names(GAPDH_df)[names(GAPDH_df) == "age_at_earliest_diagnosis_in_years.diagnoses.xena_derived"] <- "age_at_earliest_diagnosis"
names(GAPDH_df)[names(GAPDH_df) == "ajcc_pathologic_stage.diagnoses"] <- "all_stages"
names(GAPDH_df)[names(GAPDH_df) == "ajcc_pathologic_m.diagnoses"] <- "M"
names(GAPDH_df)[names(GAPDH_df) == "ajcc_pathologic_n.diagnoses"] <- "N"
names(GAPDH_df)[names(GAPDH_df) == "gender.demographic"] <- "sex"
names(GAPDH_df)[names(GAPDH_df) == "ajcc_pathologic_t.diagnoses"] <- "all_T"
colnames(GAPDH_df)

# Round ages down so you just get their year
GAPDH_df$age_at_earliest_diagnosis <- lapply(GAPDH_df$age_at_earliest_diagnosis, floor)
# Assign age category
GAPDH_df$age_range <- NA
GAPDH_df$age_range <- ifelse(GAPDH_df$age_at_earliest_diagnosis >= 80, "80+",
                             ifelse(GAPDH_df$age_at_earliest_diagnosis >= 60, "60-79",
                                    ifelse(GAPDH_df$age_at_earliest_diagnosis >= 40, "40-59",
                                           ifelse(GAPDH_df$age_at_earliest_diagnosis >= 18, "18-39", "0-17")
                                    )
                             )
)


# Group subgroups in just one group e.g stage 4A = stage4
GAPDH_df$stage <- GAPDH_df$all_stages
GAPDH_df$stage[GAPDH_df$stage == "Stage IIIA" | GAPDH_df$stage == "Stage IIIB" | GAPDH_df$stage == "Stage IIIC"] <- "Stage III"
GAPDH_df$stage[GAPDH_df$stage == "Stage IVA" | GAPDH_df$stage == "Stage IVB"] <- "Stage IV"
GAPDH_df$T <- GAPDH_df$all_T
GAPDH_df$T[GAPDH_df$T == "T2a" | GAPDH_df$T == "T2b"] <- "T2"
GAPDH_df$T[GAPDH_df$T == "T3a" | GAPDH_df$T == "T3b"] <- "T3"



# Add log prefix to transcript names
GAPDH_df <- GAPDH_df %>% 
  rename_with(function(x) paste0("log_", x), all_of(GAPDH_geneid))
colnames(GAPDH_df)


# Antilog gene count function
antilog_counts <- function(tpm){
  antilog_count <- (2^tpm) - 0.001
  return(antilog_count)
}

# Antilog
GAPDH_df["GAPDHP48"] <- lapply(GAPDH_df["log_ENSG00000232718.1"], antilog_counts)
GAPDH_df["GAPDHP62"] <- lapply(GAPDH_df["log_ENSG00000251013.1"], antilog_counts)
GAPDH_df["GAPDHP16"] <- lapply(GAPDH_df["log_ENSG00000225502.2"], antilog_counts)
GAPDH_df["GAPDHP35"] <- lapply(GAPDH_df["log_ENSG00000250982.2"], antilog_counts)
GAPDH_df["GAPDHP63"] <- lapply(GAPDH_df["log_ENSG00000218582.2"], antilog_counts)
GAPDH_df["GAPDHP72"] <- lapply(GAPDH_df["log_ENSG00000216624.3"], antilog_counts)
GAPDH_df["GAPDHP71"] <- lapply(GAPDH_df["log_ENSG00000213376.4"], antilog_counts)
GAPDH_df["GAPDHP44"] <- lapply(GAPDH_df["log_ENSG00000213352.3"], antilog_counts)
GAPDH_df["GAPDHP73"] <- lapply(GAPDH_df["log_ENSG00000226540.2"], antilog_counts)
GAPDH_df["GAPDHP69"] <- lapply(GAPDH_df["log_ENSG00000223460.1"], antilog_counts)
GAPDH_df["GAPDHP68"] <- lapply(GAPDH_df["log_ENSG00000233876.2"], antilog_counts)
GAPDH_df["GAPDHP70"] <- lapply(GAPDH_df["log_ENSG00000249489.1"], antilog_counts)
GAPDH_df["GAPDHP38"] <- lapply(GAPDH_df["log_ENSG00000249210.1"], antilog_counts)
GAPDH_df["GAPDHP59"] <- lapply(GAPDH_df["log_ENSG00000231739.2"], antilog_counts)
GAPDH_df["GAPDHP14"] <- lapply(GAPDH_df["log_ENSG00000236056.1"], antilog_counts)
GAPDH_df["GAPDHP64"] <- lapply(GAPDH_df["log_ENSG00000231072.1"], antilog_counts)
GAPDH_df["GAPDHP65"] <- lapply(GAPDH_df["log_ENSG00000235587.2"], antilog_counts)
GAPDH_df["GAPDHP76"] <- lapply(GAPDH_df["log_ENSG00000249840.2"], antilog_counts)
GAPDH_df["GAPDHP25"] <- lapply(GAPDH_df["log_ENSG00000213594.4"], antilog_counts)
GAPDH_df["GAPDHP21"] <- lapply(GAPDH_df["log_ENSG00000236993.2"], antilog_counts)
GAPDH_df["GAPDHP60"] <- lapply(GAPDH_df["log_ENSG00000248180.1"], antilog_counts)
GAPDH_df["GAPDHP58"] <- lapply(GAPDH_df["log_ENSG00000229570.2"], antilog_counts)
GAPDH_df["GAPDHP20"] <- lapply(GAPDH_df["log_ENSG00000224333.1"], antilog_counts)
GAPDH_df["GAPDHP1"] <- lapply(GAPDH_df["log_ENSG00000228232.1"], antilog_counts)
GAPDH_df["GAPDHP40"] <- lapply(GAPDH_df["log_ENSG00000248626.1"], antilog_counts)
GAPDH_df["GAPDH"] <- lapply(GAPDH_df["log_ENSG00000111640.15"], antilog_counts)
GAPDH_df["GAPDHP67"] <- lapply(GAPDH_df["log_ENSG00000219930.4"], antilog_counts)
GAPDH_df["GAPDHP22"] <- lapply(GAPDH_df["log_ENSG00000234005.3"], antilog_counts)
GAPDH_df["GAPDHP27"] <- lapply(GAPDH_df["log_ENSG00000239873.2"], antilog_counts)
GAPDH_df["GAPDHP61"] <- lapply(GAPDH_df["log_ENSG00000248415.1"], antilog_counts)
GAPDH_df["GAPDHP32"] <- lapply(GAPDH_df["log_ENSG00000226443.3"], antilog_counts)
GAPDH_df["GAPDHP55"] <- lapply(GAPDH_df["log_ENSG00000224055.1"], antilog_counts)
GAPDH_df["GAPDHP37"] <- lapply(GAPDH_df["log_ENSG00000231907.2"], antilog_counts)

# subset for antilogged counts
GAPDH_df_anti <- GAPDH_df %>% 
  select(GAPDHP48,
         GAPDHP62,
         GAPDHP16,
         GAPDHP35,
         GAPDHP63,
         GAPDHP72,
         GAPDHP71,
         GAPDHP44,
         GAPDHP73,
         GAPDHP69,
         GAPDHP68,
         GAPDHP70,
         GAPDHP38,
         GAPDHP59,
         GAPDHP14,
         GAPDHP64,
         GAPDHP65,
         GAPDHP76,
         GAPDHP25,
         GAPDHP21,
         GAPDHP60,
         GAPDHP58,
         GAPDHP20,
         GAPDHP1,
         GAPDHP40,
         GAPDH,
         GAPDHP67,
         GAPDHP22,
         GAPDHP27,
         GAPDHP61,
         GAPDHP32,
         GAPDHP55,
         GAPDHP37,
         all_stages,
         all_T,
         M,
         N,
         sex,
         T,
         age_range
  )

# Transcript survival analysis
# Filter survival data for used samples
cancer_sample_names = row.names(GAPDH_df)
cancer_sample_names

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
GAPDH_survival_df <- merge(GAPDH_df_anti, filtered_tcga_survival, by = "row.names", 
                           all.x = TRUE) # all.x will give NA for missing values

# Convert days to months
# 1 day = 0.0328767 months
days_to_months_converter <- function(no_days){
  no_months <- round(no_days*0.0328767, 1)
  return(no_months)
}

GAPDH_survival_df["OS_months"] <- lapply(GAPDH_survival_df["OS.time"], days_to_months_converter)
GAPDH_survival_df["DSS_months"] <- lapply(GAPDH_survival_df["DSS.time"], days_to_months_converter)

# Investigating median gene expression category on survival
GAPDH_survival_df$GAPDHP48_median <- NA
GAPDH_survival_df$GAPDHP48_median[GAPDH_survival_df$GAPDHP48 <= quantile(GAPDH_survival_df$GAPDHP48, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP48_median[GAPDH_survival_df$GAPDHP48 > quantile(GAPDH_survival_df$GAPDHP48, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP62_median <- NA
GAPDH_survival_df$GAPDHP62_median[GAPDH_survival_df$GAPDHP62 <= quantile(GAPDH_survival_df$GAPDHP62, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP62_median[GAPDH_survival_df$GAPDHP62 > quantile(GAPDH_survival_df$GAPDHP62, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP16_median <- NA
GAPDH_survival_df$GAPDHP16_median[GAPDH_survival_df$GAPDHP16 <= quantile(GAPDH_survival_df$GAPDHP16, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP16_median[GAPDH_survival_df$GAPDHP16 > quantile(GAPDH_survival_df$GAPDHP16, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP35_median <- NA
GAPDH_survival_df$GAPDHP35_median[GAPDH_survival_df$GAPDHP35 <= quantile(GAPDH_survival_df$GAPDHP35, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP35_median[GAPDH_survival_df$GAPDHP35 > quantile(GAPDH_survival_df$GAPDHP35, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP63_median <- NA
GAPDH_survival_df$GAPDHP63_median[GAPDH_survival_df$GAPDHP63 <= quantile(GAPDH_survival_df$GAPDHP63, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP63_median[GAPDH_survival_df$GAPDHP63 > quantile(GAPDH_survival_df$GAPDHP63, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP72_median <- NA
GAPDH_survival_df$GAPDHP72_median[GAPDH_survival_df$GAPDHP72 <= quantile(GAPDH_survival_df$GAPDHP72, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP72_median[GAPDH_survival_df$GAPDHP72 > quantile(GAPDH_survival_df$GAPDHP72, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP71_median <- NA
GAPDH_survival_df$GAPDHP71_median[GAPDH_survival_df$GAPDHP71 <= quantile(GAPDH_survival_df$GAPDHP71, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP71_median[GAPDH_survival_df$GAPDHP71 > quantile(GAPDH_survival_df$GAPDHP71, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP44_median <- NA
GAPDH_survival_df$GAPDHP44_median[GAPDH_survival_df$GAPDHP44 <= quantile(GAPDH_survival_df$GAPDHP44, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP44_median[GAPDH_survival_df$GAPDHP44 > quantile(GAPDH_survival_df$GAPDHP44, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP73_median <- NA
GAPDH_survival_df$GAPDHP73_median[GAPDH_survival_df$GAPDHP73 <= quantile(GAPDH_survival_df$GAPDHP73, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP73_median[GAPDH_survival_df$GAPDHP73 > quantile(GAPDH_survival_df$GAPDHP73, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP69_median <- NA
GAPDH_survival_df$GAPDHP69_median[GAPDH_survival_df$GAPDHP69 <= quantile(GAPDH_survival_df$GAPDHP69, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP69_median[GAPDH_survival_df$GAPDHP69 > quantile(GAPDH_survival_df$GAPDHP69, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP68_median <- NA
GAPDH_survival_df$GAPDHP68_median[GAPDH_survival_df$GAPDHP68 <= quantile(GAPDH_survival_df$GAPDHP68, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP68_median[GAPDH_survival_df$GAPDHP68 > quantile(GAPDH_survival_df$GAPDHP68, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP70_median <- NA
GAPDH_survival_df$GAPDHP70_median[GAPDH_survival_df$GAPDHP70 <= quantile(GAPDH_survival_df$GAPDHP70, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP70_median[GAPDH_survival_df$GAPDHP70 > quantile(GAPDH_survival_df$GAPDHP70, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP38_median <- NA
GAPDH_survival_df$GAPDHP38_median[GAPDH_survival_df$GAPDHP38 <= quantile(GAPDH_survival_df$GAPDHP38, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP38_median[GAPDH_survival_df$GAPDHP38 > quantile(GAPDH_survival_df$GAPDHP38, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP59_median <- NA
GAPDH_survival_df$GAPDHP59_median[GAPDH_survival_df$GAPDHP59 <= quantile(GAPDH_survival_df$GAPDHP59, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP59_median[GAPDH_survival_df$GAPDHP59 > quantile(GAPDH_survival_df$GAPDHP59, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP14_median <- NA
GAPDH_survival_df$GAPDHP14_median[GAPDH_survival_df$GAPDHP14 <= quantile(GAPDH_survival_df$GAPDHP14, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP14_median[GAPDH_survival_df$GAPDHP14 > quantile(GAPDH_survival_df$GAPDHP14, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP64_median <- NA
GAPDH_survival_df$GAPDHP64_median[GAPDH_survival_df$GAPDHP64 <= quantile(GAPDH_survival_df$GAPDHP64, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP64_median[GAPDH_survival_df$GAPDHP64 > quantile(GAPDH_survival_df$GAPDHP64, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP65_median <- NA
GAPDH_survival_df$GAPDHP65_median[GAPDH_survival_df$GAPDHP65 <= quantile(GAPDH_survival_df$GAPDHP65, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP65_median[GAPDH_survival_df$GAPDHP65 > quantile(GAPDH_survival_df$GAPDHP65, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP76_median <- NA
GAPDH_survival_df$GAPDHP76_median[GAPDH_survival_df$GAPDHP76 <= quantile(GAPDH_survival_df$GAPDHP76, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP76_median[GAPDH_survival_df$GAPDHP76 > quantile(GAPDH_survival_df$GAPDHP76, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP25_median <- NA
GAPDH_survival_df$GAPDHP25_median[GAPDH_survival_df$GAPDHP25 <= quantile(GAPDH_survival_df$GAPDHP25, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP25_median[GAPDH_survival_df$GAPDHP25 > quantile(GAPDH_survival_df$GAPDHP25, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP21_median <- NA
GAPDH_survival_df$GAPDHP21_median[GAPDH_survival_df$GAPDHP21 <= quantile(GAPDH_survival_df$GAPDHP21, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP21_median[GAPDH_survival_df$GAPDHP21 > quantile(GAPDH_survival_df$GAPDHP21, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP60_median <- NA
GAPDH_survival_df$GAPDHP60_median[GAPDH_survival_df$GAPDHP60 <= quantile(GAPDH_survival_df$GAPDHP60, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP60_median[GAPDH_survival_df$GAPDHP60 > quantile(GAPDH_survival_df$GAPDHP60, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP58_median <- NA
GAPDH_survival_df$GAPDHP58_median[GAPDH_survival_df$GAPDHP58 <= quantile(GAPDH_survival_df$GAPDHP58, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP58_median[GAPDH_survival_df$GAPDHP58 > quantile(GAPDH_survival_df$GAPDHP58, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP20_median <- NA
GAPDH_survival_df$GAPDHP20_median[GAPDH_survival_df$GAPDHP20 <= quantile(GAPDH_survival_df$GAPDHP20, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP20_median[GAPDH_survival_df$GAPDHP20 > quantile(GAPDH_survival_df$GAPDHP20, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP1_median <- NA
GAPDH_survival_df$GAPDHP1_median[GAPDH_survival_df$GAPDHP1 <= quantile(GAPDH_survival_df$GAPDHP1, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP1_median[GAPDH_survival_df$GAPDHP1 > quantile(GAPDH_survival_df$GAPDHP1, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP40_median <- NA
GAPDH_survival_df$GAPDHP40_median[GAPDH_survival_df$GAPDHP40 <= quantile(GAPDH_survival_df$GAPDHP40, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP40_median[GAPDH_survival_df$GAPDHP40 > quantile(GAPDH_survival_df$GAPDHP40, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDH_median <- NA
GAPDH_survival_df$GAPDH_median[GAPDH_survival_df$GAPDH <= quantile(GAPDH_survival_df$GAPDH, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDH_median[GAPDH_survival_df$GAPDH > quantile(GAPDH_survival_df$GAPDH, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP67_median <- NA
GAPDH_survival_df$GAPDHP67_median[GAPDH_survival_df$GAPDHP67 <= quantile(GAPDH_survival_df$GAPDHP67, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP67_median[GAPDH_survival_df$GAPDHP67 > quantile(GAPDH_survival_df$GAPDHP67, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP22_median <- NA
GAPDH_survival_df$GAPDHP22_median[GAPDH_survival_df$GAPDHP22 <= quantile(GAPDH_survival_df$GAPDHP22, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP22_median[GAPDH_survival_df$GAPDHP22 > quantile(GAPDH_survival_df$GAPDHP22, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP27_median <- NA
GAPDH_survival_df$GAPDHP27_median[GAPDH_survival_df$GAPDHP27 <= quantile(GAPDH_survival_df$GAPDHP27, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP27_median[GAPDH_survival_df$GAPDHP27 > quantile(GAPDH_survival_df$GAPDHP27, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP61_median <- NA
GAPDH_survival_df$GAPDHP61_median[GAPDH_survival_df$GAPDHP61 <= quantile(GAPDH_survival_df$GAPDHP61, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP61_median[GAPDH_survival_df$GAPDHP61 > quantile(GAPDH_survival_df$GAPDHP61, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP32_median <- NA
GAPDH_survival_df$GAPDHP32_median[GAPDH_survival_df$GAPDHP32 <= quantile(GAPDH_survival_df$GAPDHP32, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP32_median[GAPDH_survival_df$GAPDHP32 > quantile(GAPDH_survival_df$GAPDHP32, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP55_median <- NA
GAPDH_survival_df$GAPDHP55_median[GAPDH_survival_df$GAPDHP55 <= quantile(GAPDH_survival_df$GAPDHP55, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP55_median[GAPDH_survival_df$GAPDHP55 > quantile(GAPDH_survival_df$GAPDHP55, 0.5, na.rm=TRUE)] <- "HIGH"

GAPDH_survival_df$GAPDHP37_median <- NA
GAPDH_survival_df$GAPDHP37_median[GAPDH_survival_df$GAPDHP37 <= quantile(GAPDH_survival_df$GAPDHP37, 0.5, na.rm=TRUE)] <- "LOW"
GAPDH_survival_df$GAPDHP37_median[GAPDH_survival_df$GAPDHP37 > quantile(GAPDH_survival_df$GAPDHP37, 0.5, na.rm=TRUE)] <- "HIGH"

covariate_names <- c(GAPDHP62_median="GAPDHP62 median",
                     GAPDHP16_median="GAPDHP16 median",
                     GAPDHP48_median="GAPDHP48 median",
                     GAPDHP35_median="GAPDHP35 median",
                     GAPDHP63_median="GAPDHP63 median",
                     GAPDHP72_median="GAPDHP72 median",
                     GAPDHP71_median="GAPDHP71 median",
                     GAPDHP44_median="GAPDHP44 median",
                     GAPDHP73_median="GAPDHP73 median",
                     GAPDHP69_median="GAPDHP69 median",
                     GAPDHP68_median="GAPDHP68 median",
                     GAPDHP70_median="GAPDHP70 median",
                     GAPDHP38_median="GAPDHP38 median",
                     GAPDHP59_median="GAPDHP59 median",
                     GAPDHP14_median="GAPDHP14 median",
                     GAPDHP64_median="GAPDHP64 median",
                     GAPDHP65_median="GAPDHP65 median",
                     GAPDHP76_median="GAPDHP76 median",
                     GAPDHP25_median="GAPDHP25 median",
                     GAPDHP21_median="GAPDHP21 median",
                     GAPDHP60_median="GAPDHP60 median",
                     GAPDHP58_median="GAPDHP58 median",
                     GAPDHP20_median="GAPDHP20 median",
                     GAPDHP1_median="GAPDHP1 median",
                     GAPDHP40_median="GAPDHP40 median",
                     GAPDH_median="GAPDH median",
                     GAPDHP67_median="GAPDHP67 median",
                     GAPDHP22_median="GAPDHP22 median",
                     GAPDHP27_median="GAPDHP27 median",
                     GAPDHP61_median="GAPDHP61 median",
                     GAPDHP32_median="GAPDHP32 median",
                     GAPDHP55_median="GAPDHP55 median",
                     GAPDHP37_median="GAPDHP37 median"
)

# trying with median 
# OS
GAPDH_median_OS <- GAPDH_survival_df %>%
  select(Row.names, GAPDHP62_median, GAPDHP16_median, GAPDHP48_median, GAPDHP35_median, GAPDHP63_median, GAPDHP72_median, GAPDHP71_median, GAPDHP44_median, GAPDHP73_median, GAPDHP69_median, GAPDHP68_median, GAPDHP70_median, GAPDHP38_median, GAPDHP59_median, GAPDHP14_median, GAPDHP64_median, GAPDHP65_median, GAPDHP76_median, GAPDHP25_median, GAPDHP21_median, GAPDHP60_median, GAPDHP58_median, GAPDHP20_median, GAPDHP1_median, GAPDHP40_median, GAPDH_median, GAPDHP67_median, GAPDHP22_median, GAPDHP27_median, GAPDHP61_median, GAPDHP32_median, GAPDHP55_median, GAPDHP37_median, sex, OS_months, OS)
# Remove NAs
GAPDH_median_OS <- subset(GAPDH_median_OS, !is.na(OS) & !is.na(OS_months))

map(vars(GAPDHP62_median, GAPDHP16_median, GAPDHP48_median, GAPDHP35_median, GAPDHP63_median, GAPDHP72_median, GAPDHP71_median, GAPDHP44_median, GAPDHP73_median, GAPDHP69_median, GAPDHP68_median, GAPDHP70_median, GAPDHP38_median, GAPDHP59_median, GAPDHP14_median, GAPDHP64_median, GAPDHP65_median, GAPDHP76_median, GAPDHP25_median, GAPDHP21_median, GAPDHP60_median, GAPDHP58_median, GAPDHP20_median, GAPDHP1_median, GAPDHP40_median, GAPDH_median, GAPDHP67_median, GAPDHP22_median, GAPDHP27_median, GAPDHP61_median, GAPDHP32_median, GAPDHP55_median, GAPDHP37_median), function(by)
{
  analyse_multivariate(GAPDH_median_OS,
                       vars(OS_months, OS),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = covariate_names)
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="OS"),
              orderer = ~order(HR),
              labels_displayed = c("endpoint", "factor", "n"),
              ggtheme = ggplot2::theme_bw(base_size = 10))

# DSS
GAPDH_median_DSS <- GAPDH_survival_df %>%
  select(Row.names, GAPDHP62_median, GAPDHP16_median, GAPDHP48_median, GAPDHP35_median, GAPDHP63_median, GAPDHP72_median, GAPDHP71_median, GAPDHP44_median, GAPDHP73_median, GAPDHP69_median, GAPDHP68_median, GAPDHP70_median, GAPDHP38_median, GAPDHP59_median, GAPDHP14_median, GAPDHP64_median, GAPDHP65_median, GAPDHP76_median, GAPDHP25_median, GAPDHP21_median, GAPDHP60_median, GAPDHP58_median, GAPDHP20_median, GAPDHP1_median, GAPDHP40_median, GAPDH_median, GAPDHP67_median, GAPDHP22_median, GAPDHP27_median, GAPDHP61_median, GAPDHP32_median, GAPDHP55_median, GAPDHP37_median, sex, DSS_months, DSS)

# Remove NAs
GAPDH_median_DSS <- subset(GAPDH_median_DSS, !is.na(DSS) & !is.na(DSS_months))

map(vars(GAPDHP62_median, GAPDHP16_median, GAPDHP48_median, GAPDHP35_median, GAPDHP63_median, GAPDHP72_median, GAPDHP71_median, GAPDHP44_median, GAPDHP73_median, GAPDHP69_median, GAPDHP68_median, GAPDHP70_median, GAPDHP38_median, GAPDHP59_median, GAPDHP14_median, GAPDHP64_median, GAPDHP65_median, GAPDHP76_median, GAPDHP25_median, GAPDHP21_median, GAPDHP60_median, GAPDHP58_median, GAPDHP20_median, GAPDHP1_median, GAPDHP40_median, GAPDH_median, GAPDHP67_median, GAPDHP22_median, GAPDHP27_median, GAPDHP61_median, GAPDHP32_median, GAPDHP55_median, GAPDHP37_median), function(by)
{
  analyse_multivariate(GAPDH_median_DSS,
                       vars(DSS_months, DSS),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = covariate_names)
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="DSS"),
              orderer = ~order(p),
              labels_displayed = c("factor"),
              ggtheme = ggplot2::theme_bw(base_size = 13))
