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

RPS28_geneid <- c("ENSG00000233927.5",
                  "ENSG00000237039.1",
                  "ENSG00000227077.3"
)


# RPS28
RPS28_df <- TCGA_samples %>%
  filter(Ensembl_ID %in% c(RPS28_geneid))

# Rename row.names column and make the rows the gene names
rownames(RPS28_df) <- RPS28_df$Ensembl_ID
RPS28_df <- RPS28_df[,!(names(RPS28_df) %in% "Ensembl_ID")]
# Filter for used cancer samples
RPS28_df <- as.data.frame(t(RPS28_df))
RPS28_df <- tibble::rownames_to_column(RPS28_df, "sample")
RPS28_df <- RPS28_df %>%
  filter(sample %in% c(cancer_sample_names))
# Make sample the row names
rownames(RPS28_df) <- RPS28_df$sample
RPS28_df <- RPS28_df[,!(names(RPS28_df) %in% "sample")]


# Add meta
rownames(tcga_meta) <- tcga_meta$sample
tcga_meta <- tcga_meta[,!(names(tcga_meta) %in% c("sample", "samples"))]
RPS28_df <- merge(RPS28_df, tcga_meta,
                  by = 'row.names', all = TRUE)
# Remove A suffix to pair samples with survival data
names(RPS28_df)[names(RPS28_df) == "Row.names"] <- "sample"
RPS28_df$sample <- substring(RPS28_df$sample, 1, nchar(RPS28_df$sample)-1)
rownames(RPS28_df) <- RPS28_df$sample
RPS28_df <- RPS28_df[,!(names(RPS28_df) %in% "sample")]

colnames(RPS28_df)


names(RPS28_df)[names(RPS28_df) == "primary_diagnosis.diagnoses"] <- "primary_diagnosis"
names(RPS28_df)[names(RPS28_df) == "age_at_earliest_diagnosis_in_years.diagnoses.xena_derived"] <- "age_at_earliest_diagnosis"
names(RPS28_df)[names(RPS28_df) == "ajcc_pathologic_stage.diagnoses"] <- "all_stages"
names(RPS28_df)[names(RPS28_df) == "ajcc_pathologic_m.diagnoses"] <- "M"
names(RPS28_df)[names(RPS28_df) == "ajcc_pathologic_n.diagnoses"] <- "N"
names(RPS28_df)[names(RPS28_df) == "gender.demographic"] <- "sex"
names(RPS28_df)[names(RPS28_df) == "ajcc_pathologic_t.diagnoses"] <- "all_T"
colnames(RPS28_df)

# Round ages down so you just get their year
RPS28_df$age_at_earliest_diagnosis <- lapply(RPS28_df$age_at_earliest_diagnosis, floor)
# Assign age category
RPS28_df$age_range <- NA
RPS28_df$age_range <- ifelse(RPS28_df$age_at_earliest_diagnosis >= 80, "80+",
                                ifelse(RPS28_df$age_at_earliest_diagnosis >= 60, "60-79",
                                       ifelse(RPS28_df$age_at_earliest_diagnosis >= 40, "40-59",
                                              ifelse(RPS28_df$age_at_earliest_diagnosis >= 18, "18-39", "0-17")
                                       )
                                )
)


# Group subgroups in just one group e.g stage 4A = stage4
RPS28_df$stage <- RPS28_df$all_stages
RPS28_df$stage[RPS28_df$stage == "Stage IIIA" | RPS28_df$stage == "Stage IIIB" | RPS28_df$stage == "Stage IIIC"] <- "Stage III"
RPS28_df$stage[RPS28_df$stage == "Stage IVA" | RPS28_df$stage == "Stage IVB"] <- "Stage IV"
RPS28_df$T <- RPS28_df$all_T
RPS28_df$T[RPS28_df$T == "T2a" | RPS28_df$T == "T2b"] <- "T2"
RPS28_df$T[RPS28_df$T == "T3a" | RPS28_df$T == "T3b"] <- "T3"


# Add log prefix to transcript names
RPS28_df <- RPS28_df %>% 
  rename_with(function(x) paste0("log_", x), all_of(RPS28_geneid))
colnames(RPS28_df)


# Antilog gene count function
antilog_counts <- function(tpm){
  antilog_count <- (2^tpm) - 0.001
  return(antilog_count)
}

# Antilog
RPS28_df["RPS28"] <- lapply(RPS28_df["log_ENSG00000233927.5"], antilog_counts)
RPS28_df["RPS28P4"] <- lapply(RPS28_df["log_ENSG00000237039.1"], antilog_counts)
RPS28_df["RPS28_pseudogene"] <- lapply(RPS28_df["log_ENSG00000227077.3"], antilog_counts)

# subset for antilogged counts
RPS28_df_anti <- RPS28_df %>% 
  select(RPS28,
         RPS28P4,
         RPS28_pseudogene,
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
cancer_sample_names = row.names(RPS28_df)
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
RPS28_survival_df <- merge(RPS28_df_anti, filtered_tcga_survival, by = "row.names", 
                             all.x = TRUE) # all.x will give NA for missing values

# Convert days to months
# 1 day = 0.0328767 months
days_to_months_converter <- function(no_days){
  no_months <- round(no_days*0.0328767, 1)
  return(no_months)
}

RPS28_survival_df["OS_months"] <- lapply(RPS28_survival_df["OS.time"], days_to_months_converter)
RPS28_survival_df["DSS_months"] <- lapply(RPS28_survival_df["DSS.time"], days_to_months_converter)

# Investigating median gene expression category on survival
RPS28_survival_df$RPS28_median <- NA
RPS28_survival_df$RPS28_median[RPS28_survival_df$RPS28 <= quantile(RPS28_survival_df$RPS28, 0.5, na.rm=TRUE)] <- "LOW"
RPS28_survival_df$RPS28_median[RPS28_survival_df$RPS28 > quantile(RPS28_survival_df$RPS28, 0.5, na.rm=TRUE)] <- "HIGH"

RPS28_survival_df$RPS28P4_median <- NA
RPS28_survival_df$RPS28P4_median[RPS28_survival_df$RPS28P4 <= quantile(RPS28_survival_df$RPS28P4, 0.5, na.rm=TRUE)] <- "LOW"
RPS28_survival_df$RPS28P4_median[RPS28_survival_df$RPS28P4 > quantile(RPS28_survival_df$RPS28P4, 0.5, na.rm=TRUE)] <- "HIGH"

RPS28_survival_df$RPS28_pseudogene_median <- NA
RPS28_survival_df$RPS28_pseudogene_median[RPS28_survival_df$RPS28_pseudogene <= quantile(RPS28_survival_df$RPS28_pseudogene, 0.5, na.rm=TRUE)] <- "LOW"
RPS28_survival_df$RPS28_pseudogene_median[RPS28_survival_df$RPS28_pseudogene > quantile(RPS28_survival_df$RPS28_pseudogene, 0.5, na.rm=TRUE)] <- "HIGH"

covariate_names <- c(RPS28P4_median="RPS28P4 median",
                     RPS28_pseudogene_median="RPS28_pseudogene median",
                     RPS28_median="RPS28 median"
)

# OS
RPS28_median_OS <- RPS28_survival_df %>%
  select(Row.names, RPS28P4_median, RPS28_pseudogene_median, RPS28_median, sex, OS_months, OS)
# Remove NAs
RPS28_median_OS <- subset(RPS28_median_OS, !is.na(OS) & !is.na(OS_months))

map(vars(RPS28P4_median, RPS28_pseudogene_median, RPS28_median), function(by)
{
  analyse_multivariate(RPS28_median_OS,
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
RPS28_median_DSS <- RPS28_survival_df %>%
  select(Row.names, RPS28P4_median, RPS28_pseudogene_median, RPS28_median, sex, DSS_months, DSS)

# Remove NAs
RPS28_median_DSS <- subset(RPS28_median_DSS, !is.na(DSS) & !is.na(DSS_months))

map(vars(RPS28P4_median, RPS28_pseudogene_median, RPS28_median), function(by)
{
  analyse_multivariate(RPS28_median_DSS,
                       vars(DSS_months, DSS),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = covariate_names)
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="DSS"),
              orderer = ~order(p),
              labels_displayed = c("factor"),
              ggtheme = ggplot2::theme_bw(base_size = 14))
