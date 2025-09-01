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

RPL41_geneid <- c("ENSG00000256338.2",
           "ENSG00000256393.1",
           "ENSG00000227063.5",
           "ENSG00000229117.9",
           "ENSG00000264281.3"
)


# RPL41
RPL41_df <- TCGA_samples %>%
  filter(Ensembl_ID %in% c(RPL41_geneid))

# Rename row.names column and make the rows the gene names
rownames(RPL41_df) <- RPL41_df$Ensembl_ID
RPL41_df <- RPL41_df[,!(names(RPL41_df) %in% "Ensembl_ID")]
# Filter for used cancer samples
RPL41_df <- as.data.frame(t(RPL41_df))
RPL41_df <- tibble::rownames_to_column(RPL41_df, "sample")
RPL41_df <- RPL41_df %>%
  filter(sample %in% c(cancer_sample_names))
# Make sample the row names
rownames(RPL41_df) <- RPL41_df$sample
RPL41_df <- RPL41_df[,!(names(RPL41_df) %in% "sample")]


# Add meta
rownames(tcga_meta) <- tcga_meta$sample
tcga_meta <- tcga_meta[,!(names(tcga_meta) %in% c("sample", "samples"))]
RPL41_df <- merge(RPL41_df, tcga_meta,
                  by = 'row.names', all = TRUE)
# Remove A suffix to pair samples with survival data
names(RPL41_df)[names(RPL41_df) == "Row.names"] <- "sample"
RPL41_df$sample <- substring(RPL41_df$sample, 1, nchar(RPL41_df$sample)-1)
rownames(RPL41_df) <- RPL41_df$sample
RPL41_df <- RPL41_df[,!(names(RPL41_df) %in% "sample")]

colnames(RPL41_df)


names(RPL41_df)[names(RPL41_df) == "primary_diagnosis.diagnoses"] <- "primary_diagnosis"
names(RPL41_df)[names(RPL41_df) == "age_at_earliest_diagnosis_in_years.diagnoses.xena_derived"] <- "age_at_earliest_diagnosis"
names(RPL41_df)[names(RPL41_df) == "ajcc_pathologic_stage.diagnoses"] <- "all_stages"
names(RPL41_df)[names(RPL41_df) == "ajcc_pathologic_m.diagnoses"] <- "M"
names(RPL41_df)[names(RPL41_df) == "ajcc_pathologic_n.diagnoses"] <- "N"
names(RPL41_df)[names(RPL41_df) == "gender.demographic"] <- "sex"
names(RPL41_df)[names(RPL41_df) == "ajcc_pathologic_t.diagnoses"] <- "all_T"
colnames(RPL41_df)

# Round ages down so you just get their year
RPL41_df$age_at_earliest_diagnosis <- lapply(RPL41_df$age_at_earliest_diagnosis, floor)
# Assign age category
RPL41_df$age_range <- NA
RPL41_df$age_range <- ifelse(RPL41_df$age_at_earliest_diagnosis >= 80, "80+",
                                ifelse(RPL41_df$age_at_earliest_diagnosis >= 60, "60-79",
                                       ifelse(RPL41_df$age_at_earliest_diagnosis >= 40, "40-59",
                                              ifelse(RPL41_df$age_at_earliest_diagnosis >= 18, "18-39", "0-17")
                                       )
                                )
)


# Group subgroups in just one group e.g stage 4A = stage4
RPL41_df$stage <- RPL41_df$all_stages
RPL41_df$stage[RPL41_df$stage == "Stage IIIA" | RPL41_df$stage == "Stage IIIB" | RPL41_df$stage == "Stage IIIC"] <- "Stage III"
RPL41_df$stage[RPL41_df$stage == "Stage IVA" | RPL41_df$stage == "Stage IVB"] <- "Stage IV"
RPL41_df$T <- RPL41_df$all_T
RPL41_df$T[RPL41_df$T == "T2a" | RPL41_df$T == "T2b"] <- "T2"
RPL41_df$T[RPL41_df$T == "T3a" | RPL41_df$T == "T3b"] <- "T3"



# Add log prefix to transcript names
RPL41_df <- RPL41_df %>% 
  rename_with(function(x) paste0("log_", x), all_of(RPL41_geneid))
colnames(RPL41_df)


# Antilog gene count function
antilog_counts <- function(tpm){
  antilog_count <- (2^tpm) - 0.001
  return(antilog_count)
}

# Antilog
RPL41_df["RPL41P2"] <- lapply(RPL41_df["log_ENSG00000256338.2"], antilog_counts)
RPL41_df["RPL41P5"] <- lapply(RPL41_df["log_ENSG00000256393.1"], antilog_counts)
RPL41_df["RPL41P1"] <- lapply(RPL41_df["log_ENSG00000227063.5"], antilog_counts)
RPL41_df["RPL41"] <- lapply(RPL41_df["log_ENSG00000229117.9"], antilog_counts)
RPL41_df["RPL41_pseudogene"] <- lapply(RPL41_df["log_ENSG00000264281.3"], antilog_counts)

# subset for antilogged counts
RPL41_df_anti <- RPL41_df %>% 
  select(RPL41P2,
         RPL41P5,
         RPL41P1,
         RPL41,
         RPL41_pseudogene,
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
cancer_sample_names = row.names(RPL41_df)
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
RPL41_survival_df <- merge(RPL41_df_anti, filtered_tcga_survival, by = "row.names", 
                             all.x = TRUE) # all.x will give NA for missing values

# Convert days to months
# 1 day = 0.0328767 months
days_to_months_converter <- function(no_days){
  no_months <- round(no_days*0.0328767, 1)
  return(no_months)
}

RPL41_survival_df["OS_months"] <- lapply(RPL41_survival_df["OS.time"], days_to_months_converter)
RPL41_survival_df["DSS_months"] <- lapply(RPL41_survival_df["DSS.time"], days_to_months_converter)

# Investigating median gene expression category on survival
RPL41_survival_df$RPL41P2_median <- NA
RPL41_survival_df$RPL41P2_median[RPL41_survival_df$RPL41P2 <= quantile(RPL41_survival_df$RPL41P2, 0.5, na.rm=TRUE)] <- "LOW"
RPL41_survival_df$RPL41P2_median[RPL41_survival_df$RPL41P2 > quantile(RPL41_survival_df$RPL41P2, 0.5, na.rm=TRUE)] <- "HIGH"

RPL41_survival_df$RPL41P5_median <- NA
RPL41_survival_df$RPL41P5_median[RPL41_survival_df$RPL41P5 <= quantile(RPL41_survival_df$RPL41P5, 0.5, na.rm=TRUE)] <- "LOW"
RPL41_survival_df$RPL41P5_median[RPL41_survival_df$RPL41P5 > quantile(RPL41_survival_df$RPL41P5, 0.5, na.rm=TRUE)] <- "HIGH"

RPL41_survival_df$RPL41P1_median <- NA
RPL41_survival_df$RPL41P1_median[RPL41_survival_df$RPL41P1 <= quantile(RPL41_survival_df$RPL41P1, 0.5, na.rm=TRUE)] <- "LOW"
RPL41_survival_df$RPL41P1_median[RPL41_survival_df$RPL41P1 > quantile(RPL41_survival_df$RPL41P1, 0.5, na.rm=TRUE)] <- "HIGH"

RPL41_survival_df$RPL41_median <- NA
RPL41_survival_df$RPL41_median[RPL41_survival_df$RPL41 <= quantile(RPL41_survival_df$RPL41, 0.5, na.rm=TRUE)] <- "LOW"
RPL41_survival_df$RPL41_median[RPL41_survival_df$RPL41 > quantile(RPL41_survival_df$RPL41, 0.5, na.rm=TRUE)] <- "HIGH"

RPL41_survival_df$RPL41_pseudogene_median <- NA
RPL41_survival_df$RPL41_pseudogene_median[RPL41_survival_df$RPL41_pseudogene <= quantile(RPL41_survival_df$RPL41_pseudogene, 0.5, na.rm=TRUE)] <- "LOW"
RPL41_survival_df$RPL41_pseudogene_median[RPL41_survival_df$RPL41_pseudogene > quantile(RPL41_survival_df$RPL41_pseudogene, 0.5, na.rm=TRUE)] <- "HIGH"

covariate_names <- c(RPL41P5_median="RPL41P5 median",
                     RPL41P1_median="RPL41P1 median",
                     RPL41_median="RPL41 median",
                     RPL41P2_median="RPL41P2 median",
                     RPL41_pseudogene="RPL41 pseudogene median"
)


# OS
RPL41_median_OS <- RPL41_survival_df %>%
  select(Row.names, RPL41P5_median, RPL41P1_median, RPL41_median, RPL41P2_median, RPL41_pseudogene, sex, OS_months, OS)
# Remove NAs
RPL41_median_OS <- subset(RPL41_median_OS, !is.na(OS) & !is.na(OS_months))

map(vars(RPL41P5_median, RPL41P1_median, RPL41_median, RPL41P2_median, RPL41_pseudogene), function(by)
{
  analyse_multivariate(RPL41_median_OS,
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
RPL41_median_DSS <- RPL41_survival_df %>%
  select(Row.names, RPL41P5_median, RPL41P1_median, RPL41_median, RPL41P2_median, RPL41_pseudogene, sex, DSS_months, DSS)

# Remove NAs
RPL41_median_DSS <- subset(RPL41_median_DSS, !is.na(DSS) & !is.na(DSS_months))

map(vars(RPL41P5_median, RPL41P1_median, RPL41_median, RPL41P2_median, RPL41_pseudogene), function(by)
{
  analyse_multivariate(RPL41_median_DSS,
                       vars(DSS_months, DSS),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = covariate_names)
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="DSS"),
              orderer = ~order(p),
              labels_displayed = c("factor"),
              ggtheme = ggplot2::theme_bw(base_size = 14))

