
library("data.table") # to load in specific columns from large dataset
library("R.utils")
library(tidyverse)
library(dplyr) # for alt. columns 


# Clear environment
rm(list=ls())

setwd("C:/Users/leahe/Documents/LIFE4137_independent_project/xena")

TCGA_samples <- read.table(file="cancer_samples/raw_counts/TCGA-LIHC.star_counts.tsv.gz", sep='\t', header=T,   check.names = FALSE)
TCGA_cancer_meta <- read.table(file="cancer_samples/raw_counts/TCGA_allstages_all_primarydiag_noNA_metadata_348.tsv", sep='\t', header=T)

# TCGA Solid tissue normal (STN)
TCGA_STN_meta <- read.table(file="solid_tissue_normal/STN_metadata.tsv", sep='\t', header=T)

# Pair STN with TCGA cancer samples
cancer_sample_names = TCGA_cancer_meta[["samples"]]
cancer_sample_names

# Remove .01A end to get the main sample ID
trim_cancer_sample_names <- substr(cancer_sample_names, 1, nchar(cancer_sample_names) - 4)
trim_cancer_sample_names

STN_sample_names = TCGA_STN_meta[["samples"]]
STN_sample_names

is.vector(STN_sample_names)
# Remove .11A end to get the main sample ID
trim_STN_sample_names <- substr(STN_sample_names, 1, nchar(STN_sample_names) - 4)
trim_STN_sample_names

pairs <- trim_STN_sample_names %in% trim_cancer_sample_names
pairs
# All STN samples have a cancer sample pair so N=46

# Get count data with the pairs
# Filter for STN samples
STN_samples <- TCGA_samples %>%
  select(all_of(STN_sample_names))

# Some elements don't exist: `TCGA-ES-A2HS-11A`, `TCGA-G3-A25W-11A`, `TCGA-DD-A1E9-11A`, `TCGA-DD-A1ED-11A, TCGA-DD-A1EF-11A
STN_removed_samples <- c("TCGA-ES-A2HS-11A", "TCGA-G3-A25W-11A", "TCGA-G3-A25X-11A", "TCGA-DD-A1E9-11A", "TCGA-DD-A1ED-11A", "TCGA-DD-A1EF-11A", "TCGA-DD-A115-11A")

# Save metadata subset for used STN samples
used_STN_samples_meta <- TCGA_STN_meta %>%
  filter(!samples %in% STN_removed_samples)
write.table(used_STN_samples_meta,
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/used_samples_metadata/paired_STN_metadata.txt",
            sep = "\t",
            row.names = FALSE,
)

# Filter out removed elements
allowed_STN_samples <- STN_sample_names[! STN_sample_names %in% c("TCGA-ES-A2HS-11A", "TCGA-G3-A25W-11A", "TCGA-G3-A25X-11A", "TCGA-DD-A1E9-11A", "TCGA-DD-A1ED-11A", "TCGA-DD-A1EF-11A", "TCGA-DD-A115-11A")]
STN_samples <- TCGA_samples %>%
  select("Ensembl_ID", all_of(allowed_STN_samples)) # 39 samples

# Filter for cancer_samples
# Change end of the STN elements in the vector to .01A
allowed_cancer_samples <- substr(allowed_STN_samples, 1, nchar(allowed_STN_samples) - 4)
allowed_cancer_samples
allowed_cancer_samples <- paste(allowed_cancer_samples, "-01A", sep="")
allowed_cancer_samples
cancer_samples <- TCGA_samples %>%
  select("Ensembl_ID", all_of(allowed_cancer_samples))


# Merge the STN and cancer samples to get their pairs
# Alternating columns for paired stat test analysis
alt_paired_samples <- bind_cols(STN_samples, cancer_samples) %>% 
  dplyr::select(all_of(c(matrix(names(.), ncol = 40, byrow = TRUE))))
# remove first duplicated ensembl column
alt_paired_samples <- alt_paired_samples[,!(names(alt_paired_samples) %in% "Ensembl_ID...1")]
# Rename Ensembl column
names(alt_paired_samples)[names(alt_paired_samples) == "Ensembl_ID...41"] <- "Ensembl_ID"
# Save data
write.csv(alt_paired_samples, file = "paired_TCGA/alt_paired_rawcounts.csv", row.names = FALSE)


# HCC first then STN (control) samples merged
paired_samples <- merge(cancer_samples, STN_samples, by = "Ensembl_ID", 
                                     all.x = TRUE) # all.x will give NA for missing values

write.csv(paired_samples, file = "paired_TCGA/STN_cancer_pairs_rawcounts.csv", row.names = FALSE)


# Filter cancer meta data for cancer sample pair
# Convert samples column to a vector (makes a list)
used_STN_samples <- as.vector(used_STN_samples_meta["samples"])
used_STN_samples
# Trim -11A end from the sample names
used_pairs <- substr(used_STN_samples$samples, 1, nchar(used_STN_samples$samples) - 4)
used_pairs
# Replace sample name ending with -01A
used_cancer_samples <- paste(used_pairs, "-01A", sep="")
used_cancer_samples
# Filter and save meta data for the cancer samples
used_cancer_samples_meta <- TCGA_cancer_meta %>%
  filter(samples %in% used_cancer_samples)
write.table(used_cancer_samples_meta,
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/used_samples_metadata/paired_cancer_TCGA_metadata.txt",
            sep = "\t",
            row.names = FALSE,
)


# STN pair
# Count of each catergory
table(used_STN_samples_meta$primary_diagnosis.diagnoses)
table(used_STN_samples_meta$age_at_earliest_diagnosis_in_years.diagnoses.xena_derived) # ages are weird values
# Rename age column and round ages down to the year
names(used_STN_samples_meta)[names(used_STN_samples_meta) == "age_at_earliest_diagnosis_in_years.diagnoses.xena_derived"] <- "age_at_earliest_diagnosis"
used_STN_samples_meta$age_at_earliest_diagnosis <- lapply(used_STN_samples_meta$age_at_earliest_diagnosis, floor)

# Catergorise age groups and count them
used_STN_samples_meta$age_range <- NA
used_STN_samples_meta$age_range <- ifelse(used_STN_samples_meta$age_at_earliest_diagnosis >= 80, "80+",
                                      ifelse(used_STN_samples_meta$age_at_earliest_diagnosis >= 60, "60-79",
                                             ifelse(used_STN_samples_meta$age_at_earliest_diagnosis >= 40, "40-59",
                                                    ifelse(used_STN_samples_meta$age_at_earliest_diagnosis >= 18, "18-39", "0-17")
                                             )
                                      )
)
table(used_STN_samples_meta$age_range)
table(used_STN_samples_meta$ajcc_pathologic_stage.diagnoses)
table(used_STN_samples_meta$ajcc_pathologic_m.diagnoses)
table(used_STN_samples_meta$ajcc_pathologic_n.diagnoses)
table(used_STN_samples_meta$gender.demographic)
table(used_STN_samples_meta$ajcc_pathologic_t.diagnoses)
table(used_STN_samples_meta$sample_type.samples)

# Cancer pair
# Count of each catergory
table(used_cancer_samples_meta$primary_diagnosis.diagnoses)
table(used_cancer_samples_meta$age_at_earliest_diagnosis_in_years.diagnoses.xena_derived) # ages are weird values
# Rename age column and round ages down to the year
names(used_cancer_samples_meta)[names(used_cancer_samples_meta) == "age_at_earliest_diagnosis_in_years.diagnoses.xena_derived"] <- "age_at_earliest_diagnosis"
used_cancer_samples_meta$age_at_earliest_diagnosis <- lapply(used_cancer_samples_meta$age_at_earliest_diagnosis, floor)

# Catergorise age groups and count them
used_cancer_samples_meta$age_range <- NA
used_cancer_samples_meta$age_range <- ifelse(used_cancer_samples_meta$age_at_earliest_diagnosis >= 80, "80+",
                                      ifelse(used_cancer_samples_meta$age_at_earliest_diagnosis >= 60, "60-79",
                                             ifelse(used_cancer_samples_meta$age_at_earliest_diagnosis >= 40, "40-59",
                                                    ifelse(used_cancer_samples_meta$age_at_earliest_diagnosis >= 18, "18-39", "0-17")
                                             )
                                      )
)
table(used_cancer_samples_meta$age_range)
table(used_cancer_samples_meta$ajcc_pathologic_stage.diagnoses)
table(used_cancer_samples_meta$ajcc_pathologic_m.diagnoses)
table(used_cancer_samples_meta$ajcc_pathologic_n.diagnoses)
table(used_cancer_samples_meta$gender.demographic)
table(used_cancer_samples_meta$ajcc_pathologic_t.diagnoses)
table(used_cancer_samples_meta$sample_type.samples)
