library("data.table") # to load in specific columns from large dataset
library("R.utils")
library(tidyverse)


# Clear environment
rm(list=ls())

# Load in healthy xena data (gtex)
setwd("C:/Users/leahe/Documents/LIFE4137_independent_project/xena")
healthy_meta_data = read.table(file="healthy_samples/expected_counts/expected_counts_denseData.tsv", row.names=1, sep='\t', header=T)
healthy_liver_samples = healthy_meta_data[["samples"]]
class(healthy_liver_samples)

# fread used because the original file is too big
healthy_liver_data <- fread("healthy_samples/expected_counts/gtex_gene_expected_count.gz", 
                            select = c("sample", healthy_liver_samples)
)
warnings()
skipped_samples = c('GTEX-132NY-0926-SM-5P9G3', 
                    'GTEX-13N2G-0926-SM-5IFGJ',
                    'GTEX-13O1R-2026-SM-5KM3N',
                    'GTEX-145MF-0826-SM-5QGQA',
                    'GTEX-145MO-2326-SM-5NQ9K',
                    'GTEX-QESD-2026-SM-447BI',
                    'GTEX-QV44-0326-SM-4R1KD',
                    'GTEX-REY6-1226-SM-48FDR',
                    'GTEX-S4Z8-0526-SM-4AD4T',
                    'GTEX-T6MN-1226-SM-5CHSL',
                    'GTEX-TKQ2-1726-SM-4DXUP',
                    'GTEX-U8XE-1526-SM-5CHRK',
                    'GTEX-UPK5-1426-SM-4JBHH',
                    'GTEX-WFON-1726-SM-5CHT5',
                    'GTEX-WK11-1326-SM-4OOSI',
                    'GTEX-X3Y1-2726-SM-4PQZH',
                    'GTEX-YEC4-0826-SM-5P9FV',
                    'GTEX-11GSP-0626-SM-5986T',
                    'GTEX-12WSD-1426-SM-5GCN9',
                    'GTEX-12WSG-0626-SM-5FQTQ',
                    'GTEX-146FH-1526-SM-5NQBU',
                    'GTEX-14AS3-0126-SM-5Q5F4',
                    'GTEX-14DAQ-1726-SM-5S2R2',
                    'GTEX-PWOO-0826-SM-48TCL',
                    'GTEX-PX3G-0826-SM-48TZS',
                    'GTEX-Q734-0326-SM-48U15',
                    'GTEX-R53T-0326-SM-48FEC',
                    'GTEX-TSE9-1126-SM-5CHRH',
                    'GTEX-U3ZN-0226-SM-3DB8D',
                    'GTEX-W5X1-0126-SM-5CHR9',
                    'GTEX-X4EO-1126-SM-4QARQ'
)

# Save metadata subset for used samples
used_samples <- healthy_meta_data %>%
  filter(!samples %in% skipped_samples)
write.table(used_samples,
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/used_samples_metadata/gtex_healthy_metadata.txt",
            sep = "\t",
            row.names = FALSE,
)

# Count of each sex
table(used_samples$X_gender)

# Names of the healthy samples used (110 samples)
colnames(healthy_liver_data)

# Rename sample column to Ensembl_ID
names(healthy_liver_data)[names(healthy_liver_data) == "sample"] <- "Ensembl_ID"

# Remove gene version endings
healthy_liver_data$Ensembl_ID <- substr(healthy_liver_data$Ensembl_ID, 1, 15)

duplicated_healthy <- healthy_liver_data[duplicated(healthy_liver_data$Ensembl_ID),]
# No duplicates

# Check there are no duplicate samples/pseudoreplicates
healthy_samples <- colnames(healthy_liver_data)
healthy_samples
unique(healthy_samples)

TCGA_samples <- read.table(file="cancer_samples/raw_counts/TCGA-LIHC.star_counts.tsv.gz", sep='\t', header=T,  check.names = FALSE)
TCGA_cancer_meta <- read.table(file="cancer_samples/raw_counts/TCGA_allstages_all_primarydiag_noNA_metadata_348.tsv", sep='\t', header=T)

# Remove gene version endings
TCGA_samples$Ensembl_ID <- substr(TCGA_samples$Ensembl_ID, 1, 15)

# Find duplicate genes
TCGA_samples["duplicated"] <- duplicated(TCGA_samples$Ensembl_ID)
duplicated_TCGA <- TCGA_samples %>% select(Ensembl_ID, duplicated)

# 44 duplicated samples
duplicated_cancer <- TCGA_samples[duplicated(TCGA_samples$Ensembl_ID),]
unique(duplicated_cancer)
TCGA_samples <- TCGA_samples %>%
  distinct(Ensembl_ID, .keep_all = TRUE)

# Stage 1 Samples
s1_samples <- TCGA_cancer_meta %>%
  filter(ajcc_pathologic_stage.diagnoses == "Stage I")

# Stage 1
# Filter main data set for stage1 samples
stage1_sample_names = s1_samples[["samples"]]
stage1_sample_names

# Filter for cancer_samples for Stage1
stage1_cancer_samples <- TCGA_samples %>%
  select(all_of(stage1_sample_names))
# For some reason TCGA-DD-A1E9-01A`, `TCGA-DD-A3A0-01A`, `TCGA-DD-AADE-01A`, and `TCGA-DD-AAE8-01A` do not exist in cancer_samples
# Work around this by removing those samples from the selection
removed_s1 <- c("TCGA-DD-A1E9-01A", "TCGA-DD-A3A0-01A", "TCGA-DD-AADE-01A", "TCGA-DD-AAE8-01A")
allowed_stage1_sample_names <- stage1_sample_names[! stage1_sample_names %in% removed_s1]
# Filter for cancer_samples for Stage1
stage1_cancer_samples <- TCGA_samples %>%
  select("Ensembl_ID", all_of(allowed_stage1_sample_names)) 

# Stage1 sample names used in analysis (170 samples)
stage1_samples <- colnames(stage1_cancer_samples)
stage1_samples
# Check there are no duplicate samples/pseudoreplicates
unique(stage1_samples)


# Stage 2 samples
s2_samples <- TCGA_cancer_meta %>%
  filter(ajcc_pathologic_stage.diagnoses == "Stage II")

stage2_sample_names = s2_samples[["samples"]]
stage2_sample_names

# Filter for cancer_samples for Stage2
stage2_cancer_samples <- TCGA_samples %>%
  select(all_of(stage2_sample_names))
# The above errors due to columns not existing
# For some reason `TCGA.DD.AACM.01A` does not exist in cancer_samples. I think there is a formatting issue somewhere in the files
# Work around this by removing those samples from the selection
removed_s2 <- c("TCGA-DD-AACM-01A")
allowed_stage2_sample_names <- stage2_sample_names[! stage2_sample_names %in% removed_s2]
# Filter for cancer_samples for Stage2
stage2_cancer_samples <- TCGA_samples %>%
  select("Ensembl_ID", all_of(allowed_stage2_sample_names))

# Stage2 sample names used in analysis (84 samples)
stage2_samples <- colnames(stage2_cancer_samples)
stage2_samples
# Check there are no duplicate samples/pseudoreplicates
unique(stage2_samples)


# Stage 3 samples
s3_samples <- TCGA_cancer_meta %>%
  filter(ajcc_pathologic_stage.diagnoses == "Stage III" |
         ajcc_pathologic_stage.diagnoses == "Stage IIIA" |
         ajcc_pathologic_stage.diagnoses == "Stage IIIB" |
         ajcc_pathologic_stage.diagnoses == "Stage IIIC"
         )

stage3_sample_names = s3_samples[["samples"]]
stage3_sample_names

# Filter for cancer_samples for Stage3
stage3_cancer_samples <- TCGA_samples %>%
  select("Ensembl_ID", all_of(stage3_sample_names))
# For some reason `TCGA.G3.A25W.01A` does not exist in cancer_samples.
# Work around this by removing those samples from the selection, also remove -01B ending samples
removed_s3 <- c("TCGA-G3-A25W-01A", "TCGA-BC-4072-01B", "TCGA-BC-4073-01B")

allowed_stage3_sample_names <- stage3_sample_names[! stage3_sample_names %in% removed_s3]
# Filter for cancer_samples for Stage2
stage3_cancer_samples <- TCGA_samples %>%
  select("Ensembl_ID", all_of(allowed_stage3_sample_names))

# Stage3 sample names used in analysis (86 samples)
stage3_samples <- colnames(stage3_cancer_samples)
stage3_samples
# Check there are no duplicate samples/pseudoreplicates
unique(stage3_samples)
  

# Stage 4 samples
s4_samples <- TCGA_cancer_meta %>%
  filter(ajcc_pathologic_stage.diagnoses == "Stage IV" |
           ajcc_pathologic_stage.diagnoses == "Stage IVA" |
           ajcc_pathologic_stage.diagnoses == "Stage IVB"
           
  )

stage4_sample_names = s4_samples[["samples"]]
stage4_sample_names

# Filter for cancer_samples for Stage4
stage4_cancer_samples <- TCGA_samples %>%
  select("Ensembl_ID", all_of(stage4_sample_names))
# Stage4 sample names used in analysis (5 samples)
colnames(stage4_cancer_samples)


# Merge dataframes together
# Healthy vs stages
# Merge healthy with stage 1
healthy_and_stage1 <- merge(healthy_liver_data, stage1_cancer_samples, by = "Ensembl_ID", 
                            all.x = TRUE) # all.x will give NA for missing values


# Find rows with missing values to know which genes can't be compared
hs1_na_genes <- healthy_and_stage1[!complete.cases(healthy_and_stage1),]
# There are 3898 genes that do not have data for one of the two conditions

complete_healthy_and_stage1 <- na.omit(healthy_and_stage1)
# Only 56600 genes can be compared

# Save comparison count file with NO na values
write.csv(complete_healthy_and_stage1, file = "removed_versions/hs1/healthy_stage1_rawcounts.csv", row.names = FALSE)


# Merge healthy with stage2
healthy_and_stage2 <- merge(healthy_liver_data, stage2_cancer_samples, by = "Ensembl_ID", 
                            all.x = TRUE) # all.x will give NA for missing values

# Find rows with missing values to know which genes can't be compared
hs2_na_genes <- healthy_and_stage2[!complete.cases(healthy_and_stage2),]
# There are 3898 genes that do not have data for one of the two conditions

complete_healthy_and_stage2 <- na.omit(healthy_and_stage2)
# Only 56600 genes can be compared

# Save comparison count file with NO na values
write.csv(complete_healthy_and_stage2, file = "removed_versions/hs2/healthy_stage2_rawcounts.csv", row.names = FALSE)


# Merge healthy with stage 3
healthy_and_stage3 <- merge(healthy_liver_data, stage3_cancer_samples, by = "Ensembl_ID", 
                            all.x = TRUE) # all.x will give NA for missing values
# Find rows with missing values to know which genes can't be compared
hs3_na_genes <- healthy_and_stage3[!complete.cases(healthy_and_stage3),]
# There are 3898 genes that do not have data for one of the two conditions

complete_healthy_and_stage3 <- na.omit(healthy_and_stage3)
# Only 56600 genes can be compared

# Save comparison count file with NO na values
write.csv(complete_healthy_and_stage3, file = "removed_versions/hs3/healthy_stage3_rawcounts.csv", row.names = FALSE)


# Merge healthy with stage IV
healthy_and_stage4 <- merge(healthy_liver_data, stage4_cancer_samples, by = "Ensembl_ID", 
                            all.x = TRUE) # all.x will give NA for missing values
# Find rows with missing values to know which genes can't be compared
hs4_na_genes <- healthy_and_stage4[!complete.cases(healthy_and_stage4),]
# There are 3898 genes that do not have data for one of the two conditions

complete_healthy_and_stage4 <- na.omit(healthy_and_stage4)
# Only 56600 genes can be compared

# Save comparison count file with NO na values
write.csv(complete_healthy_and_stage4, file = "removed_versions/hs4/healthy_stage4_rawcounts.csv", row.names = FALSE)


# Stages vs stages
# Merge stage 1 with stage 2
stage1_and_stage2 <- merge(stage1_cancer_samples, stage2_cancer_samples, by = "Ensembl_ID", 
                           all.x = TRUE) # all.x will give NA for missing values
# Find rows with missing values to know which genes can't be compared
s1s2_na_genes <- stage1_and_stage2[!complete.cases(stage1_and_stage2),]
# All 60,660 genes can be compared between the two stages

# Save comparison count file
write.csv(stage1_and_stage2, file = "clean_script/s1s2/stage1_stage2_rawcounts.csv", row.names = FALSE)


# Merge stage 1 with stage 3
stage1_and_stage3 <- merge(stage1_cancer_samples, stage3_cancer_samples, by = "Ensembl_ID", 
                           all.x = TRUE) # all.x will give NA for missing values
# Find rows with missing values to know which genes can't be compared
s1s3_na_genes <- stage1_and_stage3[!complete.cases(stage1_and_stage3),]
# All 60,660 genes can be compared between the two stages

# Save comparison count file
write.csv(stage1_and_stage3, file = "clean_script/s1s3/stage1_stage3_rawcounts.csv", row.names = FALSE)


# Merge stageI and stage IV
stage1_and_stage4 <- merge(stage1_cancer_samples, stage4_cancer_samples, by = "Ensembl_ID", 
                           all.x = TRUE) # all.x will give NA for missing values
# Find rows with missing values to know which genes can't be compared
s1s4_na_genes <- stage1_and_stage4[!complete.cases(stage1_and_stage4),]
# All 60,660 genes can be compared between the two stages

# Save comparison count file
write.csv(stage1_and_stage4, file = "clean_script/s1s4/stage1_stage4_rawcounts.csv", row.names = FALSE)


# Merge all samples together for filtered cancer genes to match healthy samples
stages_list <- list(stage1_cancer_samples, stage2_cancer_samples, stage3_cancer_samples, stage4_cancer_samples)
merged_stages <- stages_list %>% reduce(full_join, by='Ensembl_ID')
# Merge healthy with all stages
healthy_and_stages <- merge(healthy_liver_data, merged_stages, by = "Ensembl_ID", 
                            all.x = TRUE) # all.x will give NA for missing values

# Find rows with missing values to know which genes can't be compared
healty_and_stages_na_genes <- healthy_and_stages[!complete.cases(healthy_and_stages),]
# There are 3898 genes that do not have data for one of the two conditions

complete_healthy_and_stages <- na.omit(healthy_and_stages)

# Save comparison count file with NO na values
write.csv(complete_healthy_and_stages, file = "removed_versions/healthy_vs_cancer/healthy_stages_rawcounts.csv", row.names = FALSE)


# Merge all stage samples together
complete_stages <- na.omit(merged_stages) # No genes are lost so all 60660 genes can be compared :)
# Save comparison count file for non-filtered TCGA stage samples
write.csv(complete_stages, file = "clean_script/s12_vs_s34/all_stages_rawcounts.csv", row.names = FALSE)
# Above can be used for s1_vs_s234 and stage_vs_stages comparisons


# Merge all stage samples together (Filtered)
s1234_filtered_genes <- filter(complete_stages, Ensembl_ID %in%  genes_in_hs14)  
write.csv(s1234_filtered_genes, file = "clean_script/s12_vs_s34/all_stages_rawcounts_filtered.csv", row.names = FALSE)
  

# Get metadata for the samples used (remove the samples that were not found in the TCGA_samples)
removed_samples <- append(removed_s1, removed_s2)
removed_samples <- append(removed_samples, removed_s3)
removed_samples

# Save metadata subset for used samples
used_TCGA_samples <- TCGA_cancer_meta %>%
  filter(!samples %in% removed_samples)
write.table(used_TCGA_samples,
            file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/used_samples_metadata/TCGA_metadata.txt",
            sep = "\t",
            row.names = FALSE,
)

# Count of each catergory
table(used_TCGA_samples$primary_diagnosis.diagnoses)
table(used_TCGA_samples$age_at_earliest_diagnosis_in_years.diagnoses.xena_derived) # ages are weird values
# Rename age column
names(used_TCGA_samples)[names(used_TCGA_samples) == "age_at_earliest_diagnosis_in_years.diagnoses.xena_derived"] <- "age_at_earliest_diagnosis"
used_TCGA_samples$age_at_earliest_diagnosis <- lapply(used_TCGA_samples$age_at_earliest_diagnosis, floor)

# Catergorise age groups and count them
used_TCGA_samples$age_range <- NA
used_TCGA_samples$age_range <- ifelse(used_TCGA_samples$age_at_earliest_diagnosis >= 80, "80+",
                                      ifelse(used_TCGA_samples$age_at_earliest_diagnosis >= 60, "60-79",
                                             ifelse(used_TCGA_samples$age_at_earliest_diagnosis >= 40, "40-59",
                                                    ifelse(used_TCGA_samples$age_at_earliest_diagnosis >= 18, "18-39", "0-17")
                                                           )
                                             )
                                      )
table(used_TCGA_samples$age_range)
table(used_TCGA_samples$ajcc_pathologic_stage.diagnoses)
table(used_TCGA_samples$ajcc_pathologic_m.diagnoses)
table(used_TCGA_samples$ajcc_pathologic_n.diagnoses)
table(used_TCGA_samples$gender.demographic)
table(used_TCGA_samples$ajcc_pathologic_t.diagnoses)
table(used_TCGA_samples$sample_type.samples)
