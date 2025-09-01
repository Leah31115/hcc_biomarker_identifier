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

# Uses TCGA, TARGET and GTEx combined cohourt TOIL data to look at isoform expression

# Clear environment
rm(list=ls())

setwd("C:/Users/leahe/Documents/LIFE4137_independent_project/xena")

# Using count data filtered from ada
# Filter TCGA tpm normalised counts for the paired samples
setwd("C:/Users/leahe/Documents/LIFE4137_independent_project/xena")
count_data <- read.table(file="splicing/isoform_counts.txt", sep='\t', header=T, check.names = FALSE)
cancer_meta <- read.table(file="used_samples_metadata/TCGA_metadata.txt", sep='\t', header=T)
healthy_meta <- read.table(file="used_samples_metadata/gtex_healthy_metadata.txt", sep='\t', header=T)
tcga_survival <- read.table(file="all_survival_LIHC_data.txt", sep='\t', header=T)

# Pair TCGA cancer and GTEx samples
cancer_sample_names = cancer_meta[["samples"]]
cancer_sample_names
# Need to remove the A from the end of the samples since the count data ends in 01
cancer_sample_names <- substr(cancer_sample_names, 1, 15)


healthy_sample_names = healthy_meta[["samples"]]
healthy_sample_names

test <- count_data[,1:2]

# Subset transcript counts for transcripts of interest
# Protein coding genes with upreg pseudogenes
RPL41 <- c("ENST00000617112.1",
           "ENST00000546591.5",
           "ENST00000501597.3",
           "ENST00000546485.5",
           "ENST00000358888.7",
           "ENST00000546654.1",
           "ENST00000552314.1",
           "ENST00000623936.1"
           )

GAPDH <- c("ENST00000229239.9",
           "ENST00000496049.1",
           "ENST00000396856.5",
           "ENST00000492719.5",
           "ENST00000396861.5",
           "ENST00000474249.5",
           "ENST00000466588.5",
           "ENST00000396859.5",
           "ENST00000466525.1",
           "ENST00000396858.5",
           "ENST00000619601.1"
           )

CHCHD2 <- c("ENST00000395422.3",
            "ENST00000473095.1"
            )



HMGN2 <- c("ENST00000468388.5",
           "ENST00000479815.5",
           "ENST00000460563.5",
           "ENST00000464888.5",
           "ENST00000463817.5",
           "ENST00000619352.4",
           "ENST00000467700.5",
           "ENST00000361427.5",
           "ENST00000466194.1",
           "ENST00000493418.1"
           )

# Protein coding genes

CIDEB <- c("ENST00000555817.1",
           "ENST00000258807.5",
           "ENST00000336557.9",
           "ENST00000555471.1",
           "ENST00000554411.5",
           "ENST00000556756.1"
           )


TOP2A <- c("ENST00000423485.5",
           "ENST00000581055.1",
           "ENST00000577706.1",
           "ENST00000578412.1",
           "ENST00000577541.1"
           )
  
CRHBP <- c("ENST00000274368.8",
           "ENST00000506501.1",
           "ENST00000512446.1",
           "ENST00000514258.1",
           "ENST00000503763.1"
           )

BIRC5 <- c("ENST00000301633.8",
           "ENST00000350051.7",
           "ENST00000590925.6",
           "ENST00000374948.6",
           "ENST00000590449.1",
           "ENST00000592734.5",
           "ENST00000592115.5",
           "ENST00000586192.5",
           "ENST00000587746.5",
           "ENST00000591800.1",
           "ENST00000589892.1"
           )



CDC20 <- c("ENST00000310955.10",
           "ENST00000372462.1",
           "ENST00000478882.1",
           "ENST00000482046.1"
           )

COL15A1 <- c("ENST00000467052.1",
             "ENST00000471477.1",
             "ENST00000375001.7",
             "ENST00000610452.1",
             "ENST00000496686.1"
             )

PTTG1 <- c("ENST00000352433.9",
           "ENST00000524244.5",
           "ENST00000517480.1",
           "ENST00000523659.5",
           "ENST00000520452.5",
           "ENST00000393964.1",
           "ENST00000519287.1"
)

UBE2C <- c("ENST00000352433.9",
           "ENST00000524244.5",
           "ENST00000517480.1",
           "ENST00000523659.5",
           "ENST00000520452.5",
           "ENST00000393964.1",
           "ENST00000519287.1"
           )

PITX <- c("ENST00000405520.5",
           "ENST00000617055.4",
           "ENST00000356455.8",
           "ENST00000335046.7",
           "ENST00000243893.10",
           "ENST00000352551.9",
           "ENST00000372568.4",
           "ENST00000496085.1"
           )

TMPRSS7 <- c("ENST00000452346.6",
             "ENST00000435737.5",
             "ENST00000419127.5",
             "ENST00000617607.4",
             "ENST00000460599.1"
             )

PNCK <- c("ENST00000458354.5",
          "ENST00000480693.1",
          "ENST00000484705.1",
          "ENST00000473831.2",
          "ENST00000419804.5",
          "ENST00000423545.1",
          "ENST00000425526.5",
          "ENST00000460106.1",
          "ENST00000488994.5",
          "ENST00000418241.5",
          "ENST00000462280.1",
          "ENST00000447676.6",
          "ENST00000488168.5",
          "ENST00000370142.5",
          "ENST00000466074.5",
          "ENST00000340888.7",
          "ENST00000489536.5",
          "ENST00000411968.5",
          "ENST00000466638.5",
          "ENST00000465303.5",
          "ENST00000472324.5",
          "ENST00000433470.5",
          "ENST00000475172.5",
          "ENST00000370150.5",
          "ENST00000463548.1",
          "ENST00000439087.5",
          "ENST00000434652.5",
          "ENST00000370145.8",
          "ENST00000422811.5",
          "ENST00000473680.1",
          "ENST00000466662.1",
          "ENST00000438984.5"
          )

PLA2G4D <- c("ENST00000290472.3",
             "ENST00000560132.1",
             "ENST00000560932.1"
             )

CA9 <- c("ENST00000378357.8",
         "ENST00000493245.1",
         "ENST00000485665.1")

WDR83OS <- c("ENST00000596731.5",
             "ENST00000600694.1",
             "ENST00000222190.9",
             "ENST00000598732.1"
)

SPINK13 <- c("ENST00000398450.4",
             "ENST00000512953.5",
             "ENST00000511106.5")

PTEN <- c("ENST00000371953.7",
          "ENST00000610634.1",
          "ENST00000487939.1",
          "ENST00000462694.1",
          "ENST00000618586.1",
          "ENST00000498703.1",
          "ENST00000472832.2"
          )

PPIA <- c("ENST00000481437.5",
          "ENST00000451562.5",
          "ENST00000468812.5",
          "ENST00000489459.5",
          "ENST00000355968.10",
          "ENST00000620047.4",
          "ENST00000479021.1",
          "ENST00000415933.5",
          "ENST00000494484.1",
          "ENST00000480603.1"
          )

E2F3 <- c("ENST00000613242.4",
          "ENST00000346618.7",
          "ENST00000535432.2"
          )


# goi = genes of interest
all_goi <- c(RPL41, GAPDH, CHCHD2, HMGN2, CIDEB, TOP2A, CRHBP, BIRC5, CDC20, COL15A1, WDR83OS,
             PTTG1, UBE2C, TMPRSS7, PNCK, PLA2G4D, CA9, SPINK13, PTEN, PPIA, E2F3)

# Filter for desired samples before mutating character gene counts to numeric for speed
count_data <- tibble::rownames_to_column(count_data, "Ensembl_ID")

count_data <- count_data %>%
  filter(Ensembl_ID %in% all_goi)

# Filter for healthy samples
healthy_samples <- count_data %>%
  select("Ensembl_ID", all_of(healthy_sample_names))
# 110 samples
rownames(healthy_samples) <- healthy_samples$Ensembl_ID
healthy_samples <- healthy_samples[,!(names(healthy_samples) %in% "Ensembl_ID")]

# Format healthy metadata (obtain sex only since this is shared with cancer samples)
rownames(healthy_meta) <- healthy_meta$samples
names(healthy_meta)[names(healthy_meta) == "X_gender"] <- "sex"
healthy_meta <- healthy_meta %>%
  select(sex)


# Filter for cancer_samples
cancer_samples <- count_data %>%
  select("Ensembl_ID", all_of(cancer_sample_names))

# TCGA-DD-A4NR-01` and `TCGA-DD-AACD-01`  do not exist
# Work around this by removing those samples from the selection
removed_cancer <- c("TCGA-DD-A4NR-01",
                    "TCGA-DD-AACD-01")

allowed_cancer_samples <- cancer_sample_names[! cancer_sample_names %in% removed_cancer]

# Format cancer metadata
# Trim off A suffix
cancer_meta$sample<- substring(cancer_meta$sample, 1, nchar(cancer_meta$sample)-1)
# Filter out unavailable samples
cancer_meta <- cancer_meta %>%
  filter(!sample %in% removed_cancer)
# Make rows the sample names
rownames(cancer_meta) <- cancer_meta$sample
names(cancer_meta)[names(cancer_meta) == "gender.demographic"] <- "sex"

# Filter for just sex
cancer_meta <- cancer_meta %>%
  select(sex)

# merge metadata
meta <- rbind(cancer_meta, healthy_meta)
write.csv(meta, file="C:/Users/leahe/Documents/LIFE4137_independent_project/xena/shared_meta.csv")


# Filter for cancer_samples
cancer_samples <- count_data %>%
  select("Ensembl_ID", all_of(allowed_cancer_samples)) 
# 338 samples!


# Format paired counts dataframe
# Make ensemblID column the row names
rownames(cancer_samples) <- cancer_samples$Ensembl_ID
# Remove ensemblID column
cancer_samples <- cancer_samples[,!(names(cancer_samples) %in% "Ensembl_ID")]


# Merge healthy and cancer samples
splice_df <- merge(healthy_samples, cancer_samples,
                   by = 'row.names', all = TRUE)

# Rename row.names column and make the rows the gene names
names(splice_df)[names(splice_df) == "Row.names"] <- "Ensembl_ID"
rownames(splice_df) <- splice_df$Ensembl_ID
splice_df <- splice_df[,!(names(splice_df) %in% "Ensembl_ID")]

splice_df <- tibble::rownames_to_column(splice_df, "Ensembl_ID")

# Antilog gene count function
antilog_counts <- function(tpm){
  antilog_count <- (2^tpm) - 0.001
  return(antilog_count)
}


# save splice variants
write.csv(splice_df, file = "C:/Users/leahe/Documents/LIFE4137_independent_project/xena/splicing/survival/transcripts_survival.csv") 
