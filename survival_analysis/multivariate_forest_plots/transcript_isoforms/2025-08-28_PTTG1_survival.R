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
splice_df <- read.csv(file="splicing/survival/transcripts_survival.csv", row.names=1, check.names = FALSE, header=T)
meta <- read.csv(file="shared_meta.csv", row.names=1, check.names = FALSE, header=T)
tcga_survival <- read.table(file="all_survival_LIHC_data.txt", sep='\t', header=T)

PTTG1 <- c("ENST00000352433.9",
           "ENST00000524244.5",
           "ENST00000517480.1",
           "ENST00000523659.5",
           "ENST00000520452.5",
           "ENST00000393964.1",
           "ENST00000519287.1"
)

# PTTG1
PTTG1_df <- splice_df %>%
  filter(Ensembl_ID %in% c(PTTG1))

# Rename row.names column and make the rows the gene names
rownames(PTTG1_df) <- PTTG1_df$Ensembl_ID
PTTG1_df <- PTTG1_df[,!(names(PTTG1_df) %in% "Ensembl_ID")]

# Make columns (samples) the row name
PTTG1_df <- as.data.frame(t(PTTG1_df))
# Assign sample category
PTTG1_df$sample_type <- NA
name_start <- substr(rownames(PTTG1_df), 1, 4)
PTTG1_df$sample_type <- ifelse(name_start == "GTEX", "Healthy", "Cancer")

# Add meta
PTTG1_df <- merge(PTTG1_df, meta,
                  by = 'row.names', all = TRUE)
names(PTTG1_df)[names(PTTG1_df) == "Row.names"] <- "samples"
rownames(PTTG1_df) <- PTTG1_df$samples
PTTG1_df <- PTTG1_df[,!(names(PTTG1_df) %in% "samples")]

# Change transcript IDs to transcript names
# Change transcript IDs to transcript names
colnames(PTTG1_df)[which(names(PTTG1_df) == "ENST00000352433.9")] <- "PTTG1.201"
colnames(PTTG1_df)[which(names(PTTG1_df) == "ENST00000393964.1")] <- "PTTG1.202"
colnames(PTTG1_df)[which(names(PTTG1_df) == "ENST00000520452.5")] <- "PTTG1.205"
colnames(PTTG1_df)[which(names(PTTG1_df) == "ENST00000517480.1")] <- "PTTG1.203"
colnames(PTTG1_df)[which(names(PTTG1_df) == "ENST00000519287.1")] <- "PTTG1.204"
colnames(PTTG1_df)[which(names(PTTG1_df) == "ENST00000523659.5")] <- "PTTG1.206"
colnames(PTTG1_df)[which(names(PTTG1_df) == "ENST00000524244.5")] <- "PTTG1.207"

# Add log prefix to transcript names
PTTG1_df <- PTTG1_df %>% 
  rename_with(function(x) paste0("log_", x), -c(sex, sample_type))
colnames(PTTG1_df)


# Antilog gene count function
antilog_counts <- function(tpm){
  antilog_count <- (2^tpm) - 0.001
  return(antilog_count)
}

# Antilog
PTTG1_df["PTTG1.201"] <- lapply(PTTG1_df["log_PTTG1.201"], antilog_counts)
PTTG1_df["PTTG1.207"] <- lapply(PTTG1_df["log_PTTG1.207"], antilog_counts)
PTTG1_df["PTTG1.203"] <- lapply(PTTG1_df["log_PTTG1.203"], antilog_counts)
PTTG1_df["PTTG1.206"] <- lapply(PTTG1_df["log_PTTG1.206"], antilog_counts)
PTTG1_df["PTTG1.205"] <- lapply(PTTG1_df["log_PTTG1.205"], antilog_counts)
PTTG1_df["PTTG1.202"] <- lapply(PTTG1_df["log_PTTG1.202"], antilog_counts)
PTTG1_df["PTTG1.204"] <- lapply(PTTG1_df["log_PTTG1.204"], antilog_counts)

# subset for antilogged counts
PTTG1_df_anti <- PTTG1_df %>% 
  select(PTTG1.201,
         PTTG1.207,
         PTTG1.203,
         PTTG1.206,
         PTTG1.205,
         PTTG1.202,
         PTTG1.204,
         sample_type,
         sex)


# Melt df
long_PTTG1_anti <- melt(PTTG1_df_anti, id_vars="sample_type")

# Problem with negative antilog counts and logging the y-axis
ggplot(long_PTTG1_anti, aes(x = variable, y = value, col=sample_type)) +
  geom_point(alpha=0.1, position="jitter") +
  geom_boxplot() +
  scale_y_continuous(trans = "log2") +
  xlab("PTTG1 splice variants") +
  ylab("Anti log TPM gene counts") +
  theme_bw()

# Subset for logged counts
log_PTTG1_df <- PTTG1_df %>% 
  select(log_PTTG1.201,
         log_PTTG1.207,
         log_PTTG1.203,
         log_PTTG1.206,
         log_PTTG1.205,
         log_PTTG1.202,
         log_PTTG1.204,
         sample_type,
         sex)

# Remove log prefix for graph formatting
names(log_PTTG1_df) <- sub('log_', '', names(log_PTTG1_df))


long_PTTG1_log <- melt(log_PTTG1_df, id_vars="sample_type")

# transcripts on y axis due to ensembl long names
ggplot(long_PTTG1_log, aes(x = value, y = variable, col=sample_type)) +
  geom_boxplot() +
  geom_point(alpha=0.09) +
  xlab("Log TPM gene counts") +
  ylab("PTTG1 splice variants") +
  theme_bw(base_size = 14)


# Check count distributions
# PTTG1 protein coding gene
PTTG1_distribution <- PTTG1_df %>% 
  select(log_PTTG1.201, sample_type)
# Summary stats
group_by(PTTG1_distribution, sample_type) %>%
  summarise(
    count = n(),
    median = median(log_PTTG1.201, na.rm = TRUE),
    IQR = IQR(log_PTTG1.201, na.rm = TRUE)
  )
# Healthy distribution
healthy_PTTG1 <- PTTG1_distribution %>%
  filter(sample_type == "Healthy")

ggdensity(healthy_PTTG1$log_PTTG1.201, 
          main = "Density plot of PTTG1 gene from healthy samples",
          xlab = "Tpm normalised gene count")

shapiro.test(healthy_PTTG1$log_PTTG1.201) 
# significant so non-parametric distribution

# Identifying extreme outliers
# Find Q1, Q3, and interquartile range for values
Q1 <- quantile(healthy_PTTG1$log_PTTG1.201, .25)
Q3 <- quantile(healthy_PTTG1$log_PTTG1.201, .75)
IQR <- IQR(healthy_PTTG1$log_PTTG1.201)

# subset data where points value is outside 2*1.5*IQR of Q1 and Q3
outliers <- subset(healthy_PTTG1, healthy_PTTG1$log_PTTG1.201<(Q1 - 2*1.5*IQR) | healthy_PTTG1$log_PTTG1.201>(Q3 + 2*1.5*IQR))
# 0 extreme outliers

# cancer distribution
cancer_PTTG1_t <- PTTG1_distribution %>%
  filter(sample_type == "Cancer")

ggdensity(cancer_PTTG1_t$log_PTTG1.201, 
          main = "Density plot of PTTG1 gene from cancer samples",
          xlab = "Tpm normalised gene count")

shapiro.test(cancer_PTTG1_t$log_PTTG1.201) 
# significant so non-parametric distribution

# Identifying extreme outliers
# Find Q1, Q3, and interquartile range for values
Q1 <- quantile(cancer_PTTG1_t$log_PTTG1.201, .25)
Q3 <- quantile(cancer_PTTG1_t$log_PTTG1.201, .75)
IQR <- IQR(cancer_PTTG1_t$log_PTTG1.201)

# subset data where points value is outside 2*1.5*IQR of Q1 and Q3
outliers <- subset(cancer_PTTG1_t, cancer_PTTG1_t$log_PTTG1.201<(Q1 - 2*1.5*IQR) | cancer_PTTG1_t$log_PTTG1.201>(Q3 + 2*1.5*IQR))
# 0 extreme outliers

# Generate pvalues
stat.test <- long_PTTG1_log %>%
  group_by(variable) %>%
  wilcox_test(value ~ sample_type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

# Plot
plot <- ggplot(long_PTTG1_log, aes(x = variable, y = value, col=sample_type)) +
  geom_boxplot() +
  geom_point(alpha = 0.2) +
  xlab("PTTG1 transcript isoform") +
  ylab("Log TPM gene counts") +
  theme_bw(base_size = 16) +
  theme( legend.box.margin=margin(,,,-13),
         legend.text = element_text(size=11), 
         legend.title=element_text(size=13)
  ) 
plot

# Add p-values onto the box plots
stat_test <- stat.test %>%
  add_x_position(x="variable")


plot + coord_flip() + stat_pvalue_manual(
  stat_test, tip.length = 0, y.position=16,
  symnum.args = list(cutpoints = c(0, 0.05, 1), symbols = c("*", "ns")),
  label.size = 5
)


# Transcript survival analysis
# Filter survival data for used samples
cancer_sample_names = row.names(cancer_PTTG1_t)
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

cancer_PTTG1 <- PTTG1_df_anti %>%
  filter(sample_type == "Cancer")

# Merge survival data with gene counts 
PTTG1_survival_df <- merge(cancer_PTTG1, filtered_tcga_survival, by = "row.names", 
                           all.x = TRUE) # all.x will give NA for missing values


# Convert days to months
# 1 day = 0.0328767 months
days_to_months_converter <- function(no_days){
  no_months <- round(no_days*0.0328767, 1)
  return(no_months)
}

PTTG1_survival_df["OS_months"] <- lapply(PTTG1_survival_df["OS.time"], days_to_months_converter)
PTTG1_survival_df["DSS_months"] <- lapply(PTTG1_survival_df["DSS.time"], days_to_months_converter)

# Investigating median gene expression category on survival
PTTG1_survival_df$PTTG1.201_median <- NA
PTTG1_survival_df$PTTG1.201_median[PTTG1_survival_df$PTTG1.201 <= quantile(PTTG1_survival_df$PTTG1.201, 0.5, na.rm=TRUE)] <- "LOW"
PTTG1_survival_df$PTTG1.201_median[PTTG1_survival_df$PTTG1.201 > quantile(PTTG1_survival_df$PTTG1.201, 0.5, na.rm=TRUE)] <- "HIGH"

PTTG1_survival_df$PTTG1.207_median <- NA
PTTG1_survival_df$PTTG1.207_median[PTTG1_survival_df$PTTG1.207 <= quantile(PTTG1_survival_df$PTTG1.207, 0.5, na.rm=TRUE)] <- "LOW"
PTTG1_survival_df$PTTG1.207_median[PTTG1_survival_df$PTTG1.207 > quantile(PTTG1_survival_df$PTTG1.207, 0.5, na.rm=TRUE)] <- "HIGH"

PTTG1_survival_df$PTTG1.203_median <- NA
PTTG1_survival_df$PTTG1.203_median[PTTG1_survival_df$PTTG1.203 <= quantile(PTTG1_survival_df$PTTG1.203, 0.5, na.rm=TRUE)] <- "LOW"
PTTG1_survival_df$PTTG1.203_median[PTTG1_survival_df$PTTG1.203 > quantile(PTTG1_survival_df$PTTG1.203, 0.5, na.rm=TRUE)] <- "HIGH"

PTTG1_survival_df$PTTG1.206_median <- NA
PTTG1_survival_df$PTTG1.206_median[PTTG1_survival_df$PTTG1.206 <= quantile(PTTG1_survival_df$PTTG1.206, 0.5, na.rm=TRUE)] <- "LOW"
PTTG1_survival_df$PTTG1.206_median[PTTG1_survival_df$PTTG1.206 > quantile(PTTG1_survival_df$PTTG1.206, 0.5, na.rm=TRUE)] <- "HIGH"

PTTG1_survival_df$PTTG1.205_median <- NA
PTTG1_survival_df$PTTG1.205_median[PTTG1_survival_df$PTTG1.205 <= quantile(PTTG1_survival_df$PTTG1.205, 0.5, na.rm=TRUE)] <- "LOW"
PTTG1_survival_df$PTTG1.205_median[PTTG1_survival_df$PTTG1.205 > quantile(PTTG1_survival_df$PTTG1.205, 0.5, na.rm=TRUE)] <- "HIGH"

PTTG1_survival_df$PTTG1.202_median <- NA
PTTG1_survival_df$PTTG1.202_median[PTTG1_survival_df$PTTG1.202 <= quantile(PTTG1_survival_df$PTTG1.202, 0.5, na.rm=TRUE)] <- "LOW"
PTTG1_survival_df$PTTG1.202_median[PTTG1_survival_df$PTTG1.202 > quantile(PTTG1_survival_df$PTTG1.202, 0.5, na.rm=TRUE)] <- "HIGH"

PTTG1_survival_df$PTTG1.204_median <- NA
PTTG1_survival_df$PTTG1.204_median[PTTG1_survival_df$PTTG1.204 <= quantile(PTTG1_survival_df$PTTG1.204, 0.5, na.rm=TRUE)] <- "LOW"
PTTG1_survival_df$PTTG1.204_median[PTTG1_survival_df$PTTG1.204 > quantile(PTTG1_survival_df$PTTG1.204, 0.5, na.rm=TRUE)] <- "HIGH"

covariate_names <- c(PTTG1.207_median="PTTG1.207 median",
                     PTTG1.203_median="PTTG1.203 median",
                     PTTG1.206_median="PTTG1.206 median",
                     PTTG1.201_median="PTTG1.201 median",
                     PTTG1.205_median="PTTG1.205 median",
                     PTTG1.202_median="PTTG1.202 median",
                     PTTG1.204_median="PTTG1.204 median"
)


# OS
PTTG1_median_OS <- PTTG1_survival_df %>%
  select(Row.names, PTTG1.207_median, PTTG1.203_median, PTTG1.206_median, PTTG1.201_median, PTTG1.205_median, PTTG1.202_median, PTTG1.204_median, sex, OS_months, OS)
# Remove NAs
PTTG1_median_OS <- subset(PTTG1_median_OS, !is.na(OS) & !is.na(OS_months))

map(vars(PTTG1.207_median, PTTG1.203_median, PTTG1.206_median, PTTG1.201_median, PTTG1.205_median, PTTG1.202_median, PTTG1.204_median), function(by)
{
  analyse_multivariate(PTTG1_median_OS,
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
PTTG1_median_DSS <- PTTG1_survival_df %>%
  select(Row.names, PTTG1.207_median, PTTG1.203_median, PTTG1.206_median, PTTG1.201_median, PTTG1.205_median, PTTG1.202_median, PTTG1.204_median, sex, DSS_months, DSS)

# Remove NAs
PTTG1_median_DSS <- subset(PTTG1_median_DSS, !is.na(DSS) & !is.na(DSS_months))

map(vars(PTTG1.207_median, PTTG1.203_median, PTTG1.206_median, PTTG1.201_median, PTTG1.205_median, PTTG1.202_median, PTTG1.204_median), function(by)
{
  analyse_multivariate(PTTG1_median_DSS,
                       vars(DSS_months, DSS),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = covariate_names)
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="DSS"),
              orderer = ~order(p),
              labels_displayed = c("factor"),
              ggtheme = ggplot2::theme_bw(base_size = 14))
