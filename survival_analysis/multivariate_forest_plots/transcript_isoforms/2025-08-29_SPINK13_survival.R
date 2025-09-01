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

SPINK13 <- c("ENST00000398450.4",
             "ENST00000512953.5",
             "ENST00000511106.5")
# SPINK13
SPINK13_df <- splice_df %>%
  filter(Ensembl_ID %in% c(SPINK13))

# Rename row.names column and make the rows the gene names
rownames(SPINK13_df) <- SPINK13_df$Ensembl_ID
SPINK13_df <- SPINK13_df[,!(names(SPINK13_df) %in% "Ensembl_ID")]

# Make columns (samples) the row name
SPINK13_df <- as.data.frame(t(SPINK13_df))
# Assign sample category
SPINK13_df$sample_type <- NA
name_start <- substr(rownames(SPINK13_df), 1, 4)
SPINK13_df$sample_type <- ifelse(name_start == "GTEX", "Healthy", "Cancer")

# Add meta
SPINK13_df <- merge(SPINK13_df, meta,
                  by = 'row.names', all = TRUE)
names(SPINK13_df)[names(SPINK13_df) == "Row.names"] <- "samples"
rownames(SPINK13_df) <- SPINK13_df$samples
SPINK13_df <- SPINK13_df[,!(names(SPINK13_df) %in% "samples")]

# Change transcript IDs to transcript names
colnames(SPINK13_df)[which(names(SPINK13_df) == "ENST00000398450.4")] <- "SPINK13.201"
colnames(SPINK13_df)[which(names(SPINK13_df) == "ENST00000512953.5")] <- "SPINK13.203"
colnames(SPINK13_df)[which(names(SPINK13_df) == "ENST00000511106.5")] <- "SPINK13.202"

# Add log prefix to transcript names
SPINK13_df <- SPINK13_df %>% 
  rename_with(function(x) paste0("log_", x), -c(sex, sample_type))
colnames(SPINK13_df)


# Antilog gene count function
antilog_counts <- function(tpm){
  antilog_count <- (2^tpm) - 0.001
  return(antilog_count)
}

# Antilog
SPINK13_df["SPINK13.201"] <- lapply(SPINK13_df["log_SPINK13.201"], antilog_counts)
SPINK13_df["SPINK13.202"] <- lapply(SPINK13_df["log_SPINK13.202"], antilog_counts)
SPINK13_df["SPINK13.203"] <- lapply(SPINK13_df["log_SPINK13.203"], antilog_counts)

# subset for antilogged counts
SPINK13_df_anti <- SPINK13_df %>% 
  select(SPINK13.201,
         SPINK13.202,
         SPINK13.203,
         sample_type,
         sex)

# Melt df
long_SPINK13_anti <- melt(SPINK13_df_anti, id_vars="sample_type")

# Problem with negative antilog counts and logging the y-axis
ggplot(long_SPINK13_anti, aes(x = variable, y = value, col=sample_type)) +
  geom_point(alpha=0.1, position="jitter") +
  geom_boxplot() +
  scale_y_continuous(trans = "log2") +
  xlab("SPINK13 splice variants") +
  ylab("Anti log TPM gene counts") +
  theme_bw()

# Subset for logged counts
log_SPINK13_df <- SPINK13_df %>% 
  select(log_SPINK13.201,
         log_SPINK13.202,
         log_SPINK13.203,
         sample_type,
         sex)

# Remove log prefix for graph formatting
names(log_SPINK13_df) <- sub('log_', '', names(log_SPINK13_df))

long_SPINK13_log <- melt(log_SPINK13_df, id_vars="sample_type")

# transcripts on y axis due to ensembl long names
ggplot(long_SPINK13_log, aes(x = value, y = variable, col=sample_type)) +
  geom_boxplot() +
  geom_point(alpha=0.09) +
  xlab("Log TPM gene counts") +
  ylab("SPINK13 transcript variants") +
  theme_bw(base_size = 16)


# Check count distributions
# SPINK13 protein coding gene
SPINK13_distribution <- SPINK13_df %>% 
  select(log_SPINK13.201, sample_type)
# Summary stats
group_by(SPINK13_distribution, sample_type) %>%
  summarise(
    count = n(),
    median = median(log_SPINK13.201, na.rm = TRUE),
    IQR = IQR(log_SPINK13.201, na.rm = TRUE)
  )
# Healthy distribution
healthy_SPINK13 <- SPINK13_distribution %>%
  filter(sample_type == "Healthy")

ggdensity(healthy_SPINK13$log_SPINK13.201, 
          main = "Density plot of SPINK13 gene from healthy samples",
          xlab = "Tpm normalised gene count")

shapiro.test(healthy_SPINK13$log_SPINK13.201) 
# significant so non-parametric distribution

# Identifying extreme outliers
# Find Q1, Q3, and interquartile range for values
Q1 <- quantile(healthy_SPINK13$log_SPINK13.201, .25)
Q3 <- quantile(healthy_SPINK13$log_SPINK13.201, .75)
IQR <- IQR(healthy_SPINK13$log_SPINK13.201)

# subset data where points value is outside 2*1.5*IQR of Q1 and Q3
outliers <- subset(healthy_SPINK13, healthy_SPINK13$log_SPINK13.201<(Q1 - 2*1.5*IQR) | healthy_SPINK13$log_SPINK13.201>(Q3 + 2*1.5*IQR))
# 3 extreme outliers

# cancer distribution
cancer_SPINK13_t <- SPINK13_distribution %>%
  filter(sample_type == "Cancer")

ggdensity(cancer_SPINK13_t$log_SPINK13.201, 
          main = "Density plot of SPINK13 gene from cancer samples",
          xlab = "Tpm normalised gene count")

shapiro.test(cancer_SPINK13_t$log_SPINK13.201) 
# Significant so parametric distribution

# Identifying extreme outliers
# Find Q1, Q3, and interquartile range for values
Q1 <- quantile(cancer_SPINK13_t$log_SPINK13.201, .25)
Q3 <- quantile(cancer_SPINK13_t$log_SPINK13.201, .75)
IQR <- IQR(cancer_SPINK13_t$log_SPINK13.201)

# subset data where points value is outside 2*1.5*IQR of Q1 and Q3
outliers <- subset(cancer_SPINK13_t, cancer_SPINK13_t$log_SPINK13.201<(Q1 - 2*1.5*IQR) | cancer_SPINK13_t$log_SPINK13.201>(Q3 + 2*1.5*IQR))
# 71 extreme outliers

# Generate pvalues
stat.test <- long_SPINK13_log %>%
  group_by(variable) %>%
  wilcox_test(value ~ sample_type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

# Plot
plot <- ggplot(long_SPINK13_log, aes(x = variable, y = value, col=sample_type)) +
  geom_boxplot() +
  geom_point(alpha = 0.2) +
  xlab("SPINK13 transcript isoform") +
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
cancer_sample_names = row.names(cancer_SPINK13_t)
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

cancer_SPINK13 <- SPINK13_df_anti %>%
  filter(sample_type == "Cancer")

# Merge survival data with gene counts 
SPINK13_survival_df <- merge(cancer_SPINK13, filtered_tcga_survival, by = "row.names", 
                           all.x = TRUE) # all.x will give NA for missing values


# Convert days to months
# 1 day = 0.0328767 months
days_to_months_converter <- function(no_days){
  no_months <- round(no_days*0.0328767, 1)
  return(no_months)
}

SPINK13_survival_df["OS_months"] <- lapply(SPINK13_survival_df["OS.time"], days_to_months_converter)
SPINK13_survival_df["DSS_months"] <- lapply(SPINK13_survival_df["DSS.time"], days_to_months_converter)
SPINK13_survival_df["PFI_months"] <- lapply(SPINK13_survival_df["PFI.time"], days_to_months_converter)

SPINK13

# Investigating median gene expression category on survival
SPINK13_survival_df$SPINK13.201_median <- NA
SPINK13_survival_df$SPINK13.201_median[SPINK13_survival_df$SPINK13.201 <= quantile(SPINK13_survival_df$SPINK13.201, 0.5, na.rm=TRUE)] <- "LOW"
SPINK13_survival_df$SPINK13.201_median[SPINK13_survival_df$SPINK13.201 > quantile(SPINK13_survival_df$SPINK13.201, 0.5, na.rm=TRUE)] <- "HIGH"

SPINK13_survival_df$SPINK13.202_median <- NA
SPINK13_survival_df$SPINK13.202_median[SPINK13_survival_df$SPINK13.202 <= quantile(SPINK13_survival_df$SPINK13.202, 0.5, na.rm=TRUE)] <- "LOW"
SPINK13_survival_df$SPINK13.202_median[SPINK13_survival_df$SPINK13.202 > quantile(SPINK13_survival_df$SPINK13.202, 0.5, na.rm=TRUE)] <- "HIGH"

SPINK13_survival_df$SPINK13.203_median <- NA
SPINK13_survival_df$SPINK13.203_median[SPINK13_survival_df$SPINK13.203 <= quantile(SPINK13_survival_df$SPINK13.203, 0.5, na.rm=TRUE)] <- "LOW"
SPINK13_survival_df$SPINK13.203_median[SPINK13_survival_df$SPINK13.203 > quantile(SPINK13_survival_df$SPINK13.203, 0.5, na.rm=TRUE)] <- "HIGH"

covariate_names <- c(SPINK13.202_median="SPINK13.202 median",
                     SPINK13.203_median="SPINK13.203 median",
                     SPINK13.201_median="SPINK13.201 median"
)


# OS
SPINK13_median_OS <- SPINK13_survival_df %>%
  select(Row.names, SPINK13.202_median, SPINK13.203_median, SPINK13.201_median, sex, OS_months, OS)
# Remove NAs
SPINK13_median_OS <- subset(SPINK13_median_OS, !is.na(OS) & !is.na(OS_months))

map(vars(SPINK13.202_median, SPINK13.203_median, SPINK13.201_median,sex), function(by)
{
  analyse_multivariate(SPINK13_median_OS,
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
SPINK13_median_DSS <- SPINK13_survival_df %>%
  select(Row.names, SPINK13.202_median, SPINK13.203_median, SPINK13.201_median, sex, DSS_months, DSS)

# Remove NAs
SPINK13_median_DSS <- subset(SPINK13_median_DSS, !is.na(DSS) & !is.na(DSS_months))

map(vars(SPINK13.202_median, SPINK13.203_median, SPINK13.201_median), function(by)
{
  analyse_multivariate(SPINK13_median_DSS,
                       vars(DSS_months, DSS),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = covariate_names)
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="DSS"),
              orderer = ~order(p),
              labels_displayed = c("factor"),
              ggtheme = ggplot2::theme_bw(base_size = 14))

