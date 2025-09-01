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

CDC20 <- c("ENST00000310955.10",
           "ENST00000372462.1",
           "ENST00000478882.1",
           "ENST00000482046.1"
)

# CDC20
CDC20_df <- splice_df %>%
  filter(Ensembl_ID %in% c(CDC20))

# Rename row.names column and make the rows the gene names
rownames(CDC20_df) <- CDC20_df$Ensembl_ID
CDC20_df <- CDC20_df[,!(names(CDC20_df) %in% "Ensembl_ID")]

# Make columns (samples) the row name
CDC20_df <- as.data.frame(t(CDC20_df))
# Assign sample category
CDC20_df$sample_type <- NA
name_start <- substr(rownames(CDC20_df), 1, 4)
CDC20_df$sample_type <- ifelse(name_start == "GTEX", "Healthy", "Cancer")

# Add meta
CDC20_df <- merge(CDC20_df, meta,
                  by = 'row.names', all = TRUE)
names(CDC20_df)[names(CDC20_df) == "Row.names"] <- "samples"
rownames(CDC20_df) <- CDC20_df$samples
CDC20_df <- CDC20_df[,!(names(CDC20_df) %in% "samples")]

# Change transcript IDs to transcript names
colnames(CDC20_df)[which(names(CDC20_df) == "ENST00000310955.10")] <- "CDC20.201"
colnames(CDC20_df)[which(names(CDC20_df) == "ENST00000372462.1")] <- "CDC20.202"
colnames(CDC20_df)[which(names(CDC20_df) == "ENST00000478882.1")] <- "CDC20.203"
colnames(CDC20_df)[which(names(CDC20_df) == "ENST00000482046.1")] <- "CDC20.204"

# Add log prefix to transcript names
CDC20_df <- CDC20_df %>% 
  rename_with(function(x) paste0("log_", x), -c(sex, sample_type))
colnames(CDC20_df)


# Antilog gene count function
antilog_counts <- function(tpm){
  antilog_count <- (2^tpm) - 0.001
  return(antilog_count)
}

# Antilog
CDC20_df["CDC20.201"] <- lapply(CDC20_df["log_CDC20.201"], antilog_counts)
CDC20_df["CDC20.202"] <- lapply(CDC20_df["log_CDC20.202"], antilog_counts)
CDC20_df["CDC20.203"] <- lapply(CDC20_df["log_CDC20.203"], antilog_counts)
CDC20_df["CDC20.204"] <- lapply(CDC20_df["log_CDC20.204"], antilog_counts)

# subset for antilogged counts
CDC20_df_anti <- CDC20_df %>% 
  select(CDC20.201,
         CDC20.202,
         CDC20.203,
         CDC20.204,
         sample_type,
         sex)

# Melt df
long_CDC20_anti <- melt(CDC20_df_anti, id_vars="sample_type")

# Problem with negative antilog counts and logging the y-axis
ggplot(long_CDC20_anti, aes(x = variable, y = value, col=sample_type)) +
  geom_point(alpha=0.1, position="jitter") +
  geom_boxplot() +
  scale_y_continuous(trans = "log2") +
  xlab("CDC20 splice variants") +
  ylab("Anti log TPM gene counts") +
  theme_bw()

# Subset for logged counts
log_CDC20_df <- CDC20_df %>% 
  select(log_CDC20.201,
         log_CDC20.202,
         log_CDC20.203,
         log_CDC20.204,
         sample_type,
         sex)

# Remove log prefix for graph formatting
names(log_CDC20_df) <- sub('log_', '', names(log_CDC20_df))

long_CDC20_log <- melt(log_CDC20_df, id_vars="sample_type")

# transcripts on y axis due to ensembl long names
ggplot(long_CDC20_log, aes(x = value, y = variable, col=sample_type)) +
  geom_boxplot() +
  geom_point(alpha=0.09) +
  xlab("Log TPM gene counts") +
  ylab("CDC20 transcript variants") +
  theme_bw(base_size = 16)


# Check count distributions
# CDC20 protein coding gene
CDC20_distribution <- CDC20_df %>% 
  select(log_CDC20.201, sample_type)
# Summary stats
group_by(CDC20_distribution, sample_type) %>%
  summarise(
    count = n(),
    median = median(log_CDC20.201, na.rm = TRUE),
    IQR = IQR(log_CDC20.201, na.rm = TRUE)
  )
# Healthy distribution
healthy_CDC20 <- CDC20_distribution %>%
  filter(sample_type == "Healthy")

ggdensity(healthy_CDC20$log_CDC20.201, 
          main = "Density plot of CDC20 gene from healthy samples",
          xlab = "Tpm normalised gene count")

shapiro.test(healthy_CDC20$log_CDC20.201) 
# significant so non-parametric distribution

# Identifying extreme outliers
# Find Q1, Q3, and interquartile range for values
Q1 <- quantile(healthy_CDC20$log_CDC20.201, .25)
Q3 <- quantile(healthy_CDC20$log_CDC20.201, .75)
IQR <- IQR(healthy_CDC20$log_CDC20.201)

# subset data where points value is outside 2*1.5*IQR of Q1 and Q3
outliers <- subset(healthy_CDC20, healthy_CDC20$log_CDC20.201<(Q1 - 2*1.5*IQR) | healthy_CDC20$log_CDC20.201>(Q3 + 2*1.5*IQR))
# 0 extreme outliers

# cancer distribution
cancer_CDC20_t <- CDC20_distribution %>%
  filter(sample_type == "Cancer")

ggdensity(cancer_CDC20_t$log_CDC20.201, 
          main = "Density plot of CDC20 gene from cancer samples",
          xlab = "Tpm normalised gene count")

shapiro.test(cancer_CDC20_t$log_CDC20.201) 
# not significant so parametric distribution

# Identifying extreme outliers
# Find Q1, Q3, and interquartile range for values
Q1 <- quantile(cancer_CDC20_t$log_CDC20.201, .25)
Q3 <- quantile(cancer_CDC20_t$log_CDC20.201, .75)
IQR <- IQR(cancer_CDC20_t$log_CDC20.201)

# subset data where points value is outside 2*1.5*IQR of Q1 and Q3
outliers <- subset(cancer_CDC20_t, cancer_CDC20_t$log_CDC20.201<(Q1 - 2*1.5*IQR) | cancer_CDC20_t$log_CDC20.201>(Q3 + 2*1.5*IQR))
# 0 extreme outliers

# Generate pvalues
stat.test <- long_CDC20_log %>%
  group_by(variable) %>%
  wilcox_test(value ~ sample_type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

# Plot
plot <- ggplot(long_CDC20_log, aes(x = variable, y = value, col=sample_type)) +
  geom_boxplot() +
  geom_point(alpha = 0.2) +
  xlab("CDC20 transcript isoform") +
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
cancer_sample_names = row.names(cancer_CDC20_t)
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

cancer_CDC20 <- CDC20_df_anti %>%
  filter(sample_type == "Cancer")

# Merge survival data with gene counts 
CDC20_survival_df <- merge(cancer_CDC20, filtered_tcga_survival, by = "row.names", 
                           all.x = TRUE) # all.x will give NA for missing values


# Convert days to months
# 1 day = 0.0328767 months
days_to_months_converter <- function(no_days){
  no_months <- round(no_days*0.0328767, 1)
  return(no_months)
}

CDC20_survival_df["OS_months"] <- lapply(CDC20_survival_df["OS.time"], days_to_months_converter)
CDC20_survival_df["DSS_months"] <- lapply(CDC20_survival_df["DSS.time"], days_to_months_converter)

# Investigating median gene expression category on survival
CDC20_survival_df$CDC20.201_median <- NA
CDC20_survival_df$CDC20.201_median[CDC20_survival_df$CDC20.201 <= quantile(CDC20_survival_df$CDC20.201, 0.5, na.rm=TRUE)] <- "LOW"
CDC20_survival_df$CDC20.201_median[CDC20_survival_df$CDC20.201 > quantile(CDC20_survival_df$CDC20.201, 0.5, na.rm=TRUE)] <- "HIGH"

CDC20_survival_df$CDC20.202_median <- NA
CDC20_survival_df$CDC20.202_median[CDC20_survival_df$CDC20.202 <= quantile(CDC20_survival_df$CDC20.202, 0.5, na.rm=TRUE)] <- "LOW"
CDC20_survival_df$CDC20.202_median[CDC20_survival_df$CDC20.202 > quantile(CDC20_survival_df$CDC20.202, 0.5, na.rm=TRUE)] <- "HIGH"

CDC20_survival_df$CDC20.203_median <- NA
CDC20_survival_df$CDC20.203_median[CDC20_survival_df$CDC20.203 <= quantile(CDC20_survival_df$CDC20.203, 0.5, na.rm=TRUE)] <- "LOW"
CDC20_survival_df$CDC20.203_median[CDC20_survival_df$CDC20.203 > quantile(CDC20_survival_df$CDC20.203, 0.5, na.rm=TRUE)] <- "HIGH"

CDC20_survival_df$CDC20.204_median <- NA
CDC20_survival_df$CDC20.204_median[CDC20_survival_df$CDC20.204 <= quantile(CDC20_survival_df$CDC20.204, 0.5, na.rm=TRUE)] <- "LOW"
CDC20_survival_df$CDC20.204_median[CDC20_survival_df$CDC20.204 > quantile(CDC20_survival_df$CDC20.204, 0.5, na.rm=TRUE)] <- "HIGH"

covariate_names <- c(CDC20.202_median="CDC20.202 median",
                     CDC20.203_median="CDC20.203 median",
                     CDC20.204_median="CDC20.204 median",
                     CDC20.201_median="CDC20.201 median"
)


# trying with median 
# OS
CDC20_median_OS <- CDC20_survival_df %>%
  select(Row.names, CDC20.202_median, CDC20.203_median, CDC20.204_median, CDC20.201_median, sex, OS_months, OS)
# Remove NAs
CDC20_median_OS <- subset(CDC20_median_OS, !is.na(OS) & !is.na(OS_months))

map(vars(CDC20.202_median, CDC20.203_median, CDC20.204_median, CDC20.201_median,sex), function(by)
{
  analyse_multivariate(CDC20_median_OS,
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
CDC20_median_DSS <- CDC20_survival_df %>%
  select(Row.names, CDC20.202_median, CDC20.203_median, CDC20.204_median, CDC20.201_median, sex, DSS_months, DSS)

# Remove NAs
CDC20_median_DSS <- subset(CDC20_median_DSS, !is.na(DSS) & !is.na(DSS_months))

map(vars(CDC20.202_median, CDC20.203_median, CDC20.204_median, CDC20.201_median), function(by)
{
  analyse_multivariate(CDC20_median_DSS,
                       vars(DSS_months, DSS),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = covariate_names)
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="DSS"),
              orderer = ~order(p),
              labels_displayed = c("factor"),
              ggtheme = ggplot2::theme_bw(base_size = 14))

