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

WDR83OS <- c("ENST00000596731.5",
             "ENST00000222190.9",
             "ENST00000598732.1"
)


# WDR83OS
WDR83OS_df <- splice_df %>%
  filter(Ensembl_ID %in% c(WDR83OS))

# Rename row.names column and make the rows the gene names
rownames(WDR83OS_df) <- WDR83OS_df$Ensembl_ID
WDR83OS_df <- WDR83OS_df[,!(names(WDR83OS_df) %in% "Ensembl_ID")]

# Make columns (samples) the row name
WDR83OS_df <- as.data.frame(t(WDR83OS_df))
# Assign sample category
WDR83OS_df$sample_type <- NA
name_start <- substr(rownames(WDR83OS_df), 1, 4)
WDR83OS_df$sample_type <- ifelse(name_start == "GTEX", "Healthy", "Cancer")

# Add meta
WDR83OS_df <- merge(WDR83OS_df, meta,
                  by = 'row.names', all = TRUE)
names(WDR83OS_df)[names(WDR83OS_df) == "Row.names"] <- "samples"
rownames(WDR83OS_df) <- WDR83OS_df$samples
WDR83OS_df <- WDR83OS_df[,!(names(WDR83OS_df) %in% "samples")]

# Change transcript IDs to transcript names
colnames(WDR83OS_df)[which(names(WDR83OS_df) == "ENST00000596731.5")] <- "WDR83OS.202"
colnames(WDR83OS_df)[which(names(WDR83OS_df) == "ENST00000222190.9")] <- "WDR83OS.201"
colnames(WDR83OS_df)[which(names(WDR83OS_df) == "ENST00000598732.1")] <- "WDR83OS.203"

# Add log prefix to transcript names
WDR83OS_df <- WDR83OS_df %>% 
  rename_with(function(x) paste0("log_", x), -c(sex, sample_type))
colnames(WDR83OS_df)


# Antilog gene count function
antilog_counts <- function(tpm){
  antilog_count <- (2^tpm) - 0.001
  return(antilog_count)
}

# Antilog
WDR83OS_df["WDR83OS.202"] <- lapply(WDR83OS_df["log_WDR83OS.202"], antilog_counts)
WDR83OS_df["WDR83OS.201"] <- lapply(WDR83OS_df["log_WDR83OS.201"], antilog_counts)
WDR83OS_df["WDR83OS.203"] <- lapply(WDR83OS_df["log_WDR83OS.203"], antilog_counts)


# subset for antilogged counts
WDR83OS_df_anti <- WDR83OS_df %>% 
  select(WDR83OS.202,
         WDR83OS.201,
         WDR83OS.203,
         sample_type,
         sex)

# Melt df
long_WDR83OS_anti <- melt(WDR83OS_df_anti, id_vars="sample_type")

# Problem with negative antilog counts and logging the y-axis
ggplot(long_WDR83OS_anti, aes(x = variable, y = value, col=sample_type)) +
  geom_point(alpha=0.1, position="jitter") +
  scale_y_continuous(trans = "log2") +
  geom_boxplot() +
  xlab("WDR83OS splice variants") +
  ylab("Anti log TPM gene counts") +
  theme_bw()

# Subset for logged counts
log_WDR83OS_df <- WDR83OS_df %>% 
  select(log_WDR83OS.202,
         log_WDR83OS.201,
         log_WDR83OS.203,
         sample_type,
         sex)

# Remove log prefix for graph formatting
names(log_WDR83OS_df) <- sub('log_', '', names(log_WDR83OS_df))

long_WDR83OS_log <- melt(log_WDR83OS_df, id_vars="sample_type")


ggplot(long_WDR83OS_log, aes(x = variable, y = value, col=sample_type)) +
  geom_boxplot() +
  geom_point(alpha=0.09, position="jitter") +
  xlab("WDR83OS splice variants") +
  ylab("Log TPM gene counts") +
  theme_bw(base_size = 14)

# Check count distributions
# WDR83OS protein coding gene
WDR83OS_distribution <- WDR83OS_df %>% 
  select(log_WDR83OS.202, sample_type)
# Summary stats
group_by(WDR83OS_distribution, sample_type) %>%
  summarise(
    count = n(),
    median = median(log_WDR83OS.202, na.rm = TRUE),
    IQR = IQR(log_WDR83OS.202, na.rm = TRUE)
  )
# Healthy distribution
healthy_WDR83OS <- WDR83OS_distribution %>%
  filter(sample_type == "Healthy")

ggdensity(healthy_WDR83OS$log_WDR83OS.202, 
          main = "Density plot of WDR83OS gene from healthy samples",
          xlab = "Tpm normalised gene count")

shapiro.test(healthy_WDR83OS$log_WDR83OS.202) 
#  Not significant so parametric distribution

# Identifying extreme outliers
# Find Q1, Q3, and interquartile range for values
Q1 <- quantile(healthy_WDR83OS$log_WDR83OS.202, .25)
Q3 <- quantile(healthy_WDR83OS$log_WDR83OS.202, .75)
IQR <- IQR(healthy_WDR83OS$log_WDR83OS.202)

# subset data where points value is outside 2*1.5*IQR of Q1 and Q3
outliers <- subset(healthy_WDR83OS, healthy_WDR83OS$log_WDR83OS.202<(Q1 - 2*1.5*IQR) | healthy_WDR83OS$log_WDR83OS.202>(Q3 + 2*1.5*IQR))
# No extreme outliers

# cancer distribution
cancer_WDR83OS_t <- WDR83OS_distribution %>%
  filter(sample_type == "Cancer")

ggdensity(cancer_WDR83OS_t$log_WDR83OS.202, 
          main = "Density plot of WDR83OS gene from cancer samples",
          xlab = "Tpm normalised gene count")

shapiro.test(cancer_WDR83OS_t$log_WDR83OS.202) 
# significant so non-parametric distribution

# Identifying extreme outliers
# Find Q1, Q3, and interquartile range for values
Q1 <- quantile(cancer_WDR83OS_t$log_WDR83OS.202, .25)
Q3 <- quantile(cancer_WDR83OS_t$log_WDR83OS.202, .75)
IQR <- IQR(cancer_WDR83OS_t$log_WDR83OS.202)

# subset data where points value is outside 2*1.5*IQR of Q1 and Q3
outliers <- subset(cancer_WDR83OS_t, cancer_WDR83OS_t$log_WDR83OS.202<(Q1 - 2*1.5*IQR) | cancer_WDR83OS_t$log_WDR83OS.202>(Q3 + 2*1.5*IQR))
# 0 extreme outliers

# Generate pvalues
stat.test <- long_WDR83OS_log %>%
  group_by(variable) %>%
  wilcox_test(value ~ sample_type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test


# Downstream boxplot lines
# Plot
plot <- ggplot(long_WDR83OS_log, aes(x = variable, y = value, col=sample_type)) +
  geom_boxplot() +
  geom_point(alpha = 0.2) +
  xlab("WDR83OS transcript isoform") +
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
  symnum.args = list(cutpoints = c(0, 0.05, 1), symbols = c("*", "ns"))
)


# Transcript survival analysis
# Filter survival data for used samples
cancer_sample_names = row.names(cancer_WDR83OS_t)
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

cancer_WDR83OS <- WDR83OS_df_anti %>%
  filter(sample_type == "Cancer")

# Merge survival data with gene counts 
WDR83OS_survival_df <- merge(cancer_WDR83OS, filtered_tcga_survival, by = "row.names", 
                           all.x = TRUE) # all.x will give NA for missing values


# Convert days to months
# 1 day = 0.0328767 months
days_to_months_converter <- function(no_days){
  no_months <- round(no_days*0.0328767, 1)
  return(no_months)
}

WDR83OS_survival_df["OS_months"] <- lapply(WDR83OS_survival_df["OS.time"], days_to_months_converter)
WDR83OS_survival_df["DSS_months"] <- lapply(WDR83OS_survival_df["DSS.time"], days_to_months_converter)

# Investigating median gene expression category on survival
WDR83OS_survival_df$WDR83OS.202_median <- NA
WDR83OS_survival_df$WDR83OS.202_median[WDR83OS_survival_df$WDR83OS.202 <= quantile(WDR83OS_survival_df$WDR83OS.202, 0.5, na.rm=TRUE)] <- "LOW"
WDR83OS_survival_df$WDR83OS.202_median[WDR83OS_survival_df$WDR83OS.202 > quantile(WDR83OS_survival_df$WDR83OS.202, 0.5, na.rm=TRUE)] <- "HIGH"

WDR83OS_survival_df$WDR83OS.201_median <- NA
WDR83OS_survival_df$WDR83OS.201_median[WDR83OS_survival_df$WDR83OS.201 <= quantile(WDR83OS_survival_df$WDR83OS.201, 0.5, na.rm=TRUE)] <- "LOW"
WDR83OS_survival_df$WDR83OS.201_median[WDR83OS_survival_df$WDR83OS.201 > quantile(WDR83OS_survival_df$WDR83OS.201, 0.5, na.rm=TRUE)] <- "HIGH"

WDR83OS_survival_df$WDR83OS.203_median <- NA
WDR83OS_survival_df$WDR83OS.203_median[WDR83OS_survival_df$WDR83OS.203 <= quantile(WDR83OS_survival_df$WDR83OS.203, 0.5, na.rm=TRUE)] <- "LOW"
WDR83OS_survival_df$WDR83OS.203_median[WDR83OS_survival_df$WDR83OS.203 > quantile(WDR83OS_survival_df$WDR83OS.203, 0.5, na.rm=TRUE)] <- "HIGH"


covariate_names <- c(WDR83OS.201_median="WDR83OS.201 median",
                     WDR83OS.203_median="WDR83OS.203 median",
                     WDR83OS.202_median="WDR83OS.202 median"
)

# OS
WDR83OS_median_OS <- WDR83OS_survival_df %>%
  select(Row.names, WDR83OS.201_median, WDR83OS.203_median, WDR83OS.202_median, sex, OS_months, OS)
# Remove NAs
WDR83OS_median_OS <- subset(WDR83OS_median_OS, !is.na(OS) & !is.na(OS_months))

map(vars(WDR83OS.201_median, WDR83OS.203_median, WDR83OS.202_median), function(by)
{
  analyse_multivariate(WDR83OS_median_OS,
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
WDR83OS_median_DSS <- WDR83OS_survival_df %>%
  select(Row.names, WDR83OS.201_median, WDR83OS.203_median, WDR83OS.202_median, sex, DSS_months, DSS)

# Remove NAs
WDR83OS_median_DSS <- subset(WDR83OS_median_DSS, !is.na(DSS) & !is.na(DSS_months))

map(vars(WDR83OS.201_median, WDR83OS.203_median, WDR83OS.202_median), function(by)
{
  analyse_multivariate(WDR83OS_median_DSS,
                       vars(DSS_months, DSS),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = covariate_names)
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="DSS"),
              orderer = ~order(p),
              labels_displayed = c("factor"),
              ggtheme = ggplot2::theme_bw(base_size = 14))

