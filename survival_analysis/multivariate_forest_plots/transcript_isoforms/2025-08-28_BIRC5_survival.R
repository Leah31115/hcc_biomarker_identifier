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

# Reshape
library(reshape2)

# Clear environment
rm(list=ls())

setwd("C:/Users/leahe/Documents/LIFE4137_independent_project/xena")
splice_df <- read.csv(file="splicing/survival/transcripts_survival.csv", row.names=1, check.names = FALSE, header=T)
meta <- read.csv(file="shared_meta.csv", row.names=1, check.names = FALSE, header=T)
tcga_survival <- read.table(file="all_survival_LIHC_data.txt", sep='\t', header=T)


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


# BIRC5
BIRC5_df <- splice_df %>%
  filter(Ensembl_ID %in% c(BIRC5))

# Rename row.names column and make the rows the gene names
rownames(BIRC5_df) <- BIRC5_df$Ensembl_ID
BIRC5_df <- BIRC5_df[,!(names(BIRC5_df) %in% "Ensembl_ID")]

# Make columns (samples) the row name
BIRC5_df <- as.data.frame(t(BIRC5_df))
# Assign sample category
BIRC5_df$sample_type <- NA
name_start <- substr(rownames(BIRC5_df), 1, 4)
BIRC5_df$sample_type <- ifelse(name_start == "GTEX", "Healthy", "Cancer")

# Add meta
BIRC5_df <- merge(BIRC5_df, meta,
                  by = 'row.names', all = TRUE)
names(BIRC5_df)[names(BIRC5_df) == "Row.names"] <- "samples"
rownames(BIRC5_df) <- BIRC5_df$samples
BIRC5_df <- BIRC5_df[,!(names(BIRC5_df) %in% "samples")]

# Change transcript IDs to transcript names
colnames(BIRC5_df)[which(names(BIRC5_df) == "ENST00000301633.8")] <- "BIRC5.201"
colnames(BIRC5_df)[which(names(BIRC5_df) == "ENST00000350051.7")] <- "BIRC5.202"
colnames(BIRC5_df)[which(names(BIRC5_df) == "ENST00000374948.6")] <- "BIRC5.203"
colnames(BIRC5_df)[which(names(BIRC5_df) == "ENST00000586192.5")] <- "BIRC5.204"
colnames(BIRC5_df)[which(names(BIRC5_df) == "ENST00000587746.5")] <- "BIRC5.205"
colnames(BIRC5_df)[which(names(BIRC5_df) == "ENST00000589892.1")] <- "BIRC5.206"
colnames(BIRC5_df)[which(names(BIRC5_df) == "ENST00000590449.1")] <- "BIRC5.207"
colnames(BIRC5_df)[which(names(BIRC5_df) == "ENST00000590925.6")] <- "BIRC5.208"
colnames(BIRC5_df)[which(names(BIRC5_df) == "ENST00000591800.1")] <- "BIRC5.209"
colnames(BIRC5_df)[which(names(BIRC5_df) == "ENST00000592115.5")] <- "BIRC5.210"
colnames(BIRC5_df)[which(names(BIRC5_df) == "ENST00000592734.5")] <- "BIRC5.211"

# Add log prefix to transcript names
BIRC5_df <- BIRC5_df %>% 
  rename_with(function(x) paste0("log_", x), -c(sex, sample_type))
colnames(BIRC5_df)


# Antilog gene count function
antilog_counts <- function(tpm){
  antilog_count <- (2^tpm) - 0.001
  return(antilog_count)
}

# Antilog
BIRC5_df["BIRC5.201"] <- lapply(BIRC5_df["log_BIRC5.201"], antilog_counts)
BIRC5_df["BIRC5.202"] <- lapply(BIRC5_df["log_BIRC5.202"], antilog_counts)
BIRC5_df["BIRC5.208"] <- lapply(BIRC5_df["log_BIRC5.208"], antilog_counts)
BIRC5_df["BIRC5.203"] <- lapply(BIRC5_df["log_BIRC5.203"], antilog_counts)
BIRC5_df["BIRC5.207"] <- lapply(BIRC5_df["log_BIRC5.207"], antilog_counts)
BIRC5_df["BIRC5.211"] <- lapply(BIRC5_df["log_BIRC5.211"], antilog_counts)
BIRC5_df["BIRC5.210"] <- lapply(BIRC5_df["log_BIRC5.210"], antilog_counts)
BIRC5_df["BIRC5.204"] <- lapply(BIRC5_df["log_BIRC5.204"], antilog_counts)
BIRC5_df["BIRC5.205"] <- lapply(BIRC5_df["log_BIRC5.205"], antilog_counts)
BIRC5_df["BIRC5.209"] <- lapply(BIRC5_df["log_BIRC5.209"], antilog_counts)
BIRC5_df["BIRC5.206"] <- lapply(BIRC5_df["log_BIRC5.206"], antilog_counts)


# subset for antilogged counts
BIRC5_df_anti <- BIRC5_df %>% 
  select(BIRC5.201,
         BIRC5.202,
         BIRC5.208,
         BIRC5.203,
         BIRC5.207,
         BIRC5.211,
         BIRC5.210,
         BIRC5.204,
         BIRC5.205,
         BIRC5.209,
         BIRC5.206,
         sample_type,
         sex)

# Melt df
long_BIRC5_anti <- melt(BIRC5_df_anti, id_vars="sample_type")

# Problem with negative antilog counts and logging the y-axis
ggplot(long_BIRC5_anti, aes(x = variable, y = value, col=sample_type)) +
  geom_point(alpha=0.1, position="jitter") +
  scale_y_continuous(trans = "log2") +
  geom_boxplot() +
  xlab("BIRC5 transcript variants") +
  ylab("Anti log TPM gene counts") +
  theme_bw(base_size = 16)

# Subset for logged counts
log_BIRC5_df <- BIRC5_df %>% 
  select(log_BIRC5.201,
         log_BIRC5.202,
         log_BIRC5.208,
         log_BIRC5.203,
         log_BIRC5.207,
         log_BIRC5.211,
         log_BIRC5.210,
         log_BIRC5.204,
         log_BIRC5.205,
         log_BIRC5.209,
         log_BIRC5.206,
         sample_type,
         sex)

# Remove log prefix for graph formatting
names(log_BIRC5_df) <- sub('log_', '', names(log_BIRC5_df))

long_BIRC5_log <- melt(log_BIRC5_df, id_vars="sample_type")

ggplot(long_BIRC5_log, aes(x = variable, y = value, col=sample_type)) +
  geom_boxplot() +
  geom_point(alpha=0.09, position="jitter") +
  xlab("BIRC5 transcript isoform") +
  ylab("Log TPM gene counts") +
  theme_bw(base_size = 16)

# Check count distributions
# BIRC5 protein coding gene
BIRC5_distribution <- BIRC5_df %>% 
  select(log_BIRC5.201, sample_type)
# Summary stats
group_by(BIRC5_distribution, sample_type) %>%
  summarise(
    count = n(),
    median = median(log_BIRC5.201, na.rm = TRUE),
    IQR = IQR(log_BIRC5.201, na.rm = TRUE)
  )
# Healthy distribution
healthy_BIRC5 <- BIRC5_distribution %>%
  filter(sample_type == "Healthy")

ggdensity(healthy_BIRC5$log_BIRC5.201, 
          main = "Density plot of BIRC5 gene from healthy samples",
          xlab = "Tpm normalised gene count")

shapiro.test(healthy_BIRC5$log_BIRC5.201) 
#  significant so non-parametric distribution

# Identifying extreme outliers
# Find Q1, Q3, and interquartile range for values
Q1 <- quantile(healthy_BIRC5$log_BIRC5.201, .25)
Q3 <- quantile(healthy_BIRC5$log_BIRC5.201, .75)
IQR <- IQR(healthy_BIRC5$log_BIRC5.201)

# subset data where points value is outside 2*1.5*IQR of Q1 and Q3
outliers <- subset(healthy_BIRC5, healthy_BIRC5$log_BIRC5.201<(Q1 - 2*1.5*IQR) | healthy_BIRC5$log_BIRC5.201>(Q3 + 2*1.5*IQR))
# No extreme outliers

# cancer distribution
cancer_BIRC5_t <- BIRC5_distribution %>%
  filter(sample_type == "Cancer")

ggdensity(cancer_BIRC5_t$log_BIRC5.201, 
          main = "Density plot of BIRC5 gene from cancer samples",
          xlab = "Tpm normalised gene count")

shapiro.test(cancer_BIRC5_t$log_BIRC5.201) 
# significant so non-parametric distribution

# Identifying extreme outliers
# Find Q1, Q3, and interquartile range for values
Q1 <- quantile(cancer_BIRC5_t$log_BIRC5.201, .25)
Q3 <- quantile(cancer_BIRC5_t$log_BIRC5.201, .75)
IQR <- IQR(cancer_BIRC5_t$log_BIRC5.201)

# subset data where points value is outside 2*1.5*IQR of Q1 and Q3
outliers <- subset(cancer_BIRC5_t, cancer_BIRC5_t$log_BIRC5.201<(Q1 - 2*1.5*IQR) | cancer_BIRC5_t$log_BIRC5.201>(Q3 + 2*1.5*IQR))
# 23 extreme outliers

# Generate pvalues
stat.test <- long_BIRC5_log %>%
  group_by(variable) %>%
  wilcox_test(value ~ sample_type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

# Reorder transcripts to match the multivariate forest plot
long_BIRC5_log$variable <-factor(long_BIRC5_log$variable, levels=rev(levels(long_BIRC5_log$variable)))

# Plot
plot <- ggplot(long_BIRC5_log, aes(x = variable, y = value, col=sample_type)) +
  geom_boxplot() +
  geom_point(alpha = 0.2) +
  xlab("BIRC5 transcript isoform") +
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
cancer_sample_names = row.names(cancer_BIRC5_t)
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

cancer_BIRC5 <- BIRC5_df_anti %>%
  filter(sample_type == "Cancer")

# Merge survival data with gene counts 
BIRC5_survival_df <- merge(cancer_BIRC5, filtered_tcga_survival, by = "row.names", 
                           all.x = TRUE) # all.x will give NA for missing values


# Convert days to months
# 1 day = 0.0328767 months
days_to_months_converter <- function(no_days){
  no_months <- round(no_days*0.0328767, 1)
  return(no_months)
}

BIRC5_survival_df["OS_months"] <- lapply(BIRC5_survival_df["OS.time"], days_to_months_converter)
BIRC5_survival_df["DSS_months"] <- lapply(BIRC5_survival_df["DSS.time"], days_to_months_converter)

# Investigating median gene expression category on survival
BIRC5_survival_df$BIRC5.201_median <- NA
BIRC5_survival_df$BIRC5.201_median[BIRC5_survival_df$BIRC5.201 <= quantile(BIRC5_survival_df$BIRC5.201, 0.5, na.rm=TRUE)] <- "LOW"
BIRC5_survival_df$BIRC5.201_median[BIRC5_survival_df$BIRC5.201 > quantile(BIRC5_survival_df$BIRC5.201, 0.5, na.rm=TRUE)] <- "HIGH"

BIRC5_survival_df$BIRC5.202_median <- NA
BIRC5_survival_df$BIRC5.202_median[BIRC5_survival_df$BIRC5.202 <= quantile(BIRC5_survival_df$BIRC5.202, 0.5, na.rm=TRUE)] <- "LOW"
BIRC5_survival_df$BIRC5.202_median[BIRC5_survival_df$BIRC5.202 > quantile(BIRC5_survival_df$BIRC5.202, 0.5, na.rm=TRUE)] <- "HIGH"

BIRC5_survival_df$BIRC5.208_median <- NA
BIRC5_survival_df$BIRC5.208_median[BIRC5_survival_df$BIRC5.208 <= quantile(BIRC5_survival_df$BIRC5.208, 0.5, na.rm=TRUE)] <- "LOW"
BIRC5_survival_df$BIRC5.208_median[BIRC5_survival_df$BIRC5.208 > quantile(BIRC5_survival_df$BIRC5.208, 0.5, na.rm=TRUE)] <- "HIGH"

BIRC5_survival_df$BIRC5.203_median <- NA
BIRC5_survival_df$BIRC5.203_median[BIRC5_survival_df$BIRC5.203 <= quantile(BIRC5_survival_df$BIRC5.203, 0.5, na.rm=TRUE)] <- "LOW"
BIRC5_survival_df$BIRC5.203_median[BIRC5_survival_df$BIRC5.203 > quantile(BIRC5_survival_df$BIRC5.203, 0.5, na.rm=TRUE)] <- "HIGH"

BIRC5_survival_df$BIRC5.207_median <- NA
BIRC5_survival_df$BIRC5.207_median[BIRC5_survival_df$BIRC5.207 <= quantile(BIRC5_survival_df$BIRC5.207, 0.5, na.rm=TRUE)] <- "LOW"
BIRC5_survival_df$BIRC5.207_median[BIRC5_survival_df$BIRC5.207 > quantile(BIRC5_survival_df$BIRC5.207, 0.5, na.rm=TRUE)] <- "HIGH"

BIRC5_survival_df$BIRC5.211_median <- NA
BIRC5_survival_df$BIRC5.211_median[BIRC5_survival_df$BIRC5.211 <= quantile(BIRC5_survival_df$BIRC5.211, 0.5, na.rm=TRUE)] <- "LOW"
BIRC5_survival_df$BIRC5.211_median[BIRC5_survival_df$BIRC5.211 > quantile(BIRC5_survival_df$BIRC5.211, 0.5, na.rm=TRUE)] <- "HIGH"

BIRC5_survival_df$BIRC5.210_median <- NA
BIRC5_survival_df$BIRC5.210_median[BIRC5_survival_df$BIRC5.210 <= quantile(BIRC5_survival_df$BIRC5.210, 0.5, na.rm=TRUE)] <- "LOW"
BIRC5_survival_df$BIRC5.210_median[BIRC5_survival_df$BIRC5.210 > quantile(BIRC5_survival_df$BIRC5.210, 0.5, na.rm=TRUE)] <- "HIGH"

BIRC5_survival_df$BIRC5.204_median <- NA
BIRC5_survival_df$BIRC5.204_median[BIRC5_survival_df$BIRC5.204 <= quantile(BIRC5_survival_df$BIRC5.204, 0.5, na.rm=TRUE)] <- "LOW"
BIRC5_survival_df$BIRC5.204_median[BIRC5_survival_df$BIRC5.204 > quantile(BIRC5_survival_df$BIRC5.204, 0.5, na.rm=TRUE)] <- "HIGH"

BIRC5_survival_df$BIRC5.205_median <- NA
BIRC5_survival_df$BIRC5.205_median[BIRC5_survival_df$BIRC5.205 <= quantile(BIRC5_survival_df$BIRC5.205, 0.5, na.rm=TRUE)] <- "LOW"
BIRC5_survival_df$BIRC5.205_median[BIRC5_survival_df$BIRC5.205 > quantile(BIRC5_survival_df$BIRC5.205, 0.5, na.rm=TRUE)] <- "HIGH"

BIRC5_survival_df$BIRC5.209_median <- NA
BIRC5_survival_df$BIRC5.209_median[BIRC5_survival_df$BIRC5.209 <= quantile(BIRC5_survival_df$BIRC5.209, 0.5, na.rm=TRUE)] <- "LOW"
BIRC5_survival_df$BIRC5.209_median[BIRC5_survival_df$BIRC5.209 > quantile(BIRC5_survival_df$BIRC5.209, 0.5, na.rm=TRUE)] <- "HIGH"

BIRC5_survival_df$BIRC5.206_median <- NA
BIRC5_survival_df$BIRC5.206_median[BIRC5_survival_df$BIRC5.206 <= quantile(BIRC5_survival_df$BIRC5.206, 0.5, na.rm=TRUE)] <- "LOW"
BIRC5_survival_df$BIRC5.206_median[BIRC5_survival_df$BIRC5.206 > quantile(BIRC5_survival_df$BIRC5.206, 0.5, na.rm=TRUE)] <- "HIGH"


covariate_names <- c(BIRC5.202_median="BIRC5.202 median",
                            BIRC5.208_median="BIRC5.208 median",
                            BIRC5.203_median="BIRC5.203 median",
                            BIRC5.201_median="BIRC5.201 median",
                            BIRC5.207_median="BIRC5.207 median",
                            BIRC5.211_median="BIRC5.211 median",
                            BIRC5.210_median="BIRC5.210 median",
                            BIRC5.204_median="BIRC5.204 median",
                            BIRC5.205_median="BIRC5.205 median",
                            BIRC5.209_median="BIRC5.209 median",
                            BIRC5.206_median="BIRC5.206 median"
)

# OS
BIRC5_median_OS <- BIRC5_survival_df %>%
  select(Row.names, BIRC5.202_median, BIRC5.208_median, BIRC5.203_median, BIRC5.201_median, BIRC5.207_median, BIRC5.211_median, BIRC5.210_median, BIRC5.204_median, BIRC5.205_median, BIRC5.209_median, BIRC5.206_median, sex, OS_months, OS)
# Remove NAs
BIRC5_median_OS <- subset(BIRC5_median_OS, !is.na(OS) & !is.na(OS_months))

map(vars(BIRC5.202_median, BIRC5.208_median, BIRC5.203_median, BIRC5.201_median, BIRC5.207_median, BIRC5.211_median, BIRC5.210_median, BIRC5.204_median, BIRC5.205_median, BIRC5.209_median, BIRC5.206_median), function(by)
{
  analyse_multivariate(BIRC5_median_OS,
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
BIRC5_median_DSS <- BIRC5_survival_df %>%
  select(Row.names, BIRC5.202_median, BIRC5.208_median, BIRC5.203_median, BIRC5.201_median, BIRC5.207_median, BIRC5.211_median, BIRC5.210_median, BIRC5.204_median, BIRC5.205_median, BIRC5.209_median, BIRC5.206_median, sex, DSS_months, DSS)

# Remove NAs
BIRC5_median_DSS <- subset(BIRC5_median_DSS, !is.na(DSS) & !is.na(DSS_months))

map(vars(BIRC5.202_median, BIRC5.208_median, BIRC5.203_median, BIRC5.201_median, BIRC5.207_median, BIRC5.211_median, BIRC5.210_median, BIRC5.204_median, BIRC5.205_median, BIRC5.209_median, BIRC5.206_median), function(by)
{
  analyse_multivariate(BIRC5_median_DSS,
                       vars(DSS_months, DSS),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = covariate_names)
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="DSS"),
              orderer = ~order(p),
              labels_displayed = c("factor"),
              ggtheme = ggplot2::theme_bw(base_size = 14))

