
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


CRHBP <- c("ENST00000274368.8",
           "ENST00000506501.1",
           "ENST00000512446.1",
           "ENST00000514258.1",
           "ENST00000503763.1"
)

# CRHBP
CRHBP_df <- splice_df %>%
  filter(Ensembl_ID %in% c(CRHBP))

# Rename row.names column and make the rows the gene names
rownames(CRHBP_df) <- CRHBP_df$Ensembl_ID
CRHBP_df <- CRHBP_df[,!(names(CRHBP_df) %in% "Ensembl_ID")]

# Make columns (samples) the row name
CRHBP_df <- as.data.frame(t(CRHBP_df))
# Assign sample category
CRHBP_df$sample_type <- NA
name_start <- substr(rownames(CRHBP_df), 1, 4)
CRHBP_df$sample_type <- ifelse(name_start == "GTEX", "Healthy", "Cancer")

# Add meta
CRHBP_df <- merge(CRHBP_df, meta,
                  by = 'row.names', all = TRUE)
names(CRHBP_df)[names(CRHBP_df) == "Row.names"] <- "samples"
rownames(CRHBP_df) <- CRHBP_df$samples
CRHBP_df <- CRHBP_df[,!(names(CRHBP_df) %in% "samples")]

# Change transcript IDs to transcript names
colnames(CRHBP_df)[which(names(CRHBP_df) == "ENST00000274368.8")] <- "CRHBP.201"
colnames(CRHBP_df)[which(names(CRHBP_df) == "ENST00000506501.1")] <- "CRHBP.203"
colnames(CRHBP_df)[which(names(CRHBP_df) == "ENST00000512446.1")] <- "CRHBP.204"
colnames(CRHBP_df)[which(names(CRHBP_df) == "ENST00000514258.1")] <- "CRHBP.205"
colnames(CRHBP_df)[which(names(CRHBP_df) == "ENST00000503763.1")] <- "CRHBP.202"

# Add log prefix to transcript names
CRHBP_df <- CRHBP_df %>% 
  rename_with(function(x) paste0("log_", x), -c(sex, sample_type))
colnames(CRHBP_df)


# Antilog gene count function
antilog_counts <- function(tpm){
  antilog_count <- (2^tpm) - 0.001
  return(antilog_count)
}

# Antilog
CRHBP_df["CRHBP.201"] <- lapply(CRHBP_df["log_CRHBP.201"], antilog_counts)
CRHBP_df["CRHBP.203"] <- lapply(CRHBP_df["log_CRHBP.203"], antilog_counts)
CRHBP_df["CRHBP.204"] <- lapply(CRHBP_df["log_CRHBP.204"], antilog_counts)
CRHBP_df["CRHBP.205"] <- lapply(CRHBP_df["log_CRHBP.205"], antilog_counts)
CRHBP_df["CRHBP.202"] <- lapply(CRHBP_df["log_CRHBP.202"], antilog_counts)


# subset for antilogged counts
CRHBP_df_anti <- CRHBP_df %>% 
  select(CRHBP.201,
         CRHBP.203,
         CRHBP.204,
         CRHBP.205,
         CRHBP.202,
         sample_type,
         sex)


# Melt df
long_CRHBP_anti <- melt(CRHBP_df_anti, id_vars="sample_type")

# Problem with negative antilog counts and logging the y-axis
ggplot(long_CRHBP_anti, aes(x = variable, y = value, col=sample_type)) +
  geom_point(alpha=0.1, position="jitter") +
  geom_boxplot() +
  scale_y_continuous(trans = "log2") +
  xlab("CRHBP splice variants") +
  ylab("Anti log TPM gene counts") +
  theme_bw()

# Subset for logged counts
log_CRHBP_df <- CRHBP_df %>% 
  select(log_CRHBP.201,
         log_CRHBP.203,
         log_CRHBP.204,
         log_CRHBP.205,
         log_CRHBP.202,
         sample_type,
         sex)

# Remove log prefix for graph formatting
names(log_CRHBP_df) <- sub('log_', '', names(log_CRHBP_df))

long_CRHBP_log <- melt(log_CRHBP_df, id_vars="sample_type")

# transcripts on y axis due to ensembl long names
ggplot(long_CRHBP_log, aes(x = value, y = variable, col=sample_type)) +
  geom_boxplot() +
  geom_point(alpha=0.09) +
  xlab("Log TPM gene counts") +
  ylab("CRHBP transcript isoform") +
  theme_bw(base_size = 16)


# Check count distributions
# CRHBP protein coding gene
CRHBP_distribution <- CRHBP_df %>% 
  select(log_CRHBP.201, sample_type)
# Summary stats
group_by(CRHBP_distribution, sample_type) %>%
  summarise(
    count = n(),
    median = median(log_CRHBP.201, na.rm = TRUE),
    IQR = IQR(log_CRHBP.201, na.rm = TRUE)
  )
# Healthy distribution
healthy_CRHBP <- CRHBP_distribution %>%
  filter(sample_type == "Healthy")

ggdensity(healthy_CRHBP$log_CRHBP.201, 
          main = "Density plot of CRHBP gene from healthy samples",
          xlab = "Tpm normalised gene count")

shapiro.test(healthy_CRHBP$log_CRHBP.201) 
# not-significant so non-parametric distribution

# Identifying extreme outliers
# Find Q1, Q3, and interquartile range for values
Q1 <- quantile(healthy_CRHBP$log_CRHBP.201, .25)
Q3 <- quantile(healthy_CRHBP$log_CRHBP.201, .75)
IQR <- IQR(healthy_CRHBP$log_CRHBP.201)

# subset data where points value is outside 2*1.5*IQR of Q1 and Q3
outliers <- subset(healthy_CRHBP, healthy_CRHBP$log_CRHBP.201<(Q1 - 2*1.5*IQR) | healthy_CRHBP$log_CRHBP.201>(Q3 + 2*1.5*IQR))
# 0 extreme outliers

# cancer distribution
cancer_CRHBP_t <- CRHBP_distribution %>%
  filter(sample_type == "Cancer")

ggdensity(cancer_CRHBP_t$log_CRHBP.201, 
          main = "Density plot of CRHBP gene from cancer samples",
          xlab = "Tpm normalised gene count")

shapiro.test(cancer_CRHBP_t$log_CRHBP.201) 
# significant so non-parametric distribution

# Identifying extreme outliers
# Find Q1, Q3, and interquartile range for values
Q1 <- quantile(cancer_CRHBP_t$log_CRHBP.201, .25)
Q3 <- quantile(cancer_CRHBP_t$log_CRHBP.201, .75)
IQR <- IQR(cancer_CRHBP_t$log_CRHBP.201)

# subset data where points value is outside 2*1.5*IQR of Q1 and Q3
outliers <- subset(cancer_CRHBP_t, cancer_CRHBP_t$log_CRHBP.201<(Q1 - 2*1.5*IQR) | cancer_CRHBP_t$log_CRHBP.201>(Q3 + 2*1.5*IQR))
# 0 extreme outliers

# Generate pvalues
stat.test <- long_CRHBP_log %>%
  group_by(variable) %>%
  wilcox_test(value ~ sample_type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

# Plot
plot <- ggplot(long_CRHBP_log, aes(x = variable, y = value, col=sample_type)) +
  geom_boxplot() +
  geom_point(alpha = 0.2) +
  xlab("CRHBP transcript isoform") +
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
cancer_sample_names = row.names(cancer_CRHBP_t)
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

cancer_CRHBP <- CRHBP_df_anti %>%
  filter(sample_type == "Cancer")

# Merge survival data with gene counts 
CRHBP_survival_df <- merge(cancer_CRHBP, filtered_tcga_survival, by = "row.names", 
                           all.x = TRUE) # all.x will give NA for missing values


# Convert days to months
# 1 day = 0.0328767 months
days_to_months_converter <- function(no_days){
  no_months <- round(no_days*0.0328767, 1)
  return(no_months)
}

CRHBP_survival_df["OS_months"] <- lapply(CRHBP_survival_df["OS.time"], days_to_months_converter)
CRHBP_survival_df["DSS_months"] <- lapply(CRHBP_survival_df["DSS.time"], days_to_months_converter)

CRHBP

# Investigating median gene expression category on survival
CRHBP_survival_df$CRHBP.201_median <- NA
CRHBP_survival_df$CRHBP.201_median[CRHBP_survival_df$CRHBP.201 <= quantile(CRHBP_survival_df$CRHBP.201, 0.5, na.rm=TRUE)] <- "LOW"
CRHBP_survival_df$CRHBP.201_median[CRHBP_survival_df$CRHBP.201 > quantile(CRHBP_survival_df$CRHBP.201, 0.5, na.rm=TRUE)] <- "HIGH"

CRHBP_survival_df$CRHBP.203_median <- NA
CRHBP_survival_df$CRHBP.203_median[CRHBP_survival_df$CRHBP.203 <= quantile(CRHBP_survival_df$CRHBP.203, 0.5, na.rm=TRUE)] <- "LOW"
CRHBP_survival_df$CRHBP.203_median[CRHBP_survival_df$CRHBP.203 > quantile(CRHBP_survival_df$CRHBP.203, 0.5, na.rm=TRUE)] <- "HIGH"

CRHBP_survival_df$CRHBP.204_median <- NA
CRHBP_survival_df$CRHBP.204_median[CRHBP_survival_df$CRHBP.204 <= quantile(CRHBP_survival_df$CRHBP.204, 0.5, na.rm=TRUE)] <- "LOW"
CRHBP_survival_df$CRHBP.204_median[CRHBP_survival_df$CRHBP.204 > quantile(CRHBP_survival_df$CRHBP.204, 0.5, na.rm=TRUE)] <- "HIGH"

CRHBP_survival_df$CRHBP.205_median <- NA
CRHBP_survival_df$CRHBP.205_median[CRHBP_survival_df$CRHBP.205 <= quantile(CRHBP_survival_df$CRHBP.205, 0.5, na.rm=TRUE)] <- "LOW"
CRHBP_survival_df$CRHBP.205_median[CRHBP_survival_df$CRHBP.205 > quantile(CRHBP_survival_df$CRHBP.205, 0.5, na.rm=TRUE)] <- "HIGH"

CRHBP_survival_df$CRHBP.202_median <- NA
CRHBP_survival_df$CRHBP.202_median[CRHBP_survival_df$CRHBP.202 <= quantile(CRHBP_survival_df$CRHBP.202, 0.5, na.rm=TRUE)] <- "LOW"
CRHBP_survival_df$CRHBP.202_median[CRHBP_survival_df$CRHBP.202 > quantile(CRHBP_survival_df$CRHBP.202, 0.5, na.rm=TRUE)] <- "HIGH"

covariate_names <- c(CRHBP.203_median="CRHBP.203 median",
                     CRHBP.204_median="CRHBP.204 median",
                     CRHBP.205_median="CRHBP.205 median",
                     CRHBP.201_median="CRHBP.201 median",
                     CRHBP.202_median="CRHBP.202 median"
)

# OS
CRHBP_median_OS <- CRHBP_survival_df %>%
  select(Row.names, CRHBP.203_median, CRHBP.204_median, CRHBP.205_median, CRHBP.201_median, CRHBP.202_median, sex, OS_months, OS)
# Remove NAs
CRHBP_median_OS <- subset(CRHBP_median_OS, !is.na(OS) & !is.na(OS_months))

map(vars(CRHBP.203_median, CRHBP.204_median, CRHBP.205_median, CRHBP.201_median, CRHBP.202_median), function(by)
{
  analyse_multivariate(CRHBP_median_OS,
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
CRHBP_median_DSS <- CRHBP_survival_df %>%
  select(Row.names, CRHBP.203_median, CRHBP.204_median, CRHBP.205_median, CRHBP.201_median, CRHBP.202_median, sex, DSS_months, DSS)

# Remove NAs
CRHBP_median_DSS <- subset(CRHBP_median_DSS, !is.na(DSS) & !is.na(DSS_months))

map(vars(CRHBP.203_median, CRHBP.204_median, CRHBP.205_median, CRHBP.201_median, CRHBP.202_median), function(by)
{
  analyse_multivariate(CRHBP_median_DSS,
                       vars(DSS_months, DSS),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = covariate_names)
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="DSS"),
              orderer = ~order(p),
              labels_displayed = c("factor"),
              ggtheme = ggplot2::theme_bw(base_size = 14))

