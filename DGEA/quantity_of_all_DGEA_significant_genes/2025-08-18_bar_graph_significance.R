library(tidyverse)
library(reshape2)

# Clear environment
rm(list=ls())

setwd("C:/Users/leahe/Documents/LIFE4137_independent_project/xena")
dgea <- read.csv(file="2025_08_18_bar_chart_significant_DGEA.csv", header=T, row.names = 1)

# Transpose df
dgea <- as.data.frame(t(dgea))
dgea <- tibble::rownames_to_column(dgea, "condition")

dgea <- melt(dgea, id_vars="condition")
names(dgea)[names(dgea) == 'variable'] <- 'Gene'

# Filter for healthy vs cancer comparisons (including paired sample comparisons)
hcancer <- dgea %>%
  filter(condition %in% c("HCC_vs_healthy", "StageI_vs_healthy", "StageII_vs_healthy", "StageIII_vs_healthy", "StageIV_vs_healthy", "HCC_vs_STN"))

# Reorder following a precise order
hcancer <- hcancer %>%
  mutate(condition = fct_relevel(condition, 
                            "HCC_vs_healthy", "StageI_vs_healthy", "StageII_vs_healthy", 
                            "StageIII_vs_healthy", "StageIV_vs_healthy", "HCC_vs_STN", 
                            )
         )

ggplot(hcancer, aes(x=condition, y=value, fill=Gene)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw(base_size = 14) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  ylab("Gene quantity") +
  xlab("DGEA") +
  scale_fill_manual("Genes", values = c("Significant" = "darkolivegreen3", "Upregulated" = "firebrick2", "Downregulated" = "cornflowerblue"))


# Filter for cancer stage vs stage DGEA comparisons
cancer_stages <- dgea %>%
  filter(condition %in% c("Advanced_HCC_vs_StageI", "StageII_vs_StageI", "StageIII_vs_StageI", "StageIV_vs_StageI"
                          )
         )

# Reorder following a precise order
cancer_stages <- cancer_stages %>%
  mutate(condition = fct_relevel(condition, 
                                 "Advanced_HCC_vs_StageI", "StageII_vs_StageI", "StageIII_vs_StageI", "StageIV_vs_StageI"
                                 )
         )
ggplot(cancer_stages, aes(x=condition, y=value, fill=Gene)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw(base_size = 14) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  ylab("Gene quantity") +
  xlab("DGEA") +
  scale_fill_manual("Genes", values = c("Significant" = "darkolivegreen3", "Upregulated" = "firebrick2", "Downregulated" = "cornflowerblue"))
