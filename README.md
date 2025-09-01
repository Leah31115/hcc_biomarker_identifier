# Hepatocellular carcinoma (HCC) biomarker indentification
## Contents
[Overview](#overview) <br/>
[Getting started with R](#Getting-started-with-R) <br/>
[Filter samples](#filter-samples) <br/>
### Differential gene expression analysis (DGEA)
[DGEA](#Differential-Gene-Expression-Analysis-(DGEA)) <br/>
### Survival analysis
[Univariate Kaplan-Meier curves](#Univariate-Kaplan-Meier-Curves) <br/>
[Multivariate Forest plots](#Multivariate-forest-plots) <br/>
<br/>
[Acknowledgements](#Acknowledgements) <br/>
[References](#References) <br/>

## Overview
Identifying early hepatocellular carcinoma (HCC) biomarkers for my thesis. LIFE4137 Independent Project Module, University of Nottingham. The code below was used to to investigate early biomarkers in hepatocellular carcinoma (HCC), a common form of liver cancer with high mortality rates (World Health Organisation, 2024). All data required for downstream analysis is open access and accessible via University of California Santa Cruz (UCSC)  [Xena](https://xenabrowser.net/datapages/) (Goldman et al., 2020). 

## Getting started with R
R studio (version 2024.12.0+467) can be downloaded [here](https://dailies.rstudio.com/version/2024.12.0+467/) and R (version 4.5.1) can be downloaded [here](https://posit.co/download/rstudio-desktop/) (R Core Team, 2025).
Please adjust the file paths according to your work directory when using these R scripts for analysis. To install any packages in R use the following command and replace "package_name" accordingly. 
```R
install.packages("package_name")
```
## Filter samples 
RNAseq STAR gene counts for The Cancer Genome Atlas (TCGA) cancer samples from the GDC TCGA Liver Cancer (LIHC) cohort were downloaded from UCSC Xena [here](https://xenabrowser.net/datapages/?dataset=TCGA-LIHC.star_counts.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) (file TCGA-LIHC.star_counts.tsv.gz, version 05-09-2024, 38.9 MB). This data was used for DGEA. TCGA cancer phenotype metadata, from the same cohort, was also obtained from the UCSC Xena platform [here](https://xenabrowser.net/datapages/?dataset=TCGA-LIHC.clinical.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) (file TCGA-LIHC.clinical.tsv.gz, version 09-07-2024, 92.5 KB).

Healthy sample RNAseq TOIL RSEM expected counts were obtained from the [GTEx cohort](https://xenabrowser.net/datapages/?dataset=gtex_gene_expected_count&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) in UCSC Xena (file gtex_gene_expected_count.gz, version 2016-05-19, 476 MB). This data was used for DGEA. Healthy metadata was also obtained from UCSC Xena [here](https://xenabrowser.net/datapages/?dataset=GTEX_phenotype&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) (file GTEX_phenotype.gz, version 2016-05-18, 80.5 KB).

Samples were filtered using the scripts within the [filter_samples folder](https://github.com/Leah31115/hcc_biomarker_identifier/tree/main/filter_samples). Paired Solid Tissue Normal (STN) healthy samples (ID ending -11) were paired to the cancer sample (ID ending -01) by their TCGA ID.

## Differential Gene Expression Analysis (DGEA)
Some of the code written within the UoN LIFE4136 Group Project 2 was used and expanded upon in this analysis. In which I thank my previous group members Thomas Mclaughlin and Areeba Salman for their past collaborations.
For Deseq2, a DESeq Data Set (DDS) object is made which requires a count data file as well as a condition file. The count data files for all DGEAs are generated from the [filter samples section](#Filter-samples).
The count data file consists of the gene raw counts for each sample. This file contains all gene Ensembl IDs in the first column. All remaining columns contain each sample with their gene count values supplied for each gene, as exampled below:
| EnsemblID | GTEX-1192X-1026-SM-5H12P | GTEX-11DXY-0526-SM-5EGGQ | etc. |
| --- | --- | --- | --- |
| ENSG00000000003 | 11.4424 | 9.9293 | count_value |
| ENSG00000000005 | 0 | 1 | count_value |
| etc. | count_value | count_value | count_value |

The condition file consists of two columns. Column 1 contains all samples, in the same order as what is found within the count data file. The second column contains the condition e.g sample type. An example of the format is displayed below:
| sample | condition |
| --- | --- | 
| GTEX-1192X-1026-SM-5H12P | Healthy |
| GTEX-11DXY-0526-SM-5EGGQ | Healthy |
| etc. | sample_type |

To perform DESeq2 and Gene ontology (GO) analysis, the following packages must be installed: org.Hs.eg.db, AnnotationDbi, and clusterProfiler (Carlson 2025; Pagès et al. 2025; Xu et al. 2024). Run the following, if required, before running the DGEA scripts. Installation and code for GO analysis was used from mousepixels' GO_in_R.Rmd script available [here](https://github.com/mousepixels/sanbomics_scripts/blob/main/GO_in_R.Rmd) (Mousepixels, 2022). 
# Installation for PCA/Heatmap/DESeq2
# Installation
```R
# Installation for DESeq
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("vsn")

# Installation for GO 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
```

DGEA packages used include magrittr, vsn, pheatmap, DESeq2, RColorBrewer, ggrepel, ggplot2, and tidyverse (Bache and Wickham 2022; Huber et al. 2002; Kolde 2025; Love, Huber, and Anders 2014; Neuwirth 2022; Slowikowski 2024; Wickham 2016; Wickham et al. 2019). All DGEAs were run using the scripts in the [DGEA folder](https://github.com/Leah31115/hcc_biomarker_identifier/tree/main/DGEA). Healthy samples were compared against cancer samples using the script [2025-08-11_tcga_vs_healthy_deseq.R](https://github.com/Leah31115/hcc_biomarker_identifier/blob/main/DGEA/cancer_vs_healthy/2025-08-11_tcga_vs_healthy_deseq.R). Whereas the [2025-08-11_hstages_deseq.R](https://github.com/Leah31115/hcc_biomarker_identifier/blob/main/DGEA/cancer_vs_healthy/2025-08-11_hstages_deseq.R) script generates a PCA plot for healthy Vs cancer with points coloured by stage.

Stages I/II/III/IV Vs healthy scripts are accessible in the [cancer_vs_healthy subfolder](https://github.com/Leah31115/hcc_biomarker_identifier/tree/main/DGEA/cancer_vs_healthy).

The paired samples PCA plot was generated using [2025-06-10_paired_TCGA_deseq.R](https://github.com/Leah31115/hcc_biomarker_identifier/blob/main/DGEA/paired_TCGA/2025-06-10_paired_TCGA_deseq.R) which identified a sample pair needed for removal. Thus the [2025-07-24_paired_tcga_deseq_removed_sample.R](https://github.com/Leah31115/hcc_biomarker_identifier/blob/main/DGEA/paired_TCGA/2025-07-24_paired_tcga_deseq_removed_sample.R) script runs DGEA with the undesired sample pair removed.

Cancer stage (II/III/IV) Vs Stage (IV) and Advanced HCC (Stage II-IV) Vs Stage I scripts are accessible in the [cancer_stage_vs_stage](https://github.com/Leah31115/hcc_biomarker_identifier/tree/main/DGEA/cancer_stage_vs_stage) subfolder. 

Significant and differential genes expressed (DGE) were compared from all DGEAs above and represented by bargraphs using the [2025-08-18_bar_graph_significance.R](https://github.com/Leah31115/hcc_biomarker_identifier/blob/main/DGEA/quantity_of_all_DGEA_significant_genes/2025-08-18_bar_graph_significance.R) script.

## Survival analysis
Survival data for TCGA cancer samples was downloaded from UCSC Xena [here](https://xenabrowser.net/datapages/?dataset=survival%2FLIHC_survival.txt&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) (file LIHC_survival.txt, version 2018-09-13, 23 KB). The normalised TCGA STAR counts were used for data analysis and was downloaded from [here](https://xenabrowser.net/datapages/?dataset=TCGA-LIHC.star_tpm.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) (TCGA-LIHC.star_tpm.tsv.gz, version 05-09-2024, 98.5 MB). All survival analysis scripts are located in the [survival_analysis](https://github.com/Leah31115/hcc_biomarker_identifier/tree/main/survival_analysis) folder. The rstatix, ggpubr, survminer, survival, reshape2, tidyverse, tidytidbits, and survivalAnalysis packages in RStudio were used for survival analysis (Kassambara 2023; 2025; Kassambara, Kosinski, and Biecek 2024; Therneau, Grambsch 2000; Therneau 2024; Wickham 2007; Wickham et al. 2019; Wiesweg 2022; 2025).


### Univariate Kaplan-Meier Curves
Code was provided by Dr Heshmat Borhani to investigate DSS, overall survival and progression free interval (PFI). Scripts to generate survival plots for genes of interest located in the [univariate_KM_curves](https://github.com/Leah31115/hcc_biomarker_identifier/tree/main/survival_analysis/univariate_KM_curves) subfolder.

### Multivariate forest plots
Multivariate analysis of the pseudogenes and their parent protein coding genes were investigated for their association with DSS using forest plots. 

Protein coding genes with multiple transcript isoforms were investigated for their effect on DSS. Healthy and cancer gene isoforms were obtained from the TCGA TARGET GTEx hub via UCSC Xena [here](https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_rsem_isoform_tpm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) (TcgaTargetGtex_rsem_isoform_tpm.gz, version 2016-09-02, 4.16 GB) (Vivian et al., 2017). Since this file was large (GB), the University of Nottingham's High Performance Computer (HPC) Ada was used to filter for used samples via Ubuntu. Conda for linux can be downloaded [here](https://docs.conda.io/projects/conda/en/stable/user-guide/install/linux.html). 
First the used samples file was uploaded from the local computer to the HPC
```bash
scp my_work_dir/used_samples.txt ada:/gpfs01/home/mbyle1/LIFE4137/splice/
```

The following commands were used in the Ubuntu command line to filter the file for desired samples:
```bash
# Move to working directory
cd /gpfs01/home/mbyle1/LIFE4137/splice/

# Download isoform data
wget https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_isoform_tpm.gz ./

# Create conda environment
conda create -n R4.32 -y
conda activate
conda install -n R4.32 -c conda-forge r-textshaping -y

# Filter for desired samples using R
R
install.packages(dplyr)
library(dplyr)
isoform_counts <- read.table(file="/gpfs01/home/mbyle1/LIFE4137/splice/TcgaTargetGtex_rsem_isoform_tpm", header=T, check.names = FALSE, row.names = 1)
used_samples <- read.table(file="/gpfs01/home/mbyle1/LIFE4137/splice/used_samples.txt")
used_samples <- as.vector(used_samples[,1])

isoform_counts <- isoform_counts %>%
select(, all_of(used_samples))

write.table(isoform_counts,
file = "/gpfs01/home/mbyle1/LIFE4137/splice/isoform_counts.txt",
sep = "\t",
row.names = TRUE,
col.names = TRUE
)
```

The filered isoform counts were downloaded to the local computer using the following command
```bash
scp ada:/gpfs01/home/mbyle1/LIFE4137/splice/isoform_counts.txt" ./
```

These isoforms were filtered for genes of interest using this script. The file from the filtered genes script was then used to peform survival analysis of isoforms from genes of interest. These scripts are located in the folder [transcript_isoforms](https://github.com/Leah31115/hcc_biomarker_identifier/tree/main/survival_analysis/multivariate_forest_plots/transcript_isoforms) subfolder within the multivariate_forest_plots folder.

Forest plots were also generated for the pseudogenes identified within the top 5 DEGs from all DGEAs, where any other significant pseudogene relative and its parent gene were included within the analysis. These scripts are located within the [pseudogenes](https://github.com/Leah31115/hcc_biomarker_identifier/tree/main/survival_analysis/multivariate_forest_plots/pseudogenes) subfolder within the multivariate_forest_plots folder.


# Acknowledgements
I would like to thank Dr Heshmat Borhani for providing univariate Kaplan-Meier curve survival analysis code. I thank Thomas Mclaughlin and Areeba Salman for our great teamwork when working on our past group project together in which some of the collaborated code could be used and expanded upon within this independent research project.

# References
- Bache S, Wickham H. 2022. magrittr: A Forward-Pipe Operator for R. doi:10.32614/CRAN.package.magrittr. Available at: https://doi.org/10.32614/CRAN.package.magrittr, R package version 2.0.3 available at https://CRAN.R-project.org/package=magrittr. Accessed 05/08/2025.
- Carlson M. 2025. org.Hs.eg.db: Genome wide annotation for Human. R package version 3.21.0.
- Goldman MJ, Craft B, Hastie M, Repečka K, McDade F, Kamath A, Banerjee, A, Luo Y, Rogers D, Brooks AN, et al. 2020. Visualizing and interpreting cancer genomics data via the Xena platform. Nature Biotechnology, 38(6), pp.675–678. doi: https://doi.org/10.1038/s41587-020-0546-8 .
- Huber W, Heydebreck, AV, Sültmann H, Poustka A, Vingron M. 2002. Variance stabilization applied to microarray data calibration and to the quantification of differential expression. Bioinformatics (Oxford, England), [online] 18 Suppl 1, pp.S96-104. doi:https://doi.org/10.1093/bioinformatics/18.suppl_1.s96.
- Kassambara A. 2023. rstatix: Pipe-Friendly Framework for Basic Statistical Tests. doi:10.32614/CRAN.package.rstatix. Available at: https://doi.org/10.32614/CRAN.package.rstatix. Accessed 18/08/2025. R package version 0.7.2, https://CRAN.R-project.org/package=rstatix. Accessed 18/08/2025.
- Kassambara A. 2025. ggpubr: 'ggplot2' Based Publication Ready Plots. doi:10.32614/CRAN.package.ggpubr. Available at: https://doi.org/10.32614/CRAN.package.ggpubr. Accessed 05/08/2025. R package version 0.6.0 available at https://CRAN.R-project.org/package=ggpubr. Accessed 05/08/2025.
- Kassambara A, Kosinski M, and Biecek P. 2024. survminer: Drawing Survival Curves using 'ggplot2'. doi:10.32614/CRAN.package.survminer. Available at: https://doi.org/10.32614/CRAN.package.survminer. Accessed 18/08/2025. R package version 0.5.0 https://CRAN.R-project.org/package=survminer. Accessed 18/08/2025.
- Kolde, R. 2025. pheatmap: Pretty Heatmaps. doi:10.32614/CRAN.package.pheatmap. Available at: https://doi.org/10.32614/CRAN.package.pheatmap. Accessed 05/08/2025. R package version 1.0.13 available at https://CRAN.R-project.org/package=pheatmap. Accessed 05/08/2025.
- Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), p.550. doi:https://doi.org/10.1186/s13059-014-0550-8.
-  Mousepixels. 2022. Github. sanbomics_scripts. Available at: https://github.com/mousepixels/sanbomics_scripts/blob/main/GO_in_R.Rmd. Accessed 30/08/2025.
- Neuwirth E. 2022. RColorBrewer: ColorBrewer Palettes. doi:10.32614/CRAN.package.RColorBrewer. Available at: https://doi.org/10.32614/CRAN.package.RColorBrewer, R package version 1.1-3 available at https://CRAN.R-project.org/package=RColorBrewer Accessed 05/08/2025.
- Pagès H, Carlson M, Falcon S, Li N. 2025. AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. doi:10.18129/B9.bioc.AnnotationDbi. Available at https://doi.org/10.18129/B9.bioc.AnnotationDbi. R package version 1.71.1 available at https://bioconductor.org/packages/AnnotationDbi. Accessed 01/09/2025.
-  R Core Team. 2025. R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. Available at: https://www.R-project.org [Accessed 05/08/2025].
-  Slowikowski K. 2024. ggrepel: Automatically Position Non-Overlapping Text
Labels with 'ggplot2'. doi:10.32614/CRAN.package.ggrepel. Available at: https://doi.org/10.32614/CRAN.package.ggrepel, R package version 0.9.6 https://CRAN.R-project.org/package=ggrepel. Accessed 05/08/2025.
- Therneau TM, Grambsch PM. 2000. Modeling Survival Data: Extending the Cox Model. Springer, New York. ISBN 0-387-98784-3.
- Therneau T. 2024. A Package for Survival Analysis in R. R package version 3.8-3. Available at: https://CRAN.R-project.org/package=survival. Accessed 05/08/2025.
- Vivian J, Rao AA, Nothaft FA, Ketchum C, Armstrong J, Novak A, Pfeil J, Narkizian J, Deran AD, Musselman-Brown A, et al. 2017. Toil enables reproducible, open source, big biomedical data analyses. Nature Biotechnology, 35(4), pp.314–316. doi:https://doi.org/10.1038/nbt.3772.
- Wickham H. 2007. Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. doi:https://doi.org/10.18637/jss.v021.i12.
- Wickham H. 2016. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York
- Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, et al. 2019. “Welcome to the tidyverse.”, Journal of Open Source Software, *4*(43), 1686. doi:10.21105/joss.01686. Available at: https://doi.org/10.21105/joss.01686. Accessed 05/08/2025.
- Wiesweg M. 2022. tidytidbits: A Collection of Tools and Helpers Extending the Tidyverse. doi:10.32614/CRAN.package.tidytidbits. Available at: https://doi.org/10.32614/CRAN.package.tidytidbits. Accessed 18/08/2025. R package version 0.3.2, https://CRAN.R-project.org/package=tidytidbits. Accessed 18/08/2025.
- Wiesweg M. 2025. survivalAnalysis: High-Level Interface for Survival Analysis and Associated Plots. doi:10.32614/CRAN.package.survivalAnalysis. Available at: https://doi.org/10.32614/CRAN.package.survivalAnalysis. Accessed 18/08/2025. R package version 0.4.0, https://CRAN.R-project.org/package=survivalAnalysis. Accessed 18/08/2025.
- World Health Organization. 2024. Global Cancer Burden growing, Amidst Mounting Need for Services. [online] WHO. Available at: https://www.who.int/news/item/01-02-2024-global-cancer-burden-growing--amidst-mounting-need-for-services. Accessed 06/08/2025.
- Xu S, Hu E, Cai Y, Xie Z, Luo X, Zhan L, Tang W, Wang Q, Liu B, Wang R et al. 2024. Using clusterProfiler to characterize multiomics data. Nature Protocols. [online] doi:https://doi.org/10.1038/s41596-024-01020-z. (11):3292-3320. PMID: 39019974.

