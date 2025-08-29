# Hepatocellular carcinoma (HCC) biomarker indentification
## Contents
[Overview](#overview) <br/>
[Getting started with R](#Getting-started-with-R) <br/>
[Filter samples](#filter-samples) <br/>
### Differential gene expression analysis (DGEA)
[DGEA](#Differential-Gene-Expression-Analysis-(DGEA)) <br/>
[Top differential genes expressed (DGE)](#Top-differential-genes-expressed (DGE)) <br/>
### Survival analysis
[Univariate Kaplan-Meier curves](#Univariate-Kaplan-Meier-Curves) <br/>
[Multivariate Forest plots](#Multivariate-forest-plots) <br/>
<br/>
[Acknowledgements](#Acknowledgements) <br/>
[References](#References) <br/>

## Overview
Identifying early hepatocellular carcinoma (HCC) biomarkers for my thesis. LIFE4137 Independent Project Module, University of Nottingham. The code below was used to to investigate early biomarkers in hepatocellular carcinoma (HCC), a common form of liver cancer with high mortality rates. All data required for downstream analysis is open access and accessible via University of California Santa Cruz (UCSC)  [Xena](https://xenabrowser.net/datapages/). 

## Getting started with R
R studio (version 2025.05.1+513) and R (version 4.5.1) can be downloaded [here](https://posit.co/download/rstudio-desktop/).
Please adjust the file paths according to your work directory when using these R scripts for analysis. To install any packages in R use the following command and replace "package_name" accordingly. 
```R
install.packages("package_name")
```
## Filter samples 
RNAseq STAR gene counts for The Cancer Genome Atlas (TCGA) cancer samples from the GDC TCGA Liver Cancer (LIHC) cohort were downloaded from UCSC Xena [here](https://xenabrowser.net/datapages/?dataset=TCGA-LIHC.star_counts.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) (file TCGA-LIHC.star_counts.tsv.gz, version 05-09-2024, 38.9 MB). TCGA cancer phenotype metadata, from the same cohort, was also obtained from the UCSC Xena platform [here](https://xenabrowser.net/datapages/?dataset=TCGA-LIHC.clinical.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) (file TCGA-LIHC.clinical.tsv.gz, version 09-07-2024, 92.5 KB).

Healthy sample RNAseq TOIL RSEM expected counts were obtained from the [GTEx cohort](https://xenabrowser.net/datapages/?dataset=gtex_gene_expected_count&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) in UCSC Xena (file gtex_gene_expected_count.gz, version 2016-05-19, 476 MB). Healthy metadata was also obtained from UCSC Xena [here](https://xenabrowser.net/datapages/?dataset=GTEX_phenotype&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) (file GTEX_phenotype.gz, version 2016-05-18, 80.5 KB). 

## Differential Gene Expression Analysis (DGEA)
For Deseq2, a DESeq Data Set (DDS) object is made which requires a count data file as well as a condition file. The count data files for all DGEAs are generated from #Filter-samples
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

Healthy samples were compared against cancer samples using the script:. To generate a PCA plot for healthy Vs cancer with points coloured by stage this script was used.

Stages I/II/III/IV Vs healthy scripts are accessible in this folder.

Paired TCGA samples used healthy and solid tissue normal (STN) samples. Samples were paired by TCGA ID using the following script. All paired samples were subject to DGEA using this script.

Cancer stage (II/III/IV) Vs Stage (IV) and Advanced HCC (Stage II-IV) Vs Stage I scripts are accessible in this folder. 

Significant and differential genes expressed (DGE) were compared from all DGEAs above and represented by bargraphs using this script

## Top differential genes expressed (DGE)
From all of the DGEAs (here), the top 5 DGEs were filtered for their up/down-regulation and their p-values here. 

## Survival analysis
Survival data for TCGA cancer samples was downloaded from UCSC Xena [here](https://xenabrowser.net/datapages/?dataset=survival%2FLIHC_survival.txt&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) (file LIHC_survival.txt, version 2018-09-13, 23 KB)

### Univariate Kaplan-Meier Curves
Code was provided by Dr Heshmat Borhani to investigate DSS, overall survival and progression free interval (PFI).

### Multivariate forest plots
Multivariate analysis of the pseudogenes and their parent protein coding genes were investigated for their association with DSS using forest plots. 

Protein coding genes with multiple transcript isoforms were investigated for their effect on DSS. Healthy and cancer gene isoforms were obtained from the TCGA TARGET GTEx hub via UCSC Xena [here](https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_rsem_isoform_tpm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) (TcgaTargetGtex_rsem_isoform_tpm.gz, version 2016-09-02, 4.16 GB). Since this file was large (GB), the University of Nottingham's High Performance Computer (HPC) Ada was used to filter for used samples via Ubuntu. Used samples were obtained using the script XXXXXXXXX. Conda for linux can be downloaded [here](https://docs.conda.io/projects/conda/en/stable/user-guide/install/linux.html). 
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

These isoforms were filtered for genes of interest using this script. The file from the filtered genes script was then used to peform survival analysis of isoforms from genes of interest. These scripts are located in this folder.



# Acknowledgements
I would like to thank Dr Heshmat Borhani for providing univariate Kaplan-Meier curve survival analysis code.

# References
