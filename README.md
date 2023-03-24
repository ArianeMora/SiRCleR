# SiRCleR: Signature Regulatory Clustering package
## Overview
The **S**ignature **R**egulatory **Cl**ust**e**ring (SiRCle) method integrates DNA methylation, RNA-seq and proteomics data at the gene level to deconvolute the association between dysregulation within and across possible regulatory layers (DNA methylation, transcription and/or translation).
Based on logical regulatory rules, `sircleRCM`, genes are grouped in SiRCle clusters based on the layer (DNA methylation, transcription and/or translation) where dys-regulation first occurs.\
Using the output of `sircleRCM`, the SiRCle clusters, one can find the primary biological processes altered by applying Over Representation Analysis (ORA) (**`sircleORA`**) and the drivers behind it using Transcription Factor (TF) analysis (**`sircleTF`**).\
Lastly, to compare patient’s subsets (e.g. based on stage), we found that integrating across the data layers prior to performing differential analysis and biological enrichment better captures the biological signal. Hence we use a variational autoencoder (VAE) to learn gene-wise relationships across the three data layers to obtain an integrated value for each gene (**`sircleVAE`**). Unsing the integarted value we next perform a Mann-Whitney U test to identify genes with a significant integrated difference between the patient's groups.

![](https://www.biorxiv.org/content/biorxiv/early/2022/07/04/2022.07.02.498058/F6.large.jpg?width=800&height=600&carousel=1)

**`sircleRCM` functions create the SiRCle Regulatory Clustering Model (RCM) based on logical regulatory rules, which in turn can be used for further downstream analysis.**

[![PyPI](https://img.shields.io/pypi/v/scircm)](https://pypi.org/project/scircm/)


## Tutorial
See R tutorial in the **`vignettes`** folder, which includes a tutorial and data to try the functions.\
If you want to read more about how SiRCle works, please check out our [preprint](https://www.biorxiv.org/content/10.1101/2022.07.02.498058v1)

## Install
**`SiRCleR`** is an R package.
1. Install Rtools if you haven't done this yet, using the appropriate version (e.g.[windows](https://cran.r-project.org/bin/windows/Rtools/) or [macOS](https://cran.r-project.org/bin/macosx/tools/)).
2. Install the latest development version from GitHub with: **`SiRCleR package`** direcly in R:
    ```
    #install.packages("devtools")
    devtools::install_github("https://github.com/ArianeMora/SiRCleR/tree/v1.0.1")
    library(sircle)
    ```
### Dependencies 
If you are using the visualisations for the over representation analysis you will need to install the following tools and cite them.\
1. CRAN packages
```
install.packages('ggnewscale')
```
2. Biocmanager packages
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("enrichplot")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler")
```
While we have done our best to ensure all the dependencies are documented, if they aren't please let us know and we will try to resolve them.

### Windows specifications
Note if you are running **Windows** you might have an issue with long paths, which you can resolve in the registry on Windows 10:\
`Computer Configuration > Administrative Templates > System > Filesystem > Enable Win32 long paths`\
(If you have a different version of Windows, just google "Long paths fix" and your Windows version)

## Run
See the **vignettes** folder for a proper tutorial with data included that you can run!

#### Making your CpGs map to a single gene version
This is only in python at the moment if you're interested in this post an issue and I'll work on adding this to the R version :) 

## Signature Regulatory Clustering model 

The general table of how we define regulatory clusters.

| Methylation      | RNAseq    | Proteomics | Regulation driver_1          | Regulation driver_2     | Regulation_Grouping1 | Regulation_Grouping2 | Regulation_Grouping3 |
|------------------|-----------|------------|------------------------------|-------------------------|----------------------|----------------------|----------------------|
| Hypermethylation | DOWN      | DOWN       | Methylation increase (MDS)   | None                    | MDS                  | MDS                  | MDS                  |
| Hypermethylation | UP        | DOWN       | mRNA increase (TPDE)         | Protein decrease (TMDS) | TPDE+TMDS            | TPDE+TMDS            | TMDS                 |
| Hypermethylation | UP        | UP         | mRNA increase (TPDE)         | None                    | TPDE                 | TPDE                 | TPDE                 |
| Hypermethylation | DOWN      | UP         | Methylation increase (MDS)   | Protein increase (TMDE) | MDS+TMDE             | TMDE                 | TMDE                 |
| Hypermethylation | No Change | UP         | mRNA increase (TPDE)         | Protein increase (TMDE) | TPDE+TMDE            | TMDE                 | TMDE                 |
| Hypermethylation | No Change | DOWN       | mRNA increase (TPDE)         | Protein decrease (TMDS) | TPDE+TMDS            | TMDS                 | TMDS                 |
| Hypermethylation | UP        | No Change  | mRNA increase (TPDE)         | Protein decrease (TMDS) | TPDE+TMDS            | TPDE+TMDS            | TMDS                 |
| Hypermethylation | DOWN      | No Change  | Methylation increase (MDS)   | Protein increase (TMDE) | MDS+TMDE             | MDS+TMDE             | TMDE                 |
| Hypermethylation | No Change | No Change  | Methylation increase (ncRNA) | None                    | MDS-ncRNA            | MDS_ncRNA            | MDS_ncRNA            |
| Hypomethylation  | DOWN      | DOWN       | mRNA decrease (TPDS)         | None                    | TPDS                 | TPDS                 | TPDS                 |
| Hypomethylation  | UP        | DOWN       | Methylation decrease (MDE)   | Protein decrease (TMDS) | MDE+TMDS             | TMDS                 | TMDS                 |
| Hypomethylation  | UP        | UP         | Methylation decrease (MDE)   | None                    | MDE                  | MDE                  | MDE                  |
| Hypomethylation  | DOWN      | UP         | mRNA decrease (TPDS)         | Protein increase (TMDE) | TPDS+TMDE            | TPDS+TMDE            | TMDE                 |
| Hypomethylation  | No Change | UP         | mRNA decrease (TPDS)         | Protein increase (TMDE) | TPDS+TMDE            | TMDE                 | TMDE                 |
| Hypomethylation  | No Change | DOWN       | mRNA decrease (TPDS)         | Protein decrease (TMDS) | TPDS+TMDS            | TMDS                 | TMDS                 |
| Hypomethylation  | UP        | No Change  | Methylation decrease (MDE)   | Protein decrease (TMDS) | MDE+TMDS             | MDE+TMDS             | TMDS                 |
| Hypomethylation  | DOWN      | No Change  | mRNA decrease (TPDS)         | Protein increase (TMDE) | TPDS+TMDE            | TPDS+TMDE            | TMDE                 |
| Hypomethylation  | No Change | No Change  | Methylation decrease (ncRNA) | None                    | MDE+ncRNA            | MDE_ncRNA            | MDE_ncRNA            |
| No Change        | DOWN      | UP         | mRNA decrease (TPDS)         | Protein increase (TMDE) | TPDS+TMDE            | TPDS+TMDE            | TMDE                 |
| No Change        | UP        | DOWN       | mRNA increase (TPDE)         | Protein decrease (TMDS) | TPDE+TMDS            | TPDE+TMDS            | TMDS                 |
| No Change        | DOWN      | DOWN       | mRNA decrease (TPDS)         | None                    | TPDS                 | TPDS                 | TPDS                 |
| No Change        | UP        | UP         | mRNA increase (TPDE)         | None                    | TPDE                 | TPDE                 | TPDE                 |
| No Change        | No Change | UP         | Protein increase (TMDE)      | None                    | TMDE                 | TMDE                 | TMDE                 |
| No Change        | No Change | DOWN       | Protein decrease (TMDS)      | None                    | TMDS                 | TMDS                 | TMDS                 |
| No Change        | UP        | No Change  | mRNA increase (TPDE)         | Protein decrease (TMDS) | TPDE+TMDS            | TPDE+TMDS            | TMDS                 |
| No Change        | DOWN      | No Change  | mRNA decrease (TPDS)         | Protein increase (TMDE) | TPDS+TMDE            | TPDS+TMDE            | TMDE                 |
| No Change        | No Change | No Change  | NoChange                     | NoChange                | NoChange             | NoChange             | NoChange             |

Please post questions and issues related to sci-rcm on the `Issues`  section of the GitHub repository.



