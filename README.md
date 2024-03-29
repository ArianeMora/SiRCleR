
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![GitHub
issues](https://img.shields.io/github/issues/ArianeMora/SiRCleR)](https://github.com/ArianeMora/SiRCleR/issues)
<!-- badges: end -->

## Overview

The **S**ignature **R**egulatory **Cl**ust**e**ring (SiRCle) method
integrates DNA methylation, RNA-seq and proteomics data at the gene
level to deconvolute the association between dysregulation within and
across possible regulatory layers (DNA methylation, transcription and/or
translation). Based on logical regulatory rules, `sircleRCM`, genes are
grouped in SiRCle clusters based on the layer (DNA methylation,
transcription and/or translation) where dys-regulation first occurs.  
Using the output of `sircleRCM`, the SiRCle clusters, one can find the
primary biological processes altered by applying Over Representation
Analysis (ORA) (**`sircleORA`**) and the drivers behind it using
Transcription Factor (TF) analysis (**`sircleTF`**).  
Lastly, to compare patient’s subsets (e.g. based on stage), we found
that integrating across the data layers prior to performing differential
analysis and biological enrichment better captures the biological
signal. Hence we use a variational autoencoder (VAE) to learn gene-wise
relationships across the three data layers to obtain an integrated value
for each gene (**`sircleVAE`**). Unsing the integarted value we next
perform a Mann-Whitney U test to identify genes with a significant
integrated difference between the patient’s groups.  
  
**`sircleRCM` functions create the SiRCle Regulatory Clustering Model
(RCM) based on logical regulatory rules, which in turn can be used for
further downstream analysis:**

- Over Representation Analysis (**`sircleORA`**) on the individual
  SiRCle clusters using GO-term pathways or your own pathway list of
  choice, which returns the results DFs and Visualisations of the
  results
- Transcription Factor (**`sircleTF`**) analysis using over
  representation analysis (ORA) on the individual SiRCle clusters using
  amy TF regulon input (e.g. dorothea or FIMO) and returns the results
  and TF visualisation
- Variational Autoencoder statistics currently only runs with the python
  package
  [![PyPI](https://img.shields.io/pypi/v/scircm)](https://pypi.org/project/scircm/)

## Tutorial

See R tutorial in the **`vignettes`** folder, which includes a tutorial
and data to try the functions.  
If you want to read more about how SiRCle works, please check out our
[preprint](https://www.biorxiv.org/content/10.1101/2022.07.02.498058v1)

## Install

**`SiRCleR`** is an R package. 1. Install Rtools if you haven’t done
this yet, using the appropriate version
(e.g.[windows](https://cran.r-project.org/bin/windows/Rtools/) or
[macOS](https://cran.r-project.org/bin/macosx/tools/)). 2. Install the
latest development version from GitHub with: **`SiRCleR package`**
direcly in R:
`devtools::install_github("https://github.com/ArianeMora/SiRCleR")     library(SiRCleR)`
\### Dependencies If you are using the visualisations for the over
representation analysis you will need to install the following tools and
cite them.  
1. CRAN packages

    install.packages('ggnewscale')

2.  Biocmanager packages

<!-- -->

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("enrichplot")
    BiocManager::install("org.Hs.eg.db")
    BiocManager::install("clusterProfiler")

While we have done our best to ensure all the dependencies are
documented, if they aren’t please let us know and we will try to resolve
them.

### Windows specifications

Note if you are running **Windows** you might have an issue with long
paths, which you can resolve in the registry on Windows 10:  
`Computer Configuration > Administrative Templates > System > Filesystem > Enable Win32 long paths`  
(If you have a different version of Windows, just google “Long paths
fix” and your Windows version)

## Run

See the **vignettes** folder for a proper tutorial including data you
can run the functions with. FYI: currently under development!

## Questions & Issues

Please post questions and issues related to **`SiRCleR`** functions on
the `Issues` section of this GitHub repository.

## Reproducibility

If you want to reproduce the results of of our
[publication](https://doi.org/10.1101/2022.07.02.498058), please use the
python package version found here:
<https://doi.org/10.1101/2022.07.02.498058>

## Citation

Ariane Mora, Christina Schmidt, Brad Balderson, Christian Frezza &
Mikael Bodén 2022. SiRCle (Signature Regulatory Clustering) model
integration reveals mechanisms of phenotype regulation in renal cancer.
BioRxiv. <https://doi.org/10.1101/2022.07.02.498058>
