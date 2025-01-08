
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![GitHub
issues](https://img.shields.io/github/issues/ArianeMora/SiRCleR)](https://github.com/ArianeMora/SiRCleR/issues)
<!-- badges: end -->

# Short Introduction

The **S**ignature **R**egulatory **Cl**ust**e**ring (***SiRCle***)
method has first been released as a biorxiv preprint and since then we
have enhanced the method further as part of our revision in genome
medicine.  
The ***SiRCle*** method was first developed in 2022 to integrate DNA
methylation, RNA-seq and proteomics data at the gene level to
deconvolute the association between dysregulation within and across
possible regulatory layers (DNA methylation, transcription and/or
translation). The results of the study were published on biorxiv in a
shared effort by [Mora & Schmidt et.
al.](https://www.biorxiv.org/content/10.1101/2022.07.02.498058v1). As
part of the revisions in Genome Medicine, we have updated some parts of
the SiRCle clustering method and extended our methods section in the
manuscript. We also added new data and a PanCan analysis. All of this
can be found in [Mora&Schmidt et
al.](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-024-01415-3)
\[@Mora_Schmidt2024\].  You can find the detailed code for both
manuscripts on our
[website](https://arianemora.github.io/SiRCle_multiomics_integration/#)
or in our [SiRCleManuscript GitHub
repository](https://github.com/ArianeMora/SiRCle_multiomics_integration).  
  
![The SiRCle method integrates DNA methylation, RNA-seq and proteomics
data at the gene level to deconvolute the association between
dysregulation within and across possible regulatory layers (DNA
methylation, transcription and/or translation).Based on logical
regulatory rules, `sircleRCM`, genes are grouped in SiRCle clusters
based on the layer (DNA methylation, transcription and/or translation)
where dys-regulation first occurs. Using the output of `sircleRCM`, the
SiRCle clusters, one can find the primary biological processes altered
by applying Over Representation Analysis (ORA) (**`sircleORA`**) and the
drivers behind it using Transcription Factor (TF) analysis
(**`sircleTF`**). Lastly, to compare patient’s subsets (e.g. based on
stage), we found that integrating across the data layers prior to
performing differential analysis and biological enrichment better
captures the biological signal. Hence we use a variational autoencoder
(VAE) to learn gene-wise relationships across the three data layers to
obtain an integrated value for each gene (**`sircleVAE`**). Unsing the
integarted value we next perform a Mann-Whitney U test to identify genes
with a significant integrated difference between the patient’s
groups.](https://www.biorxiv.org/content/biorxiv/early/2022/07/04/2022.07.02.498058/F6.large.jpg?width=800&height=600&carousel=1)  
  
***SiRCle*** uses logical regulatory rules, `sircleRCM` to group genes
in SiRCle clusters based on the data layer (DNA methylation,
transcription and/or translation) where dys-regulation first occurs.  
Using the output of `sircleRCM`, the SiRCle clusters, one can find the
primary biological processes altered by applying Over Representation
Analysis (ORA) and the drivers behind it using Transcription Factor (TF)
analysis.  
Lastly, to compare patient’s subsets (e.g. based on stage), we found
that integrating across the data layers prior to performing differential
analysis and biological enrichment better captures the biological
signal. Hence we use a variational autoencoder (VAE) to learn gene-wise
relationships across the three data layers to obtain an integrated value
for each gene (**`sircleVAE`**). Using the integrated value we next
perform a Mann-Whitney U test to identify genes with a significant
integrated difference between the patient’s groups.  

# Tutorials

We have generated several vignettes showcasing the the usage of SiRCleR
using the publicly available datasets as in the manuscript [Mora&Schmidt
et
al.](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-024-01415-3),
which are included as example data within SiRCleR. You can find those
tutorial on the top under the “Tutorials” button, where you can follow a
specific user case example and learn more about how to choose the
settings including the `Background Method`, `input data thresholds` and
`Regulation Grouping`. You can also follow the links below:  
\*
[ChooseSettings](https://ArianeMora.github.io/SiRCleR/articles/ChooseSettings.html)  
\*
[DataAnalysisExample](https://ArianeMora.github.io/SiRCleR/articles/DataAnalysisExample.html)  
  
If you are seeking to also apply the variational autoencoder on the
SiRCle clustering output, please visit the python package or use
reticulate in R. You can find the python package on the [Python
SiRCle](https://github.com/ArianeMora/scircm).  

# Install

**`SiRCleR`** is an R package. 1. Install Rtools if you haven’t done
this yet, using the appropriate version
(e.g.[windows](https://cran.r-project.org/bin/windows/Rtools/) or
[macOS](https://cran.r-project.org/bin/macosx/tools/)). 2. Install the
latest development version from GitHub with: **`SiRCleR package`**
direcly in R:
`devtools::install_github("https://github.com/ArianeMora/SiRCleR")     library(SiRCleR)`

## Dependencies

While we have done our best to ensure all the dependencies are
documented, if they aren’t please let us know and we will try to resolve
them.

## Windows specifications

Note if you are running **Windows** you might have an issue with long
paths, which you can resolve in the registry on Windows 10:  
`Computer Configuration > Administrative Templates > System > Filesystem > Enable Win32 long paths`  
(If you have a different version of Windows, just google “Long paths
fix” and your Windows version)

## Questions & Issues

Please post questions and issues related to **`SiRCleR`** functions on
the `Issues` section of this GitHub repository.

## Reproducibility

If you want to reproduce the results of of our
[publication](https://doi.org/10.1101/2022.07.02.498058), please use the
python package version found here:
<https://doi.org/10.1101/2022.07.02.498058>

# Citation

If you use this please cite our manuscript [Mora&Schmidt et
al.](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-024-01415-3).
Ariane Mora & Christina Schmidt, Brad Balderson, Christian Frezza &
Mikael Bodén. 2024. “SiRCle (Signature Regulatory Clustering) Model
Integration Reveals Mechanisms of Phenotype Regulation in Renal Cancer.”
Genome Medicine 16 (1): 144.
<https://doi.org/10.1186/s13073-024-01415-3>.
