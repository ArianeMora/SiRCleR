# SiRCleR: Signature Regulatory ClusteRing package
**Creates the SiRCle Regulatory Clustering Model (RCM) based on logical regulatory rules, which can be used for downstream analysis:**
* SiRCle Regulatory Clustering Model Visualisation (**RCMVis**)
* Over Representation Analysis (**ORA**) on the individual SiRCle clusters
* Fishers Exact Test (**FET**) for data overlay on the individual SiRCle clusters
* Transcription Factor (**TF**) analysis based on Dorothea regulons or motif analysis on the individual SiRCle clusters and TF visualisation (**TFVis**)
* Variational Autoencoder (**VAE**) statistics

[![PyPI](https://img.shields.io/pypi/v/scircm)](https://pypi.org/project/scircm/)


## Tutorial

See R tutorial in the **vignettes** folder, which includes a tutorial and data to try the functions.

If you want to read more about how SiRCle works, please check out our paper: https://www.biorxiv.org/content/10.1101/2022.07.02.498058v1 

## Install

### R version
First install Rtools if you haven't done this yet. There are different versions (windows: https://cran.r-project.org/bin/windows/Rtools/, macOS: https://cran.r-project.org/bin/macosx/tools/)

### R version

Just install direcly in R:

```
#install.packages("devtools")
library(devtools)
install_github("ArianeMora/SiRCleR")
library(sircle)
```
### Windows 
Note if you are running **Windows** you might have an issue with long paths, which you can resolve in the registry on Windows 10 (Computer Configuration > Administrative Templates > System > Filesystem > Enable Win32 long paths). If you have a different version of Windows, just google "Long paths fix" and your Windows version.

### Other dependencies 
If you are using the visualisations for the over representation analysis you will need to install the following tools and cite them.

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

Also while we have done our best to ensure all the dependencies are documented, if they aren't please let us know! And we will try to resolve them.

## Run
See the **vignettes** folder for a proper tutorial with data included that you can run!

#### Making your CpGs map to a single gene version
This is only in python at the moment if you're interested in this post an issue and I'll work on adding this to the R version :) 

## Signature Regulatory Clustering model 

The general table of how we define regulatory clusters.

Please post questions and issues related to sci-rcm on the `Issues`  section of the GitHub repository.



