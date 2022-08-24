# sci-RegulatoryClusteringModel
[![PyPI](https://img.shields.io/pypi/v/scircm)](https://pypi.org/project/scircm/)

## Tutorial

See R tutorial in this folder: https://github.com/ArianeMora/scircm/blob/main/examples/

## Install

#### R version

If you don't have conda, you'll need to do the below, first make sure you have reticulate installed. 
```
install.packages('reticulate')
```
Create a new environment and install scircm.
```
virtualenv_create(
      envname = "scircm",
      python = 3.8,
      packages = "scircm",
      system_site_packages = getOption("reticulate.virtualenv.system_site_packages",
                                       default = FALSE)
    )
```

Note we expect python 3.8 so if things don't work first time, check you're running python 3.8 and then try again :) 

Now you can install it in R:

```
#install.packages("devtools")
library(devtools)
install_github("ArianeMora/SiRCleR")
library(sircle)
```

#### Command line version
Optionally create a new conda env in your terminal (Mac and linux):
```
conda create --name scircm python=3.8
conda activate scircm
```
Go to your terminal and type:
``` 
pip install scircm
```

## Run
See the examples folder for a proper tutorial with data included that you can run!

#### Quick version
```
from scircm import SciRCM
# FORMAT must be csv :) 
protFile = f'path to the output from protein differential abundance file'
rnaFile = f'path to the output from differential expression analysis file'
methFile = f'path to the output from methylation DCpG analysis file'
geneId <- 'ensembl_gene_id'

sircleFileName <- "SircleR-RCM.csv"

#logFC_rna = column name in your RNA file that has your RNA logFC (same for the protein and CpG)
#padj_rna = column name in your RNA file that has your padj value (same for protein and CpG)
#NOTE: these need to be unique from one another since we merge the datasets, if they aren't, you need
#to update your csv files.
#Lastly: ensembl_gene_id this is the gene ID column, All must use the same identifier, and this must be
#labelled the same in each file, if it isn't, update your column names before running.

# Run the sircle RCM this may take some time
rcm <- sircleRCM(rnaFile, methFile, protFile, geneId,  "logFC_rna", "padj_rna", "CpG_Beta_diff", "padj_meth", "logFC_protein", "padj_protein",
                 outputFileName = sircleFileName, 
                 condaEnvName="scircm")

# Plot the sircle function
sirclePlot(sircleFileName, regLabels="Regulation_Grouping_3") 
# That DF now has your rcm clustering results, how easy was that :D
```

#### Making your CpGs map to a single gene version
This is only in python at the moment. Please post any issue if you want to do this in R and we'll write a wrapper :)  (https://github.com/ArianeMora/scircm)


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

Please post questions and issues related to sci-rcm on the `Issues <https://github.com/ArianeMora/scircm/issues>`_  section of the GitHub repository.

### Developers

```
bas_scircm <- BasiliskEnvironment(envname="bas_scircm",
    pkgname="sircle",
    packages=c("scircm==1.0.1")
)

res <- basiliskRun(env=bas_scircm, fun=function(args) {
    out <- reticulate::import("bas_scircm")
    # Do something with pandas
    return(some_r_object)
})
```