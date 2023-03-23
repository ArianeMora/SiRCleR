## ---------------------------
##
## Script name: sircleRCMVis
##
## Purpose of script: Uses the sircleRCM output to generate visualisations.
##
## Author: Christina Schmidt and Ariane Mora
##
## Date Created: 2021-01-18
##
## Copyright (c) Christina Schmidt and Ariane Mora
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

#' sircleRCM_MRP
#'
#' Computes the regulatory clustering model (RCM) using the SiRCle regulatory rules for DNA-methylation, RNAseq and Proteomics data layers (MRP).
#'
#' @param rnaFile Filename for your RNAseq data (results from DeSeq2 and also your normalised expression counts)
#' @param methFile  Filename for your DNA methylation data (results from differential methylation analysis)
#' @param proteinFile Filename/path of you Protein data (results from DeSeq2 and also your normalised expression counts)
#' @param geneID Column name of geneId this MUST BE THE SAME in each of your protein, RNAseq and DNAmethylation files (we join on this)
#' @param rnaValueCol Column name of RNA value in rnaFile 
#' @param rnaPadjCol Column name of RNA p adjusted value in rnaFile 
#' @param methValueCol Column name of Methylation difference value in methFile 
#' @param methPadjCol Column name of Methylation p adjusted value in methFile 
#' @param proteinValueCol  Column name of protein log fold change in proteinFile
#' @param proteinPadjCol Column name of protein p adjusted value in proteinFile
#' @param rnaPadjCutoff  \emph{Optional: }Padjusted cutoff for RNAseq data \strong{Default=0.05}
#' @param rnaLogFCCutoff \emph{Optional: } LogFoldchange cutoff for RNAseq data \strong{Default=0.5}
#' @param proteinPadjCutoff \emph{Optional: } Padjusted cutoff for Protein data \strong{Default=0.05}
#' @param proteinValueCutoff \emph{Optional: } LogFoldchange cutoff for Protein data \strong{Default=0.3}
#' @param methPadjCutoff \emph{Optional: } Padjusted cutoff for DNA methylation data \strong{Default=0.05}
#' @param methDiffCutoff \emph{Optional: } DNA Methylation difference cutoff for DNA methylation \strong{Default=10}
#' @param backgroundMethod \emph{Optional: } Background method (NEED Description for each one here) \strong{Default="P|(M&R)"}
#' @param outputFileName \emph{Optional: } Output filename \strong{Default=SiRCle_RCM.csv}
#' @return rcm an instance of the rcm package
#' @export
#'

#' sirclePlot
#'
#' Uses scircm to and plots a circle plot of the sircle regulatory groups.
#'
#' @param filename filename of the RCM
#' @param regLabels Labels of the regulatory file
#' @param figType  \emph{Optional: } File ending for the saved figures \strong{Default="pdf"}
#' @return
#' @export
#'
sirclePlot <- function(filename, regLabels="RegulatoryLabels", fileType="pdf") {
  ## ------------ Setup and installs ----------- ##
  packages <- c("packcircles", "tidyr", "ggplot2")
  install.packages(setdiff(packages, rownames(installed.packages())))

  library(packcircles)
  library(tidyr)
  library(ggplot2)

  ## ------------ Run ----------- ##
  scircm_Output <- read.csv(filename)
  #Prepare the Dataframe:
  scircm_Output[,"GeneNumber"]  <- as.numeric(1)
  scircm_Output <- scircm_Output [ , c(regLabels, "GeneNumber")]
  scircm_Output <-aggregate(scircm_Output[,"GeneNumber"], by=list(scircm_Output[[regLabels]]), FUN=sum)
  names(scircm_Output)[names(scircm_Output) == "x"] <- "GeneNumber"
  names(scircm_Output)[names(scircm_Output) == "Group.1"] <- regLabels
  scircm_Output <- subset(scircm_Output, !scircm_Output[[regLabels]] == "None")%>%
    unite(col=ClusterName,c(regLabels, GeneNumber), sep = " ", remove = FALSE, na.rm = FALSE)
  #Prepare the Plot:
  packing <- circleProgressiveLayout(scircm_Output$GeneNumber, sizetype='area')
  DataPlot <-cbind(scircm_Output, packing)
  packing <- circleLayoutVertices(packing, npoints=50)
  #Make the Plot:
  sircle <- ggplot() +
    geom_polygon(data = packing, aes(x, y, group = id, fill=as.factor(id)), colour = "black", alpha = 0.6) +
    ## Add text in the center of each bubble + control its size
    geom_text(data = DataPlot, aes(x, y, size=GeneNumber, label = ClusterName)) +
    scale_size_continuous(range = c(1,4)) +
    ## General theme:
    theme_void() +
    theme(legend.position="none") +
    coord_equal()
  ggsave(file=paste("SiRCle_plot.", fileType, sep="" ), plot=sircle, width=10, height=8)
}

