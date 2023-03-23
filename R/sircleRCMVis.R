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



#' sircleAlluvial
#'
#' This function creates an alluvial plot from an input CSV to create an alluvial plot for the RCM_MRP workflow.
#'
#' @param filename Path to the input file
#' @return alluvialPlt The plot generated
#' @export

sircleAlluvial_MRP <- function(filename, outputFilename="SiRCle_Alluvial.pdf") {
  ## ------------ Setup and installs ----------- ##
  packages <- c("alluvial")
  install.packages(setdiff(packages, rownames(installed.packages())))
  library(alluvial)
  Alluvial_DF <- read.csv(filename)
  Alluvial_DF <-Alluvial_DF[,c(1:3,6,4)]
  Alluvial_DF[,"Frequency"]  <- as.numeric("1")
  #Now we can make the plot. First lets try to use the alluvial command:
  alluvialPlt <- alluvial(Alluvial_DF[,c(1:5)], freq=Alluvial_DF$Frequency,
                          col = case_when(Alluvial_DF$Regulation == "Methylation Driven Supression (Repression of Gene Transcription)"  ~ '#AA0E37',
                                          Alluvial_DF$Regulation == "Transcription and Processing Driven Enhancement and/or Translation and Post-translational Modification Driven Supression" ~ '#5A5A5E',
                                          Alluvial_DF$Regulation == "Transcription and Processing Driven Enhancement" ~ '#034D89',
                                          Alluvial_DF$Regulation == "Translation and post-translational Modification Driven Enhancement" ~ '#067734',
                                          Alluvial_DF$Regulation == "Transcription and Processing Driven Supression" ~ '#4998F3',
                                          Alluvial_DF$Regulation == "Translation and Post-translational Modification Driven Supression" ~ '#46D37F',
                                          Alluvial_DF$Regulation == "Methylation Driven Enhancement (Enhancement of Gene Transcription)" ~ '#ED4367',
                                          Alluvial_DF$Regulation == "Transcription and Processing Driven Supression and/or Translation and post-translational Modification Driven Enhancement" ~ '#D2D2D6',
                                          Alluvial_DF$Regulation == "Methylation Driven Supression (Repression of Gene Transcription) and Translation and Post-translational Modification Driven Enhancement"~ '#FC81DC',
                                          Alluvial_DF$Regulation == "Methylation Driven Enhancement (Enhancement of Gene Transcription) and Translation and Post-translational Modification Driven Supression"~ '#991A7E',
                                          Alluvial_DF$Regulation == "Methylation Driven Supression of non-coding RNA"~ '#FDFB89',
                                          Alluvial_DF$Regulation == "Methylation Driven Enhancement of non-coding RNA"~ '#FCCE3A',
                                          TRUE ~ 'grey'),
                          border = case_when(Alluvial_DF$Regulation == "Methylation Driven Supression (Repression of Gene Transcription)"  ~ '#AA0E37',
                                             Alluvial_DF$Regulation == "Transcription and Processing Driven Enhancement and/or Translation and Post-translational Modification Driven Supression" ~ '#5A5A5E',
                                             Alluvial_DF$Regulation == "Transcription and Processing Driven Enhancement" ~ '#034D89',
                                             Alluvial_DF$Regulation == "Translation and post-translational Modification Driven Enhancement" ~ '#067734',
                                             Alluvial_DF$Regulation == "Transcription and Processing Driven Supression" ~ '#4998F3',
                                             Alluvial_DF$Regulation == "Translation and Post-translational Modification Driven Supression" ~ '#46D37F',
                                             Alluvial_DF$Regulation == "Methylation Driven Enhancement (Enhancement of Gene Transcription)" ~ '#ED4367',
                                             Alluvial_DF$Regulation == "Transcription and Processing Driven Supression and/or Translation and post-translational Modification Driven Enhancement" ~ '#D2D2D6',
                                             Alluvial_DF$Regulation == "Methylation Driven Supression (Repression of Gene Transcription) and Translation and Post-translational Modification Driven Enhancement"~ '#FC81DC',
                                             Alluvial_DF$Regulation == "Methylation Driven Enhancement (Enhancement of Gene Transcription) and Translation and Post-translational Modification Driven Supression"~ '#991A7E',
                                             Alluvial_DF$Regulation == "Methylation Driven Supression of non-coding RNA"~ '#FDFB89',
                                             Alluvial_DF$Regulation == "Methylation Driven Enhancement of non-coding RNA"~ '#FCCE3A',
                                             TRUE ~ 'grey'),
                          hide = Alluvial_DF$Frequency == 0,
                          cex = 0.5,
                          gap.width = 0.7,
                          xw=0.1,
                          cw=0.2
  )
  ggsave(file=outputFilename, plot=alluvialPlt, width=10, height=8)
  return(alluvialPlt)
}

