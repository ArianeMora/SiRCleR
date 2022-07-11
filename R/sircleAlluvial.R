## ---------------------------
##
## Script name: sircleAlluvial
##
## Purpose of script: Creates an alluvial plot for the RCM workflow.
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



#' sircleAlluvial
#'
#' This function creates an alluvial plot from an input CSV.
#'
#' @param filename Path to the input file
#' @return alluvialPlt The plot generated
#' @export

sircleAlluvial <- function(filename, outputFilename="SiRCle_Alluvial.pdf") {
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
