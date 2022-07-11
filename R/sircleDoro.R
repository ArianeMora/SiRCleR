## ---------------------------
##
## Script name: sircleDoro
##
## Purpose of script: Runs and saves Dorothea results for SiRCle RCM
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

#' sircleDoro
#'
#' Uses dorothea and viper to predict the activity and activity change of transcription factors.
#'
#' @param filename Path to the input file
#' @param regulons Regulons from dorothea)
#' @param matRNAseq Matrix of RNAseq values
#' @param calcDiffBetweenCond Boolean to determine whether or not to calculate the difference in predicted activities
#' @param cond1 A condition that we can filter the columns of the matrix by so the column must contain the value in Condition 1 (i.e. "KO")
#' @param cond2 A condition that we can filter the columns of the matrix by so the column must contain the value in Condition 2 (i.e. "WT")
#' @return tfDf A data frame with the predicted activities
#' @export

sircleDoro <- function(regulons, matRNAseq, calcDiffBetweenCond=F, cond1=NULL, cond2=NULL, method="rank", outputFileName="SiRCle-dorethea_pred_tfactivities.csv") {
  ## ------------ Setup and installs ----------- ##
  packages <- c("dorothea", "viper", "dplyr")
  install.packages(setdiff(packages, rownames(installed.packages())))
  library(dorothea)
  library(viper)

  ## ------------ Run ----------- ##
  # Compute the TF activities
  tf_activities <- run_viper(matRNAseq, regulons, options = list(method = method,
                                                                 minsize = 1,
                                                                 nes = TRUE,
                                                                 eset.filter = FALSE,
                                                                 cores = 1,
                                                                 verbose = FALSE))
  tfDf <- as.data.frame(tf_activities)
  # compute the mean change between the TF activities

  if (calcDiffBetweenCond) {

    meanEv <- rowMeans(dplyr::select(tfDf, matches(cond1)))
    meanKO <- rowMeans(dplyr::select(tfDf, matches(cond2)))
    meanTFChange <- meanKO - meanEv # The mean change in activity as predicted by doro
    tfDf$meanTFChange <- meanTFChange
  }

  tfDf$TF <- rownames(tfDf)
  write.csv(tfDf, outputFileName, row.names = FALSE)
  return(tfDf)
}
