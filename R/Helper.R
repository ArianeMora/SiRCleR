## ---------------------------
##
## Script name: HelperFunctions
##
## Purpose of script: Metabolomics (raw ion counts) pre-processing, normalisation and outlier detection
##
## Author: Dimitrios Prymidis and Christina Schmidt
##
## Date Created: 2023-06-14
##
## Copyright (c) Dimitrios Prymidis and Christina Schmidt
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------


#' Imports toy data into environment
#'
#' @title Toy Data Import
#' @description Import and process .csv file to create toy data.
#' @importFrom utils read.csv
#' @return A data frame containing the toy data.
#' @export
toy_data <- function() {
  # Read the .csv files
  Methylation <- system.file("data_example", "rna_DE_Stage IV_sircle_renamed-cols.csv", package = "SiRCleR")
  Methylation<- read.csv(Methylation, check.names=FALSE)

  RNA <- system.file("data_example", "rna_DE_Stage IV_sircle_renamed-cols.csv", package = "SiRCleR")
  RNA<- read.csv( RNA, check.names=FALSE)

  Protein <-system.file("data_example", "prot_DE_Stage IV_sircle.csv", package = "SiRCleR")
  Protein<- read.csv( Protein, check.names=FALSE)

  # Return the toy data into environment
  assign("Meth_DE", Methylation, envir=.GlobalEnv)
  assign("RNA_DE", RNA, envir=.GlobalEnv)
  assign("Prot_DE",  Protein, envir=.GlobalEnv)
}



