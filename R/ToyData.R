## ---------------------------
##
## Script name: HelperFunctions
##
## Purpose of script: General helper functions to check function input and save results
##
## Author: Christina Schmidt
##
## Date Created: 2023-06-14
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

################################################################################################
### ### ### Example Data ### ### ###
################################################################################################

#' Access built-in example data
#'
#' @param Dataset Character: name of a built-in dataset:
#'     \itemize{
#'         \item{\code{"DNA_Methylation"}: }
#'         \item{\code{"RNAseq"}: }
#'         \item{\code{"Proteomics"}: }
#'         \item{\code{"Hallmarks"}: }
#'     }
#'
#' @return A data frame containing the toy data.
#'
#' @description Import and process .csv file to create toy data DF.
#'
#' @examples
#' Proteomics <- SiRCleR::ToyData("Proteomics")
#'
#' @importFrom readr read_csv cols
#' @importFrom magrittr %>% extract2
#' @importFrom tibble column_to_rownames
#' @importFrom logger log_trace
#'
#' @export
#'
ToyData <- function(Dataset) {
  #Available Datasets:
  datasets <- list(
    DNA_Methylation = "DNAMethylation.csv",
    RNAseq = "RNAseq.csv",
    Proteomics = "Proteomics.csv",
    Hallmarks = "Hallmarks.csv"
  )

  rncols <- c("Code", "Metabolite")

  #Load dataset:
  if (!Dataset %in% names(datasets)) {
    message <- sprintf("No such dataset: `%s`. Available datasets: %s", Dataset, paste(names(datasets), collapse = ", "))
    logger::log_trace(paste("Error ", message, sep=""))
    stop(message)
  }

  datasets %>%
  magrittr::extract2(Dataset) %>%
    system.file("data/", ., package = "SiRCleR") %>%
    readr::read_csv(col_types = readr::cols())
}

