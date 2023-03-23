## ---------------------------
##
## Script name: sircleRCM
##
## Purpose of script: Performs the main functions of the sircle RCM package.
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



#' sircleVAE
#'
#' Uses scircm to compute the regulatory clustering model.
#' @param rcmFile Output from running sircleRCM
#' @param patientSampleFile File which has columns with the sample information
#' @param methFile File which has the data from the DNA methylation, this should be already filtered to be a single value (e.g. CpG beta value) for each gene, with the gene ID as the rownames of this file, and the columns associated with a specific samples' values
#' @param methSampleFile File with the sample data corresponding to the samples in the methFile (i.e. is the sample tumour or normal, what other meta info you have)
#' @param rnaFile File which has the data from the RNAseq, expect it to be normalised (e.g. log2(TMM + 1), should use the same gene IDs as in the methylation file
#' @param rnaSampleFile File with the sample data corresponding to the samples in the rnaFile (i.e. is the sample tumour or normal, what other meta info you have)
#' @param proteinFile File which has the data from the proteomics, this should be already be normalised and associated to genes, with the gene ID as the rownames of this file, and the columns associated with a specific samples' values
#' @param proteinSampleFile File with the sample data corresponding to the samples in the proteinFile (i.e. is the sample tumour or normal, what other meta info you have)
#' @param geneId  Gene ID column name in you dataframe
#' @param sampleLabelColumn  Column name in each of the sample files which corresponds to the sample label in the respective file (e.g. the label which has the column names in e.g. proteinFile, this has to be the same in all of your sample DFs)
#' @param conditionColumn  Column name in each of the sample files which corresponds to the sample label in the respective file (e.g. the label which has the column names in e.g. proteinFile, this has to be the same in all of your sample DFs)
#' @param patientIdColumn  Column name in each of the sample files which corresponds to the sample label in the respective file (e.g. the label which has the column names in e.g. proteinFile, this has to be the same in all of your sample DFs)
#' @param configVAEFile Filename of the config.json (see scivae for description.)
#' @param RegulatoryLabel Regulatory label of interest (e.g. "Regulation_Grouping_1", "Regulation_Grouping_2", "Regulation_Grouping_3" )
#' @param fromSaved \emph{Optional: } Whether you're reloading from a saved run \strong{Default=F}
#' @param casesForTraining \emph{Optional: } If you're training this is necessary it is a list of case IDs (in the column patientIDcolumn) \strong{Default=NULL}
#' @param outputFolder \emph{Optional: } path to the folder wih the data \strong{Default=NULL}
#' @param runName \emph{Optional: } Name of the run, useful for reloading saved versions \strong{Default="VAE"}
#' @param normalise \emph{Optional: } How to normalise the data \strong{Default="rows"}
#' @param missingMethod \emph{Optional: } How to fill in the missing data. \strong{Default="mean"}
#' @param verbose \emph{Optional: } How to fill in the missing data. \strong{Default=T}
#' @param envName \emph{Optional: } Name of your previously setup python virtual environment \strong{Default=NULL}
#' @param condaEnvName \emph{Optional: } Name of your previously setup python conda environment \strong{Default=NULL}
#' @param envPath \emph{Optional: } Path as a string to your previously setup python \strong{Default=NULL}
#' @return sv a trained model for performing statistics
#' @export
#'
makeSircleStatModel <- function(rcmFile, patientSampleFile,
                                methFile, methSampleFile, rnaFile, rnaSampleFile, proteinFile, proteinSampleFile,
                                sampleLabelColumn, conditionColumn, patientIdColumn, configVAEFile, RegulatoryLabel, fromSaved=F, casesForTraining=NULL,
                                outputFolder="", runName="VAE", normalise="rows",
                                missingMethod="mean", verbose=T, outputFileName="SiRCle_RCM.csv",
                                logfile="logfileRCM.csv", envName=NULL, condaEnvName=NULL, envPath=NULL) {
  setupEnv = F
  ## ------------ Setup and installs ----------- ##
  packages <- c("tidyverse", "reticulate", "dplyr")
  install.packages(setdiff(packages, rownames(installed.packages())))

  library(tidyverse)
  library(dplyr)
  library(reticulate)
  scircm <<- import("scircm")    # Make global

  ## ------------ Run the RCM Stats ----------- ##
  sv <- scircm$RCMStats(rcmFile,
                        patientSampleFile,
                        methFile,
                        methSampleFile,
                        rnaFile,
                        rnaSampleFile,
                        proteinFile,
                        proteinSampleFile,
                        outputFolder,
                        conditionColumn,
                        sampleLabelColumn,
                        patientIdColumn,
                        configVAEFile, RegulatoryLabel,
                        runName, normalise, verbose, missingMethod)

  if (fromSaved) {
    ## ------------ Re-load saved VAEs ----------- ##
    sv$load_saved_vaes()
    sv$load_saved_encodings(paste0(sv$output_folder, 'encoded_df_', sv$run_name, '.csv'))
    sv$load_saved_inputs(paste0(sv$output_folder, 'vae_input_df_', sv$run_name, '.csv'))
    sv$load_saved_raws(paste0(sv$output_folder, 'raw_input_df_', sv$run_name, '.csv'))
  } else {
    sv$train_vae(casesForTraining) #ToDo: allow users to set their own config. cases=matching_cases, config=config)
    sv$save()  # Save the information we have generated.
  }

  return(sv)
}


#' runSircleStats
#'
#' Runs stats on two groups of patients by integrating their data and computing stats on the latent dimensions.
#'
#' @param sv a pretrained instance of the stats model
#' @param condLabel The condition label (e.g. column) that exists in your patient file and the sample files (e.g. gender)
#' @param cond0 Condition 0 (i.e, female)
#' @param cond1 e.g. your condition 1 value (e.g. male)
#' @return df: a DF with the stats in it
#' @export
#'
runSircleStats <- function(sv, condLabel, cond0, cond1) {
  df <- sv$run_vae_stats(condLabel, cond0, cond1)
  df <- as.data.frame(df)
  return(df)
}

#' sircleVAE
#'
#' Uses scircm to compute the regulatory clustering model.
#' @param rcm Instance of RCM after running sircleRCM
#' @param df Dataframe with the values we want as our input to the VAE (these should all be numeric except for the gene ID column!)
#' @param geneId  Gene ID column name in you dataframe
#' @param configVAEFile Filename of the config.json (see scivae for description.)
#' @param plotVAE  \emph{Optional: } Whether or not you want to plot (and save) the figures for the VAE \strong{Default=T}
#' @param plottingCols  \emph{Optional: } Columns for plotting in your dataset - if you don't select these, and you choose to plot all your columns will be used! \strong{Default=NULL}
#' @param figType  \emph{Optional: } File ending for the saved figures \strong{Default="pdf"}
#' @return rcm an instance of the rcm package
#' @export
#'
sircleVAE <- function(rcm, df, geneId, configVAEFile, plotVAE=T, plottingCols=NULL, figType="pdf") {
  # Now actually run the VAE
  rcm$compute_vae(df, geneId, configVAEFile)

  if (plotVAE) {
    rcmMarkers <- list(rcm$get_genes_in_reg_grp("MDS"),
                       rcm$get_genes_in_reg_grp("TPDE_TMDS"),
                       rcm$get_genes_in_reg_grp("TPDE"),
                       rcm$get_genes_in_reg_grp("TMDE"),
                       rcm$get_genes_in_reg_grp("TPDS_TMDE"),
                       rcm$get_genes_in_reg_grp("TPDS"),
                       rcm$get_genes_in_reg_grp("TMDS"),
                       rcm$get_genes_in_reg_grp("MDE"),
                       rcm$get_genes_in_reg_grp("MDE_TMDS"),
                       rcm$get_genes_in_reg_grp("MDS_TMDE")
    )
    rcmLabels <- c("MDS", "TPDE_TMDS", "TPDE", "TMDE", "TPDS_TMDE", "TPDS", "TMDS", "MDE", "MDE_TMDS")
    vaeVis <- scivae$Vis(rcm$vae, NULL, NULL)
    vaeVis$plot_values_on_scatters(df, geneId, rcmLabels, rcmMarkers, fig_type=figType) # Change to PDF if you want a PDF
    vaeVis$plot_feature_scatters(df, geneId, columns=plottingCols, fig_type=figType)
    rcm$plot_venn(fig_type=figType)
    vaeVis$plot_node_hists(fig_type=figType)
    vaeVis$plot_node_feature_correlation(df, geneId, columns=plottingCols, fig_type=figType)
    vaeVis$plot_node_correlation(fig_type=figType)
    vaeVis$plot_feature_correlation(df, geneId, columns=plottingCols, fig_type=figType)
    vaeVis$plot_input_distribution(df, geneId, columns=plottingCols, fig_type=figType)
  }

}
