## ---------------------------
##
## Script name: sircleTF
##
## Purpose of script: Runs and saves motif analysis results for SiRCle RCM
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

#' sircleMotif
#'
#' Uses dorothea and viper to predict the activity and activity change of transcription factors.
#'
#' @param geneId Column name with the gene identifier
#' @param fimoFilename Filename/path of the output from fimo (usually fimo.tsv)
#' @param sircleRCMFilename Filename/path of the output from sircleRCM
#' @param TFPadjCol Column name of the padjusted value in the sircleRCMFilename file (i.e. your RNAseq, or your protein p column)
#' @param TFValueCol Column name of  rna value for the TF (could also use the proteomics logfc)
#' @param pDisplayCutoff Cutoff for the clusters (only returns results less than this, usually 0.05)
#' @param TFInDataset Whether or not you require the TF to be expressed in the dataframe
#' @param fimoPcol Column name of p or q value used to filter fimo results(q-value for adjusted or p-value for un-adj)
#' @param regLabels Column name of the regulatory labels
#' @param outputFilename Filename for the output file.
#' @return scimo an instance of the scimo package
#' @export

sircleMotif <- function(geneId, fimoFilename, sircleRCMFilename, TFPadjCol, TFValueCol,
                        pDisplayCutoff=0.05, TFInDataset=T, fimoPcol="p-value", outputDir="",
                        regLabels="RegulatoryLabels", outputFilename="SiRCle-motif.csv") {

  scimo <- scimotf$SciMotf(fimoFilename,  # fimo.tsv (here you can put in the different files 2500bp, 100bp)
                           sircleRCMFilename,  # output from RCM in the first section
                           regLabels,           # Column name of the regulatory labels from the earlier section
                           geneId,                     # gene identifier we have been using
                           TFPadjCol,                      # p value for the TF (could also use the value from the proteomics)
                           TFValueCol,                     # rna value for the TF (could also use the proteomics logfc)
                           outputDir,                          # output directory
                           tf_in_dataset=TFInDataset,  # Whether or not you require the TF to be expressed in the dataframe
                           fimo_pcol=fimoPcol,           # p or q value used to filter fimo results
                           cluster_pcutoff=pDisplayCutoff)          # cutoff for the clusters (only returns results less than this)
  m_df <- scimo$run()  # Run the motif finding for each cluster
  scimo$u$save_df(m_df, outputFilename)
  return(scimo)
}

#' sircleMotifAddInfo
#'
#' Adds info to the sircleMotif dataframe (i.e. TF activtieis or protein values.)
#'
#' @param scimo The return value from sircleMotif (this is an instance of the scimo package.)
#' @param motifDf Dataframe output by the previous analysis
#' @param otherFilename Filename of the file from which the new data come from
#' @param TFNameColumn Column name of the TF within the otherfile (we use the TFs in this column to match with those in motif DF so they need to be in the same format!)
#' @param TFVlaueColumn Column name with the TF values within the other file
#' @param valueLabel Label to append to the TFValue column name so that you know where the source came from.
#' @param outputFilename Filename of the new file to be saved
#' @return dm_df A data frame with the predicted activities and the added column
#' @export

sircleMotifAddInfo <- function(scimo, motifDf, otherFilename, TFNameColumn, TFValueColumn, valueLabel, outputFilename="SiRCle-motif-updated.csv") {

  # Check if they want to add dorothea as well
  # We can also optionally add in the information from dorethea above to get the "predicted" tf change between KO and WT
  dm_df <- scimo$add_tf_predictions(motifDf,                              # output from scimo
                                          otherFilename,  # output from dorethea above
                                          TFNameColumn,                              # column with the TF name
                                          TFValueColumn,                    # column with the TF value we're interested in.
                                          valueLabel,
                                          outputFilename)                # ouput file name
  scimo$u$save_df(dm_df, outputFilename)
  dm_df <- as.data.frame(read.csv(outputFilename))
  return(dm_df)
}
