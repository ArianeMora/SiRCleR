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

#' sircleMotifFimo
#'
#' Uses FIMO to predict TF binding and then finds enrichment in clusters.
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
sircleMotifFimo <- function(geneId, fimoFilename, sircleRCMFilename, TFPadjCol, TFValueCol,
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

#' sircleMotif
#'
#' Uses known TF interactions from DoRoThea and then finds enrichment of TF-Target relationships in clusters.
#'
#' @param clusterGeneId Column name with the gene identifier (note this must match the TF ids used in DoRoTHea i.e. gene name)
#' @param doroFile Filename/path of the output from Dorothea
#' @param sircleRCMFilename Filename/path of the output from sircleRCM
#' @param TFPadjCol Column name of the padjusted value in the sircleRCMFilename file (i.e. your RNAseq, or your protein p column)
#' @param TFValueCol Column name of  rna value for the TF (could also use the proteomics logfc)
#' @param pDisplayCutoff Cutoff for the clusters (only returns results less than this, usually 0.05)
#' @param TFInDataset Whether or not you require the TF to be expressed in the dataframe
#' @param doroLevel level(s) from Dorothea (i.e. A - D) see DoRoThea for more info.
#' @param clusters clusters to test for enrichment.
#' @param outputFilename Filename for the output file.
#' @param ouputDir Directory for saving
#' @param plotOn Whether or not to plot the results
#' @return DF of the results
#' @export
sircleMotif <- function(clusterGeneId, doroFile, sircleRCMFilename, TFPadjCol, TFValueCol, TargetPadjCol, TargetValueCol,
                        RegulatoryLabel, TFInDataset=T, doroLevel=c("A"), clusters=c("MDS", "MDE", "TPDE", "TPDS"),
                        outputFilename="SiRCle-motif.csv", outputDir="", plotOn=T) {
  motf <<- import("scimotf")    # Make global
  scimo <- motf$SciMotf_Doro(doro_file=doroFile,  # fimo.tsv (here you can put in the different files 2500bp, 100bp)
                             cluster_file=sircleRCMFilename,  # output from RCM in the first section
                             cluster_gene_id=clusterGeneId,                     # gene identifier we have been using
                             padj_protein=TFPadjCol,                      # p value for the TF (could also use the value from the proteomics)
                             logfc_protein=TFValueCol,                     # rna value for the TF (could also use the proteomics logfc)
                             padj_rna=TargetPadjCol,
                             logfc_rna=TargetValueCol,
                             output_dir=outputDir,                          # output directory
                             tf_in_dataset=TFInDataset  # Whether or not you require the TF to be expressed in the dataframe
  )          # cutoff for the clusters (only returns results less than this)
  m_df <- scimo$run(as.list(doroLevel), rcm_clusters=as.list(clusters))  # Run the motif finding for each cluster
  write.csv(m_df, outputFilename)
  if (plotOn == T) {
    motf$plot_cluster_tf(outputFilename, save_fig=T)
  }
  return(as.data.frame(m_df))
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
