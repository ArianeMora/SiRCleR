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

  ## ------------ Setup the environment ----------- ##
  if (! is.null(condaEnvName)) {
    use_condaenv(condaEnvName, required = TRUE)
    setupEnv = T
  }
  if (! is.null(envName)) {
    use_virtualenv(envName, required = TRUE)
    setupEnv = T
  }
  if (! is.null(envPath)) {
    use_python(envPath, required = TRUE)
    setupEnv = T
  }
  if (! setupEnv) {
    print("WARNING: NEXT TIME USE envNAME --> we are creating you a new virtual environment called sircle. THIS TAKES TIME, so next time pass envName=sircle to skip this step. ")
    print("WARNING: YOUR SYSTEM NEED PYTHON > 3.6. If it doesn't nothing will work, please update to python 3.6 or consider installing conda!")
    virtualenv_create(
      envname = "sircle",
      python = NULL,
      packages = "scircm",
      system_site_packages = getOption("reticulate.virtualenv.system_site_packages",
                                       default = FALSE)
    )
    use_python(envPath, required = TRUE)
    setupEnv = T
  }
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

#' sircleRCM
#'
#' Uses scircm to compute the regulatory clustering model.
#'
#' @param rnaFile Filename for your RNAseq data (results from DeSeq2 and also your normalised expression counts)
#' @param methFile  Filename for your DNA methylation data (results from differential methylation analysis)
#' @param proteinFile Filename/path of you Protein data (results from DeSeq2 and also your normalised expression counts)
#' @param geneId Column name of geneId this MUST BE THE SAME in each of your protein, RNAseq and DNAmethylation files (we join on this)
#' @param rnaValueCol Column name of RNA value in rnaFile (usually logFoldChange_r)MUST BE UNIQUE BETWEEN methFile and proteinFile
#' @param rnaPadjCol Column name of RNA p adjusted value in rnaFile (usually padj_r) MUST BE UNIQUE BETWEEN methFile and proteinFile
#' @param methValueCol Column name of Methylation difference value in methFile (usually meth.diff) MUST BE UNIQUE BETWEEN rnaFile and proteinFile
#' @param methPadjCol Column name of Methylation p adjusted value in methFile (usually padj_m) MUST BE UNIQUE BETWEEN rnaFile and proteinFile
#' @param proteinValueCol  Column name of protein log fold change in proteinFile (usually logFC_p) MUST BE UNIQUE BETWEEN rnaFile and methFile
#' @param proteinPadjCol Column name of protein p adjusted value in proteinFile (usually padh_p) MUST BE UNIQUE BETWEEN rnaFile and proteinFile
#' @param proteinCols A list of columns that you want to impute missing values for (optional) jsut uses minimum value (for VAE)
#' @param rnaPadjCutoff  \emph{Optional: }Padjusted cutoff for RNAseq data \strong{Default=0.05}
#' @param rnaLogFCCutoff \emph{Optional: } LogFoldchange cutoff for RNAseq data \strong{Default=0.5}
#' @param proteinPadjCutoff \emph{Optional: } Padjusted cutoff for Protein data \strong{Default=0.05}
#' @param proteinValueCutoff \emph{Optional: } LogFoldchange cutoff for Protein data \strong{Default=0.3}
#' @param methPadjCutoff \emph{Optional: } Padjusted cutoff for DNA methylation data \strong{Default=0.05}
#' @param methDiffCutoff \emph{Optional: } DNA Methylation difference cutoff for DNA methylation \strong{Default=10}
#' @param backgroundMethod \emph{Optional: } Background method (NEED Description for each one here) \strong{Default="P|(M&R)"}
#' @param fileSep \emph{Optional: } Separator for files i.e. expecting CSV's however if they are all TSVs you can change this to "\t" \strong{Default=","}
#' @param nonCodingGeneList \emph{Optional: } List of genes that are annotated to be non-coding related \strong{Default=NULL}
#' @param outputFileName \emph{Optional: } Output filename \strong{Default=SiRCle_RCM.csv}
#' @param logfile \emph{Optional: } Name of the logfile \strong{Default=logfileRCM.csv}
#' @param envName \emph{Optional: } Name of your previously setup python virtual environment \strong{Default=NULL}
#' @param condaEnvName \emph{Optional: } Name of your previously setup python conda environment \strong{Default=NULL}
#' @param envPath \emph{Optional: } Path as a string to your previously setup python \strong{Default=NULL}
#' @return rcm an instance of the rcm package
#' @export
#'

sircleRCM <- function(rnaFile, methFile, proteinFile, geneId,
                      rnaValueCol, rnaPadjCol, methValueCol, methPadjCol, proteinValueCol, proteinPadjCol, proteinCols=NULL,
                      rnaPadjCutoff=0.05, rnaLogFCCutoff=0.5, proteinPadjCutoff=0.05, proteinValueCutoff=0.3,
                      methPadjCutoff=0.05, methDiffCutoff=10, backgroundMethod="P|(M&R)", fileSep=",",
                      nonCodingGeneList=NULL, outputFileName="SiRCle_RCM.csv",
                      logfile="logfileRCM.csv", envName=NULL, condaEnvName=NULL, envPath=NULL) {
  setupEnv = F
  ## ------------ Setup and installs ----------- ##
  packages <- c("tidyverse", "reticulate", "dplyr")
  install.packages(setdiff(packages, rownames(installed.packages())))

  library(tidyverse)
  library(dplyr)
  library(reticulate)

  ## ------------ Setup the environment ----------- ##
  if (! is.null(condaEnvName)) {
    use_condaenv(condaEnvName, required = TRUE)
    setupEnv = T
  }
  if (! is.null(envName)) {
    use_virtualenv(envName, required = TRUE)
    setupEnv = T
  }
  if (! is.null(envPath)) {
    use_python(envPath, required = TRUE)
    setupEnv = T
  }
  if (! setupEnv) {
    print("WARNING: NEXT TIME USE envNAME --> we are creating you a new virtual environment called sircle. THIS TAKES TIME, so next time pass envName=sircle to skip this step. ")
    print("WARNING: YOUR SYSTEM NEED PYTHON > 3.6. If it doesn't nothing will work, please update to python 3.6 or consider installing conda!")
    virtualenv_create(
      envname = "sircle",
      python = NULL,
      packages = "scircm",
      system_site_packages = getOption("reticulate.virtualenv.system_site_packages",
                                       default = FALSE)
    )
    use_python(envPath, required = TRUE)
    setupEnv = T
  }
  scimotf <<- import("scimotf")  # Make global
  scircm <<- import("scircm")    # Make global
  scivae <<- import("scivae")    # Make global

  ## ------------ Run the RCM ----------- ##
  rcm <- scircm$SciRCM(methFile, rnaFile, proteinFile,
                       rnaValueCol, rnaPadjCol, methValueCol, methPadjCol,
                       proteinValueCol, proteinPadjCol, geneId, sep=fileSep,
                       rna_padj_cutoff=rnaPadjCutoff, prot_padj_cutoff=proteinPadjCutoff, meth_padj_cutoff=methPadjCutoff,
                       rna_logfc_cutoff=rnaLogFCCutoff, prot_logfc_cutoff=proteinValueCutoff, meth_diff_cutoff=methDiffCutoff,
                       non_coding_genes=nonCodingGeneList, bg_type=backgroundMethod, logfile=logfile
  )
  # Check if the user wants to impute the protein columns
  if (!is.null(proteinCols)) {
    grp_values <- rcm$run(fill_protein=T, protein_cols=unlist(proteinCols))
  } else {
    rcm$run()
  }
  df <- rcm$get_df()
  # This changes it so we can use it in R again
  rcm$save_df(outputFileName)

  return(rcm)
}

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

#' sircleGSEA
#'
#' Uses scircm to compute the regulatory clustering model.
#' @param rcm Instance of RCM after running sircleRCM
#' @param df Dataframe with the results from the sircleRCM
#' @param geneId  Gene ID column name in you dataframe
#' @param numNodes Number of latent nodes in your VAE
#' @param pathway_list_go A list of pathways for the GSEA (see example for details)
#' @param fgseaPCutoff  \emph{Optional: } Padjusted cutoff for GSEA data (the other results aren't saved) \strong{Default=0.2}
#' @param  groupLabels \emph{Optional: } Labels from the groups  \strong{Default=c('MDE', 'MDS', 'TMDE', 'TMDS', 'TPDE', 'TPDE_TMDS', 'TPDS', 'MDS_TMDE', 'TPDS_TMDE', "MDE-ncRNA", "MDS-ncRNA")}
#' @return
#' @export
#'
sircleGSEA <- function(rcm, df, geneId, numNodes, userPathways, fgseaPCutoff=1.0, fileLabel="GSEA",
                       groupLabels=c('MDE', 'MDS', 'TMDE', 'TMDS', 'TPDE', 'TPDE_TMDS', 'TPDS', 'MDS_TMDE', 'TPDS_TMDE', "MDE-ncRNA", "MDS-ncRNA")) {
  ## ------------ Setup and installs ----------- ##
  packages <- c("fgsea")
  install.packages(setdiff(packages, rownames(installed.packages())))
  library(fgsea)

  for (i in groupLabels) {
    for (n in 0:numNodes) {
      rankRegGrp <- rcm$rank_rcm_by_vae(df, i, as.integer(n), geneId)
      vaeRankValues <- rankRegGrp[[paste("node", as.integer(n), sep="_")]]

      names(vaeRankValues) <- as.character(rankRegGrp[[geneId]])
      #run the GSEA analysis
      fgseaVAE <- fgsea(pathways = userPathways, vaeRankValues, nperm=10000)
      fgseaVAE <- fgseaVAE[fgseaVAE$padj <= fgseaPCutoff, ]

      if (nrow(fgseaVAE) > 1) {
        write_csv(fgseaVAE[,c(-8)], paste("SiRCle-vaeRank-nperm10000-",i, "-Node", n + 1, "_", fileLabel,".csv", sep=""))
      }
    }
  }
}

#' sircleFET
#'
#' Uses scircm to compute the regulatory clustering model.
#' @param rcm Instance of RCM after running sircleRCM
#' @param geneId GeneId column name in your dataframe
#' @param allGenes Background of all genes
#' @param genesOfInterest Genes that you are interested in (i.e. from some gene list)
#' @param groupLabels \emph{Optional: } Labels from the groups  \strong{Default=c('MDE', 'MDS', 'TMDE', 'TMDS', 'TPDE', 'TPDE_TMDS', 'TPDS', 'MDS_TMDE', 'TPDS_TMDE', "MDE-ncRNA", "MDS-ncRNA")}
#' @return
#' @export
#'
sircleFET <- function(rcm, allGenes, genesOfInterest, outputFilename="SiRCle-FET.tsv",
                      groupLabels=c('MDE', 'MDS', 'TMDE', 'TMDS', 'TPDE', 'TPDE_TMDS', 'TPDS', 'MDS_TMDE', 'TPDS_TMDE', "MDE-ncRNA", "MDS-ncRNA")) {
  oddsRatios <- list()
  pValues <- list()
  genesWithInBG <- list()
  genesWithInClutser <- list()
  genesWithOUTInBG <- list()
  genesWithOUTInClutser <- list()
  genesIn <- list()

  i <- 1
  for(g in groupLabels) {
    grpGenes <- rcm$get_genes_in_reg_grp(g, geneId)
    hasInCluster <- genesOfInterest[genesOfInterest %in% grpGenes]

    numGenes_WITH_InCluster <- length(hasInCluster)
    numGenes_WITHOUT_InCluster <- length(grpGenes) - numGenes_WITH_InCluster
    numGenes_WITH_InBG <- length(genesOfInterest) - numGenes_WITH_InCluster
    numGenes_WITHOUT_InBG <- length(allGenes) - numGenes_WITH_InBG

    dat <- matrix(c(numGenes_WITH_InCluster, numGenes_WITHOUT_InCluster, numGenes_WITH_InBG, numGenes_WITHOUT_InBG), ncol=2)
    names <- list(c("in cluster", "not in cluster"), c("has value", "does not have value"))
    Convictions <- matrix(dat, nrow = 2, dimnames = names)
    test <- fisher.test(dat)
    # Now we want to save some stuff
    print(paste(g, test$estimate, test$p.value))
    genesWithInBG[i] <- numGenes_WITH_InBG
    genesWithInClutser[i] <- numGenes_WITH_InCluster
    genesWithOUTInBG[i] <- numGenes_WITHOUT_InBG
    genesWithOUTInClutser[i] <- numGenes_WITHOUT_InCluster
    genesIn[i] <- paste(hasInCluster, collapse=", ")  # Keep track of the genes in the cluster
    oddsRatios[i] <- test$estimate
    pValues[i] <- test$p.value
    i <- i + 1

  }

  pAdjValues <- p.adjust(pValues, method="BH")

  FET_DF <- data.frame(groupLabels)
  FET_DF$oddsRatios <- oddsRatios
  FET_DF$padj <- pAdjValues
  FET_DF$pValues <- pValues
  FET_DF$numGenes_WITH_InCluster <- genesWithInClutser
  FET_DF$numGenes_WITH_InBG <- genesWithInBG
  FET_DF$numGenes_WITHOUT_InCluster <- genesWithOUTInClutser
  FET_DF$numGenes_WITHOUT_InBG <- genesWithOUTInBG
  FET_DF$genes <- genesIn
  # Save the results
  rcm$u$save_df(FET_DF, outputFilename, sep="\t")
  FET_DF <- as.data.frame(read_tsv(outputFilename))
  return(FET_DF)
}
