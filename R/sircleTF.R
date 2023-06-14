## ---------------------------
##
## Script name: sircleTF
##
## Purpose of script: Runs and saves transcription factor (TF) analysis based on motifs or dorothea regulons on the results of SiRCle RCM.
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
#' 
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



#' sircleORA_TF
#'
#' Uses enricher to run ORA on each of the clusters as defined by the regulatory labels using a TF regulon of choice
#'
#' @param filename Path to the input file
#' @param regLabels \emph{Optional: } regLabels The label of the column with the regulatory labels \strong{default: "RegulatoryLabels"}
#' @param emptyRegLabel \emph{Optional: } emptyRegLabel The label of the empty regulatory group \strong{default: ""}
#' @paramRemoveBackgroundGenes\emph{Optional: } If TRUE, genes that fall into background based on the choosen Background method for SiRCle RCM are removed from the universe. \strong{default: "TRUE"}
#' @param enricher_geneID Provide column name for the gene ID. Needs to match the the gene ID in the pathway file provided
#' @param enricher_Pathways Provide TF regulon file. TF regulon file must include column "term" with the TF name, column "gene" with the target gene name and column "Description" with e.g. TF description that will be depicted on the plots.
#' @param enricher_PathwayName \emph{Optional: } Name of the pathway list used \strong{default: ""}
#' @param fileType \emph{Optional: } fileType Output file type for the figures one of: "svg", "pdf", "png" \strong{default: "pdf"}
#' @param minGSSize \emph{Optional: } minimum group size in ORA \strong{default: 10}
#' @param maxGSSize \emph{Optional: } maximum group size in ORA \strong{default: 1000}
#' @param Plot_p.adj \emph{Optional: } q value cutoff from ORA that should be plotted \strong{default: 0.2}
#' @param Plot_Percentage \emph{Optional: } Percentage of genes that are detected of a pathway \strong{default: 10}
#'
#' @return
#' @export

sircleORA_TF <- function(filename, regLabels="RegulatoryLabels", emptyRegLabel="", RemoveBackgroundGenes="TRUE", enricher_geneID, enricher_Pathways, enricher_PathwayName="", fileType="pdf", minGSSize=10, maxGSSize=1000 , Plot_p.adj=0.2, Plot_Percentage=10, outputFolder=''){
  ## ------------ Setup and installs ----------- ##
  packages <- c("clusterProfiler", "enrichplot", "ggupset")
  install.packages(setdiff(packages, rownames(installed.packages())))
  ## ------------ Run ----------- ##
  # open the data
  if(RemoveBackgroundGenes=="TRUE"){
    df <- read.csv(filename)
    df <- subset(df, ! df$BG_Method == "FALSE")
  } else{
    df <- read.csv(filename)
  }
  
  #Select universe and SiRCle clusters
  allGenes <- as.character(df[[enricher_geneID]]) #
  clusterGenes <- subset(df, ! df[[regLabels]] == emptyRegLabel)
  grps_labels <- unlist(unique(clusterGenes[regLabels]))
  #Load Pathways
  Pathway <- enricher_Pathways
  Term2gene <- Pathway[,c("term", "gene")]# term and geneID (e.g. ENTEREZ, must match enricher_geneID)
  term2name <- Pathway[,c("term", "Description")]# geneID and description
  #Add the number of genes present in each pathway
  Pathway$Count <- 1
  Pathway_Mean <- aggregate(Pathway$Count, by=list(term=Pathway$term), FUN=sum)
  names(Pathway_Mean)[names(Pathway_Mean) == "x"] <- "TF_targets"
  Pathway <- merge(x= Pathway[,-4], y=Pathway_Mean,by="term", all.x=TRUE)
  #Run ORA
  for(g in grps_labels) {
    grpGenes <- subset(df, df[[regLabels]] == g)
    print(g)
    print(dim(grpGenes))
    clusterGo <- clusterProfiler::enricher(gene=as.character(grpGenes[[enricher_geneID]]),
                                            pvalueCutoff = 1,
                                            pAdjustMethod = "BH",
                                            universe = allGenes,
                                            minGSSize=minGSSize,
                                            maxGSSize=maxGSSize,
                                            qvalueCutoff = 1,
                                            #gson  = NULL,
                                            TERM2GENE=Term2gene ,
                                            TERM2NAME = term2name)
    clusterGoSummary <- data.frame(clusterGo)%>%
      dplyr::rename("TF"="ID",
                    "DetectedTargets"="geneID")
    
    if (!(dim(clusterGoSummary)[1] == 0)){
      #Add pathway information (% of genes in pathway detected)
      clusterGoSummary <- merge(x= clusterGoSummary[,-2], y=Pathway[,-2],by.x="TF",by.y="term", all=TRUE)
      clusterGoSummary$Count[is.na(clusterGoSummary$Count)] <- 0
      clusterGoSummary$Percentage_of_TFtargets_detected <-round(((clusterGoSummary$Count/clusterGoSummary$TF_targets)*100),digits=2)
      clusterGoSummary <- clusterGoSummary[!duplicated(clusterGoSummary$ID),]
      clusterGoSummary <- clusterGoSummary[order(clusterGoSummary$p.adjust),]
      clusterGoSummary <- clusterGoSummary[,c(1,9,2:8, 10:11)]
      #Safe file
      write_csv(clusterGoSummary, paste(outputFolder, 'ClusterGoSummary_', enricher_PathwayName, '-', g, '.csv', sep=""))#Export the ORA results as .csv
      #Make Selection of terms that should be displayed on the plots
      clusterGoSummary_Select <- clusterGoSummary %>%
        subset(p.adjust <= Plot_p.adj & Percentage_of_TFtargets_detected >= Plot_Percentage)
      rownames(clusterGoSummary_Select)<-clusterGoSummary_Select$ID
      #Make the Plots
      if (!(dim(clusterGoSummary_Select)[1] == 0)) {#exclude df's that have no observations
        clusterGo@result <- clusterGoSummary_Select[,1:9]
        #1. Dotplot:
        Dotplot <-  enrichplot::dotplot(clusterGo, showCategory=nrow(clusterGoSummary_Select)) +
          ggtitle(paste("Dotplot ", g, enricher_PathwayName, sep=" "))
        ggsave(file=paste(outputFolder, "SiRCle-ORA_Dotplot_", g,"_" ,enricher_PathwayName, ".", fileType, sep=""), plot=Dotplot, width=10, height=8)
        #2. Emapplot
        x2 <- enrichplot::pairwise_termsim(clusterGo)
        Emapplot <-  enrichplot::emapplot(x2, pie_scale=1.5, layout = "nicely", showCategory=nrow(clusterGoSummary_Select))+
          ggtitle(paste("Emapplot", g, enricher_PathwayName, sep=" "))
        ggsave(file=paste(outputFolder, "SiRCle-ORA_Emapplot_", g,"_", enricher_PathwayName, ".", fileType, sep=""), plot=Emapplot, width=10, height=8)
        #4. Upsetplot:
        UpsetPlot <-  enrichplot::upsetplot(clusterGo, showCategory=nrow(clusterGoSummary_Select))+
          ggtitle(paste("UpsetPlot", g, enricher_PathwayName, sep=" "))
        ggsave(file=paste(outputFolder, "SiRCle-ORA_UpsetPlot_", g,"_", enricher_PathwayName, ".", fileType, sep=""), plot=UpsetPlot, width=10, height=8)
      }
    }
  }
}
