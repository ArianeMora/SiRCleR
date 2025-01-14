## ---------------------------
##
## Script name: sircleORA
##
## Purpose of script: Runs and saves ORA results for SiRCle RCM
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

##' sircleORA
#'
#' Uses enricher to run ORA on each of the clusters as defined by the regulatory labels using prior knowledge such as gene-sets
#'
#' @param InputData DF with column gene names amd regulatory labels (i.e. sircle clusters). Usually output from SiRCle RCM functions.
#' @param regLabels \emph{Optional: } Column name of the column with the regulatory labels in InputData. \strong{default: "RegulatoryLabels"}
#' @param geneID Provide column name for the gene ID in InputData. Needs to match the the gene ID in the pathway file provided
#' @param PriorKnowledge DF including the prior knowledge (PK), ie. pathway file. File must include column "term" with the pathway or TF-regulon name, column "gene" with the gene name and column "Description" with description that will be depicted on the plots.
#' @param PKName \emph{Optional: } Name of the prior knowledge used, i.e. pathway name such as KEGG \strong{default: ""}
#' @param emptyRegLabel \emph{Optional: } Column name of the empty regulatory group \strong{default: ""}
#' @param RemoveBackgroundGenes\emph{Optional: } If TRUE, genes that fall into background based on the choosen Background method for SiRCle RCM are removed from the universe. \strong{default: "TRUE"}
#' @param minGSSize \emph{Optional: } minimum group size in ORA \strong{default: 10}
#' @param maxGSSize \emph{Optional: } maximum group size in ORA \strong{default: 1000}
#' @param Plot_FileType \emph{Optional: } Output file type for the figures one of: "svg", "pdf", "png" \strong{default: "pdf"}
#' @param Plot_p.adj \emph{Optional: } q value cutoff from ORA that should be plotted \strong{default: 0.2}
#' @param Plot_Percentage \emph{Optional: } Percentage of genes that are detected of a pathway \strong{default: 10}
#' @param OutputFileName \emph{Optional: } Output filename \strong{default: ""}
#'
#' @return DF with the results of the over representation analysis
#' @export

sircleORA <- function(InputData,
                      geneID,
                      regLabels="RegulatoryLabels",
                      PriorKnowledge,
                      PKName="",
                      emptyRegLabel="",
                      RemoveBackgroundGenes="TRUE",
                      Plot_FileType="pdf",
                      minGSSize=10,
                      maxGSSize=1000 ,
                      Plot_p.adj=0.2,
                      Plot_Percentage=10,
                      OutputFileName=''){
  ## ------------ Setup and installs ----------- ##
  packages <- c("clusterProfiler", "enrichplot", "ggupset")

  ## ------------ Create Folder ----------##
  #Create Folder, Safe file :
  SiRCleRCM_results_folder = paste(getwd(), "/SiRCleRCM",  sep="")
  if (!dir.exists(SiRCleRCM_results_folder)){dir.create(SiRCleRCM_results_folder)}#
  SiRCleRCM_ORA_results_folder = paste(SiRCleRCM_results_folder, "/SiRCleRCM_ORA_",Sys.Date(), sep="")
  if (!dir.exists(SiRCleRCM_ORA_results_folder)){dir.create(SiRCleRCM_ORA_results_folder)}  # check and create folder

  ## ------------ Run ----------- ##
  # open the data
  if(RemoveBackgroundGenes=="TRUE"){
    df <- InputData
    df <- subset(df, ! df$BG_Method == "FALSE")
  } else{
    df <- InputData
  }

  #Select universe and SiRCle clusters
  allGenes <- as.character(df[[geneID]]) #
  clusterGenes <- subset(df, ! df[[regLabels]] == emptyRegLabel)
  grps_labels <- unlist(unique(clusterGenes[regLabels]))
  #Load Pathways
  Pathway <- PriorKnowledge
  Term2gene <- Pathway[,c("term", "gene")]# term and geneID (e.g. ENTEREZ, must match geneID)
  term2name <- Pathway[,c("term", "Description")]# geneID and description
  #Add the number of genes present in each pathway
  Pathway$Count <- 1
  Pathway_Mean <- aggregate(Pathway$Count, by=list(term=Pathway$term), FUN=sum)
  names(Pathway_Mean)[names(Pathway_Mean) == "x"] <- "Genes_in_Pathway"
  Pathway <- merge(x= Pathway[,-4], y=Pathway_Mean,by="term", all.x=TRUE)
  #Run ORA
  ResList <- list()
  for(g in grps_labels) {
    grpGenes <- subset(df, df[[regLabels]] == g)
    message( "Number of genes in ",g, ": ", nrow(grpGenes), sep="")

    clusterGo <- clusterProfiler::enricher(gene=as.character(grpGenes[[geneID]]),
                                            pvalueCutoff = 1,
                                            pAdjustMethod = "BH",
                                            universe = allGenes,
                                            minGSSize=minGSSize,
                                            maxGSSize=maxGSSize,
                                            qvalueCutoff = 1,
                                            #gson  = NULL,
                                            TERM2GENE=Term2gene ,
                                            TERM2NAME = term2name)
    #ResList[[paste0(g, "S4")]] <- clusterGo
    clusterGoSummary <- data.frame(clusterGo)
    ResList[[g]] <- clusterGoSummary
    if (!(dim(clusterGoSummary)[1] == 0)){
      #Add pathway information (% of genes in pathway detected)
      clusterGoSummary <- merge(x= clusterGoSummary[,-2], y=Pathway[,-2],by.x="ID",by.y="term", all=TRUE)
      clusterGoSummary$Count[is.na(clusterGoSummary$Count)] <- 0
      clusterGoSummary$Percentage_of_Pathway_detected <-round(((clusterGoSummary$Count/clusterGoSummary$Genes_in_Pathway)*100),digits=2)
      clusterGoSummary <- clusterGoSummary[!duplicated(clusterGoSummary$ID),]
      clusterGoSummary <- clusterGoSummary[order(clusterGoSummary$p.adjust),]
      clusterGoSummary <- clusterGoSummary[,c(1, 12, 2:11, 13:14)]
      #Safe file
      write.csv(clusterGoSummary, paste("SiRCleRCM/SiRCleRCM_ORA_", Sys.Date(), "/", OutputFileName,"ClusterGoSummary_", PKName, '-', g, ".csv", sep=""), row.names = FALSE)#Export the ORA results as .csv
      #Add file to list
      ResList[[g]] <- clusterGoSummary
      #Make Selection of terms that should be displayed on the plots
      clusterGoSummary_Select <- clusterGoSummary %>%
        subset(p.adjust <= Plot_p.adj & Percentage_of_Pathway_detected >= Plot_Percentage)
      rownames(clusterGoSummary_Select)<-clusterGoSummary_Select$ID
      #Make the Plots
      if (!(dim(clusterGoSummary_Select)[1] == 0)){#exclude df's that have no observations
        clusterGo@result <- clusterGoSummary_Select[,c(1:12)]
        #1. Dotplot:
        Dotplot <-  enrichplot::dotplot(clusterGo, showCategory=nrow(clusterGoSummary_Select)) +
          ggplot2::ggtitle(paste("Dotplot Cluster", g, PKName, sep=" "))
        ggplot2::ggsave(file=paste("SiRCleRCM/SiRCleRCM_ORA_", Sys.Date(), "/", OutputFileName,"Dotplot_", g,"_" ,PKName, ".", Plot_FileType, sep=""), plot=Dotplot, width=10, height=8)
        plot(Dotplot)
        #2. Emapplot
        if(nrow(clusterGoSummary_Select)<200){
          x2 <- enrichplot::pairwise_termsim(clusterGo)
          Emapplot <-  enrichplot::emapplot(x2, pie_scale=1.5, layout.params = list(layout = "nicely"), showCategory=nrow(clusterGoSummary_Select))+
            ggplot2::ggtitle(paste("Emapplot Cluster", g, PKName, sep=" "))
          ggplot2::ggsave(file=paste("SiRCleRCM/SiRCleRCM_ORA_", Sys.Date(), "/", OutputFileName,"Emapplot_", g,"_", PKName, ".", Plot_FileType, sep=""), plot=Emapplot, width=10, height=8)
          plot(Emapplot)
        }else{
          message("Emapplot not plotted for ", g,  " due to too many terms", sep="")
        }
      }else{
        message("No terms made the thresholds set for ", g, sep="")
      }
    }
  }
  assign(paste("ResORA_",OutputFileName, sep=""), ResList, envir=.GlobalEnv)
}


