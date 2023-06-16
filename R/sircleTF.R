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


#' sircleORA_TF
#'
#' Uses enricher to run ORA on each of the clusters as defined by the regulatory labels using a TF regulon or FIMO input of choice.
#'
#' @param filename DF input file (= Output of rcm function).
#' @param regLabels \emph{Optional: } regLabels The label of the column with the regulatory labels \strong{default: "RegulatoryLabels"}
#' @param emptyRegLabel \emph{Optional: } emptyRegLabel The label of the empty regulatory group \strong{default: ""}
#' @param RemoveBackgroundGenes\emph{Optional: } If TRUE, genes that fall into background based on the chosen Background method for SiRCle RCM are removed from the universe. \strong{default: "TRUE"}
#' @param enricher_geneID Provide column name for the gene ID. Needs to match the the gene ID in the pathway file provided
#' @param enricher_Pathways Provide TF regulon file. TF regulon file must include column "term" with the TF name, column "gene" with the target gene name and column "Description" with e.g. TF description that will be depicted on the plots.
#' @param enricher_PathwayName \emph{Optional: } Name of the pathway list used \strong{default: ""}
#' @param fileType \emph{Optional: } fileType Output file type for the figures one of: "svg", "pdf", "png" \strong{default: "pdf"}
#' @param minGSSize \emph{Optional: } minimum group size in ORA \strong{default: 10}
#' @param maxGSSize \emph{Optional: } maximum group size in ORA \strong{default: 1000}
#' @param Plot_p.adj \emph{Optional: } q value cutoff from ORA that should be plotted \strong{default: 0.2}
#' @param Plot_Percentage \emph{Optional: } Percentage of genes that are detected of a TF regulon \strong{default: 10}
#'
#' @return sircleTF an instance of the sircle package
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
    clusterGoSummary <- data.frame(clusterGo)

    if (!(dim(clusterGoSummary)[1] == 0)){
      clusterGoSummary <- clusterGoSummary%>%
      dplyr::rename("TF"="ID",
                    "DetectedTargets"="geneID")
      #Add pathway information (% of genes in pathway detected)
      clusterGoSummary <- merge(x= clusterGoSummary[,-2], y=Pathway[,-2],by.x="TF",by.y="term", all=TRUE)
      clusterGoSummary$Count[is.na(clusterGoSummary$Count)] <- 0
      clusterGoSummary$Percentage_of_TFtargets_detected <-round(((clusterGoSummary$Count/clusterGoSummary$TF_targets)*100),digits=2)
      clusterGoSummary <- clusterGoSummary[!duplicated(clusterGoSummary$TF),]
      clusterGoSummary <- clusterGoSummary[order(clusterGoSummary$p.adjust),]
      clusterGoSummary <- clusterGoSummary[,c(1,9,2:8, 10:11)]
      #Safe file
      write_csv(clusterGoSummary, paste(outputFolder, 'ClusterGoSummary_', enricher_PathwayName, '-', g, '.csv', sep=""))#Export the ORA results as .csv
      #Make Selection of terms that should be displayed on the plots
      clusterGoSummary_Select <- clusterGoSummary %>%
        subset(p.adjust <= Plot_p.adj & Percentage_of_TFtargets_detected >= Plot_Percentage)
      rownames(clusterGoSummary_Select)<-clusterGoSummary_Select$TF
      clusterGoSummary_Select <-clusterGoSummary_Select%>%
        dplyr::rename("ID"="TF",
                    "geneID"="DetectedTargets")
      clusterGoSummary_Select$Description <- clusterGoSummary_Select$ID
      clusterGoSummary_Select<-clusterGoSummary_Select[,c(1,12,2:9)]
      #Make the Plots
      if (!(dim(clusterGoSummary_Select)[1] == 0)) {#exclude df's that have no observations
        clusterGo@result <- clusterGoSummary_Select
        #1. Dotplot:
        Dotplot <-  enrichplot::dotplot(clusterGo, showCategory=nrow(clusterGoSummary_Select)) +
          ggtitle(paste("Dotplot ", g, enricher_PathwayName, sep=" "))
        ggsave(file=paste(outputFolder, "SiRCle-ORA_Dotplot_", g,"_" ,enricher_PathwayName, ".", fileType, sep=""), plot=Dotplot, width=10, height=8)
        plot(Dotplot)
        #2. Emapplot
        x2 <- enrichplot::pairwise_termsim(clusterGo)
        Emapplot <-  enrichplot::emapplot(x2, pie_scale=1.5, layout = "nicely", showCategory=nrow(clusterGoSummary_Select))+
          ggtitle(paste("Emapplot", g, enricher_PathwayName, sep=" "))
        ggsave(file=paste(outputFolder, "SiRCle-ORA_Emapplot_", g,"_", enricher_PathwayName, ".", fileType, sep=""), plot=Emapplot, width=10, height=8)
         plot(Emapplot)
        #4. Upsetplot:
        UpsetPlot <-  enrichplot::upsetplot(clusterGo, showCategory=nrow(clusterGoSummary_Select))+
          ggtitle(paste("UpsetPlot", g, enricher_PathwayName, sep=" "))
        ggsave(file=paste(outputFolder, "SiRCle-ORA_UpsetPlot_", g,"_", enricher_PathwayName, ".", fileType, sep=""), plot=UpsetPlot, width=10, height=8)
        plot(UpsetPlot)
      }
    }
  }
}
