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

#' sircleORA
#'
#' Uses encrichGo to run ORA on each of the clusters as defined by the regulatory labels using GO terms pathways
#'
#' @param filename Input file
#' @param entrezId Column name for the entrez ID
#' @param regLabels \emph{Optional: } regLabels The label of the column with the regulatory labels \strong{default: "RG3_Protein"}
#' @param emptyRegLabel \emph{Optional: } emptyRegLabel The label of the empty regulatory group \strong{default: ""}
#' @param RemoveBackgroundGenes\emph{Optional: } If TRUE, genes that fall into background based on the choosen Background method for SiRCle RCM are removed from the universe. \strong{default: "TRUE"}
#' @param fileType \emph{Optional: } fileType Output file type for the figures one of: "svg", "pdf", "png" \strong{default: "pdf"}
#' @param minGSSize \emph{Optional: } minimum group size in ORA \strong{default: 10}
#' @param qvalueCutoff \emph{Optional: } q value cutoff from ORA \strong{default: 0.2}
#' @param OutputFileName Output filename
#' @param Organism \emph{Optional: } Choose either "Human" or "Mouse" \strong{default: Human"}
#'
#' @return sircleORA an instance of the sircle package
#' @export

sircleORA <- function(filename, entrezId, regLabels="RG3_Protein", emptyRegLabel="", RemoveBackgroundGenes="TRUE", fileType="pdf",
                           minGSSize=10, qvalueCutoff=0.2, pvalueCutoff=0.05, showCatagory=30, OutputFileName='', Organism="Human"){
  ## ------------ Setup and installs ----------- ##
  packages <- c("org.Hs.eg.db", "org.Mm.eg.db", "clusterProfiler", "svglite", "enrichplot")
  install.packages(setdiff(packages, rownames(installed.packages())))
  library(svglite)

  ## ------------ Create Folder ----------##
  #Create Folder, Safe file :
  SiRCleRCM_results_folder = paste(getwd(), "/SiRCleRCM",  sep="")
  if (!dir.exists(SiRCleRCM_results_folder)){dir.create(SiRCleRCM_results_folder)}#
  SiRCleRCM_ORA_results_folder = paste(SiRCleRCM_results_folder, "/SiRCleRCM_ORA_",Sys.Date(), sep="")
  if (!dir.exists(SiRCleRCM_ORA_results_folder)){dir.create(SiRCleRCM_ORA_results_folder)}  # check and create folder

  ## ------------ Run ----------- ##
  # open the data
  if(RemoveBackgroundGenes=="TRUE"){
    df <-filename
    df <- subset(df, ! df$BG_Method == "FALSE")
  } else{
    df <- read.csv(filename)
  }

  #Organism
  if(Organism=="Human"){
    Org <- org.Hs.eg.db::org.Hs.eg.db
  } else if(Organism=="Mouse"){
    Org <- org.Mm.eg.db::org.Mm.eg.db
  } else{
    stop("You have to choose an Organism.")
  }

  #define input
  allGenes <- as.character(df[[entrezId]]) #
  clusterGenes <- subset(df, ! df[[regLabels]] == emptyRegLabel)
  grps_labels <- unlist(unique(clusterGenes[regLabels]))
  for(g in grps_labels) {
    grpGenes <- subset(df, df[[regLabels]] == g)
    print(g)
    print(dim(grpGenes))
    clusterGo <- clusterProfiler::enrichGO(gene = as.character(grpGenes[[entrezId]]),
                          universe = allGenes,
                          keyType = "ENTREZID",
                          OrgDb =  Org,
                          ont = "ALL",
                          pAdjustMethod = "BH",
                          qvalueCutoff = qvalueCutoff,
                          pvalueCutoff = pvalueCutoff,
                          minGSSize = minGSSize,
                          readable = TRUE)
    # We have a cutoff of all of them, and then only visualise the ones that the user wants...

    clusterGoSummary <- data.frame(clusterGo)
    write.csv(clusterGoSummary, paste("SiRCleRCM/SiRCleRCM_ORA_", Sys.Date(), "/", OutputFileName,"_",g, ".csv", sep=""), row.names = FALSE)#Export the ORA results as .csv

    if (!(dim(clusterGoSummary)[1] == 0)){#exclude df's that have no observations
      Dotplot <- enrichplot::dotplot(clusterGo, showCategory=showCatagory) +
        ggtitle(paste("Dotplot ", g, sep=""))
      ggsave(file=paste("SiRCleRCM/SiRCleRCM_ORA_", Sys.Date(), "/", OutputFileName,"_Dotplot_", g, ".", fileType, sep=""), plot=Dotplot, width=10, height=8)
      plot(Dotplot)

      x2 <- enrichplot::pairwise_termsim(clusterGo)
      Emapplot <- enrichplot::emapplot(x2, pie_scale=1.5, layout = "nicely")+
        ggtitle(paste("Emapplot ", g, sep=""))
      ggsave(file=paste("SiRCleRCM/SiRCleRCM_ORA_", Sys.Date(), "/", OutputFileName,"_Emapplot_", g, ".", fileType, sep="" ), plot=Emapplot, width=10, height=8)
      plot(Emapplot)

      Heatplot <- enrichplot::heatplot(clusterGo,showCategory=showCatagory) +
        theme(axis.text.x =element_text(size=5), axis.text.y =element_text(size=8,face="bold"), axis.title=element_text(size=12,face="bold"))+
        ggtitle(paste("Heatplot ", g, sep=""))
      ggsave(file=paste("SiRCleRCM/SiRCleRCM_ORA_", Sys.Date(), "/", OutputFileName,"_Heatplot_", g, ".", fileType, sep="" ), plot=Heatplot, width=10, height=8)

    }
  }
}


#' sircleORA_Enrich
#'
#' Uses enricher to run ORA on each of the clusters as defined by the regulatory labels using a selfdefined pathway list
#'
#' @param filename Input file
#' @param regLabels \emph{Optional: } regLabels The label of the column with the regulatory labels \strong{default: "RegulatoryLabels"}
#' @param emptyRegLabel \emph{Optional: } emptyRegLabel The label of the empty regulatory group \strong{default: ""}
#' @param RemoveBackgroundGenes\emph{Optional: } If TRUE, genes that fall into background based on the choosen Background method for SiRCle RCM are removed from the universe. \strong{default: "TRUE"}
#' @param enricher_geneID Provide column name for the gene ID. Needs to match the the gene ID in the pathway file provided
#' @param enricher_Pathways Provide pathway file. Pathway file must include column "term" with the pathway name, column "gene" with the gene name and column "Description" with pathway description that will be depicted on the plots.
#' @param enricher_PathwayName \emph{Optional: } Name of the pathway list used \strong{default: ""}
#' @param fileType \emph{Optional: } fileType Output file type for the figures one of: "svg", "pdf", "png" \strong{default: "pdf"}
#' @param minGSSize \emph{Optional: } minimum group size in ORA \strong{default: 10}
#' @param maxGSSize \emph{Optional: } maximum group size in ORA \strong{default: 1000}
#' @param Plot_p.adj \emph{Optional: } q value cutoff from ORA that should be plotted \strong{default: 0.2}
#' @param Plot_Percentage \emph{Optional: } Percentage of genes that are detected of a pathway \strong{default: 10}
#' @param OutputFileName Output filename
#'
#' @return sircleORA an instance of the sircle package
#' @export

sircleORA_Enrich <- function(filename, regLabels="RegulatoryLabels", emptyRegLabel="", RemoveBackgroundGenes="TRUE", enricher_geneID, enricher_Pathways, enricher_PathwayName="", fileType="pdf", minGSSize=10, maxGSSize=1000 , Plot_p.adj=0.2, Plot_Percentage=10, OutputFileName=''){
  ## ------------ Setup and installs ----------- ##
  packages <- c("clusterProfiler", "enrichplot", "ggupset")
  install.packages(setdiff(packages, rownames(installed.packages())))

  ## ------------ Create Folder ----------##
  #Create Folder, Safe file :
  SiRCleRCM_results_folder = paste(getwd(), "/SiRCleRCM",  sep="")
  if (!dir.exists(SiRCleRCM_results_folder)){dir.create(SiRCleRCM_results_folder)}#
  SiRCleRCM_ORA_results_folder = paste(SiRCleRCM_results_folder, "/SiRCleRCM_ORA_",Sys.Date(), sep="")
  if (!dir.exists(SiRCleRCM_ORA_results_folder)){dir.create(SiRCleRCM_ORA_results_folder)}  # check and create folder

  ## ------------ Run ----------- ##
  # open the data
  if(RemoveBackgroundGenes=="TRUE"){
    df <- filename
    df <- subset(df, ! df$BG_Method == "FALSE")
  } else{
    df <- filename
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
  names(Pathway_Mean)[names(Pathway_Mean) == "x"] <- "Genes_in_Pathway"
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
      #Add pathway information (% of genes in pathway detected)
      clusterGoSummary <- merge(x= clusterGoSummary[,-2], y=Pathway[,-2],by.x="ID",by.y="term", all=TRUE)
      clusterGoSummary$Count[is.na(clusterGoSummary$Count)] <- 0
      clusterGoSummary$Percentage_of_Pathway_detected <-round(((clusterGoSummary$Count/clusterGoSummary$Genes_in_Pathway)*100),digits=2)
      clusterGoSummary <- clusterGoSummary[!duplicated(clusterGoSummary$ID),]
      clusterGoSummary <- clusterGoSummary[order(clusterGoSummary$p.adjust),]
      clusterGoSummary <- clusterGoSummary[,c(1,9,2:8, 10:11)]
      #Safe file
      write.csv(clusterGoSummary, paste("SiRCleRCM/SiRCleRCM_ORA_", Sys.Date(), "/", OutputFileName,"ClusterGoSummary_", enricher_PathwayName, '-', g, ".csv", sep=""), row.names = FALSE)#Export the ORA results as .csv
      #Make Selection of terms that should be displayed on the plots
      clusterGoSummary_Select <- clusterGoSummary %>%
        subset(p.adjust <= Plot_p.adj & Percentage_of_Pathway_detected >= Plot_Percentage)
      rownames(clusterGoSummary_Select)<-clusterGoSummary_Select$ID
      #Make the Plots
      if (!(dim(clusterGoSummary_Select)[1] == 0)){#exclude df's that have no observations
        clusterGo@result <- clusterGoSummary_Select[,1:9]
        #1. Dotplot:
        Dotplot <-  enrichplot::dotplot(clusterGo, showCategory=nrow(clusterGoSummary_Select)) +
          ggtitle(paste("Dotplot ", g, enricher_PathwayName, sep=" "))
        ggsave(file=paste("SiRCleRCM/SiRCleRCM_ORA_", Sys.Date(), "/", OutputFileName,"Dotplot_", g,"_" ,enricher_PathwayName, ".", fileType, sep=""), plot=Dotplot, width=10, height=8)
        plot(Dotplot)
        #2. Emapplot
        x2 <- enrichplot::pairwise_termsim(clusterGo)
        Emapplot <-  enrichplot::emapplot(x2, pie_scale=1.5, layout = "nicely", showCategory=nrow(clusterGoSummary_Select))+
          ggtitle(paste("Emapplot", g, enricher_PathwayName, sep=" "))
        ggsave(file=paste("SiRCleRCM/SiRCleRCM_ORA_", Sys.Date(), "/", OutputFileName,"Emapplot_", g,"_", enricher_PathwayName, ".", fileType, sep=""), plot=Emapplot, width=10, height=8)
        plot(Emapplot)
        #4. Upsetplot:
        UpsetPlot <-  enrichplot::upsetplot(clusterGo, showCategory=nrow(clusterGoSummary_Select))+
          ggtitle(paste("UpsetPlot", g, enricher_PathwayName, sep=" "))
        ggsave(file=paste("SiRCleRCM/SiRCleRCM_ORA_", Sys.Date(), "/", OutputFileName,"UpsetPlot_", g,"_", enricher_PathwayName, ".", fileType, sep=""), plot=UpsetPlot, width=10, height=8)
        plot(UpsetPlot)
      }
    }
  }
}


