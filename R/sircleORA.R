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

#' sircleORAHuman
#'
#' Uses encrichGo to run ORA on each of the clusters as defined by the regulatory labels (ON HUMAN)
#'
#' @param filename Path to the input file
#' @param entrezId Column name for the entrez ID
#' @param regLabels \emph{Optional: } regLabels The label of the column with the regulatory labels \strong{default: "RegulatoryLabels"}
#' @param emptyRegLabel \emph{Optional: } emptyRegLabel The label of the empty regulatory group \strong{default: ""}
#' @param fileType \emph{Optional: } fileType Output file type for the figures one of: "svg", "pdf", "png" \strong{default: "pdf"}
#' @param minGSSize \emph{Optional: } minimum group size in ORA \strong{default: 10}
#' @param qvalueCutoff \emph{Optional: } q value cutoff from ORA \strong{default: 0.2}
#'
#' @return
#' @export

sircleORAHuman <- function(filename, entrezId, regLabels="RegulatoryLabels", emptyRegLabel="", fileType="pdf", minGSSize=10, qvalueCutoff=0.2) {
  ## ------------ Setup and installs ----------- ##
  packages <- c("org.Hs.eg.db", "clusterProfiler", "svglite", "enrichplot")
  install.packages(setdiff(packages, rownames(installed.packages())))
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(svglite)
  library(enrichplot)
  ## ------------ Run ----------- ##
  # open the data
  df <- read.csv(filename)
  allGenes <- as.character(df[[entrezId]]) #
  clusterGenes <- subset(df, ! df[[regLabels]] == emptyRegLabel)
  grps_labels <- unlist(unique(clusterGenes[regLabels]))
  for(g in grps_labels) {
    grpGenes <- subset(df, df[[regLabels]] == g)
    print(g)
    print(dim(grpGenes))
    clusterGo <- enrichGO(gene = as.character(grpGenes[[entrezId]]),
                          universe = allGenes,
                          keyType = "ENTREZID",
                          OrgDb = org.Hs.eg.db,
                          ont = "ALL",
                          pAdjustMethod = "BH",
                          qvalueCutoff = qvalueCutoff,
                          minGSSize = minGSSize,
                          readable = TRUE)
    clusterGoSummary <- data.frame(clusterGo)
    write_csv(clusterGoSummary, paste('ClusterGoSummary_', g, '.csv', sep=""))#Export the ORA results as .csv

    if (!(dim(clusterGoSummary)[1] == 0)) {#exclude df's that have no observations
      Dotplot <- dotplot(clusterGo, showCategory=30) +
        ggtitle(paste("Dotplot ", g, sep=""))
      ggsave(file=paste("SiRCle-ORA_Dotplot_Human_", g, ".", fileType, sep=""), plot=Dotplot, width=10, height=8)
      x2 <- pairwise_termsim(clusterGo)

      Emapplot <- emapplot(x2, pie_scale=1.5, layout = "nicely")+
        ggtitle(paste("Emapplot ", g, sep=""))
      ggsave(file=paste("SiRCle-ORA_Emapplot_Human_", g, ".", fileType, sep="" ), plot=Emapplot,width=10, height=8)

      Heatplot <- heatplot(clusterGo,showCategory=30) +
        theme(axis.text.x =element_text(size=5), axis.text.y =element_text(size=8,face="bold"), axis.title=element_text(size=12,face="bold"))+
        ggtitle(paste("Heatplot ", g, sep=""))
      ggsave(file=paste("SiRCle-ORA_Heatplot_Human_", g, ".", fileType, sep="" ), plot=Heatplot,width=10, height=8)

      }
  }
}

#' sircleORAMouse
#'
#' Uses encrichGo to run ORA on each of the clusters as defined by the regulatory labels (on MOUSE)
#'
#' @param rcm Instance of RCM after running sircleRCM
#' @param filename Path to the input file
#' @param envPath \emph{Optional: } regLabels The label of the column with the regulatory labels \strong{default: "RegulatoryLabels"}
#' @param envPath \emph{Optional: } fileType Output file type for the figures one of: "svg", "pdf", "png" \strong{default: "pdf"}
#' @return
#' @export
sircleORAMouse<- function(filename, regLabels="RegulatoryLabels", fileType="pdf") {
  ## ------------ Setup and installs ----------- ##
  packages <- c("org.Mm.eg.db", "clusterProfiler", "svglite", "enrichplot")

  install.packages(setdiff(packages, rownames(installed.packages())))
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(svglite)

  ## ------------ Run ----------- ##

  # open the data
  df <- read.csv(filename)
  allGenes <- as.character(df[[entrezId]])# I had to add as.character
  clusterGenes <- subset(df, ! df[[regLabels]] == emptyRegLabel)
  #oraDir <- "ora_figs" #This did not work for me, and I just had to make the folder myself before running the code
  grps_labels <- unlist(unique(clusterGenes[regLabels]))
  for(g in grps_labels) {
    grpGenes <- unlist(rcm$get_genes_in_reg_grp(g, entrezId))#I had to add unlist
    clusterGo <- enrichGO(gene = grpGenes,
                          universe = allGenes,
                          keyType = "ENTREZID",
                          OrgDb = org.Mm.eg.db,
                          ont = "ALL",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.1,
                          readable = TRUE)
    clusterGoSummary <- data.frame(clusterGo)
    write_csv(clusterGoSummary, paste('ClusterGoSummary_', g, '.csv', sep=""))#Export the ORA results as .csv
    if (!(dim(clusterGoSummary)[1] == 0)) {#exclude df's that have no observations

      Dotplot <- dotplot(clusterGo, showCategory=30) +
        ggtitle(paste("Dotplot ", g, sep=""))
      ggsave(file=paste("SiRCle-ORA_Dotplot_Mouse_", g, ".", fileType, sep=""), plot=Dotplot, width=10, height=8)

      Emapplot <- emapplot(clusterGo, pie_scale=1.5, showCategory=30, layout = "nicely")+
        ggtitle(paste("Emapplot ", g, sep=""))
      ggsave(file=paste("SiRCle-ORA_Emapplot_Mouse_", g, ".", fileType, sep="" ), plot=Emapplot,width=10, height=8)

      Heatplot <- heatplot(clusterGo,showCategory=30) +
        theme(axis.text.x =element_text(size=5), axis.text.y =element_text(size=8,face="bold"),  axis.title=element_text(size=12,face="bold"))+
        ggtitle(paste("Heatplot ", g, sep=""))
      ggsave(file=paste("SiRCle-ORA_Heatplot_Mouse_", g, ".", fileType, sep="" ),
             plot=Heatplot,width=10, height=8)
    }
  }
}
