## ---------------------------
##
## Script name: sircleFET
##
## Purpose of script: Performs Fishers Exact Test (FET) on sircleRCM output.
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