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

#' sircleTFVis
#'
#' Creates visualisations for the TF analysis.
#'
#' @param filename Column name with the gene identifier
#' @param colorCol \emph{Optional:} Column name to use for colouring \strong{Default: "q-value"}
#' @param fillNAVal \emph{Optional:} Value to fill None values with \strong{Default: 0.0}
#' @param plotBar \emph{Optional:} Whether or not to plot the bar charts \strong{Default: T}
#' @param plotEman \emph{Optional:} Whether or not to plot the connection plot \strong{Default: T}
#' @param titleStr \emph{Optional:} String for the title \strong{Default: ""}
#' @return scimo an instance of the scimo package
#' @export

sircleTFVis <- function(filename, colorCol="q-value", fillNAVal=NA, plotBar=T, plotEman=T, titleStr="") {
  ## ------------ Setup and installs  ----------- ##
  packages <- c("ggraph", "igraph", "tidygraph", "purrr", "rlang", "dplyr")

  install.packages(setdiff(packages, rownames(installed.packages())))
  library(ggraph)
  library(igraph)
  library(tidygraph)
  library(purrr)
  library(rlang)
  #Understand the data and amount of TFs in each regulatory cluster:
  TF_reg <- read.csv(filename)

  # Add the "size" to each TF, in order to be able to sum them up -->  This will be 1 for each
  TF_reg[,"size"]  <- as.numeric(1)

  #Get the sum of the TFs detected in each cluster:
  TF_reg_inClust_All <- aggregate(TF_reg[, "size"], by=list(TF_reg$cluster), FUN=sum)#Built the sum

  #There are TFs that are detected within multiple clusters. Now I will exclude all TFs that are duplicated, so that we get an overview of the TFs that are unique to each cluster:
  doublons <-  as.character(TF_reg[duplicated(TF_reg$motif), "motif"])
  doublons <- as.data.frame(doublons)
  matches <- match(TF_reg$motif, doublons$doublons)
  matches <- as.data.frame(matches)
  matches <- cbind(TF_reg, matches)
  TF_reg_Unique <- matches[is.na(matches$matches),]
  TF_reg_inClust_Unique <- aggregate(TF_reg_Unique[,"size"], by=list(TF_reg_Unique$cluster), FUN=sum)#Built the sum

  # Make a combined data frame:
  TF <- merge(TF_reg_inClust_All, TF_reg_inClust_Unique,by="Group.1", all=TRUE)
  names(TF)[names(TF) == "x.x"] <- "TF_in_Cluster"#Rename the rowname
  names(TF)[names(TF) == "x.y"] <- "TF_in_Cluster_Unique"#Rename the rowname
  names(TF)[names(TF) == "Group.1"] <- "cluster"#Rename the rowname

  ##################################
  # Visualise the unique TFs making a bar graph/dot plot:
  TF_reg_Unique_Percentage <- TF_reg_Unique[,c(1:2, 13:14)]%>%
    separate(col="motif", into=c("motif", "id"), sep="_")
  names(TF_reg_Unique_Percentage)[names(TF_reg_Unique_Percentage) == "X..coverage"] <- "Percentage_in_Cluster"#Rename the rowname
  # Add in padj
  TF_reg_Unique_Percentage$color <- TF_reg_Unique[[colorCol]]
  TF_reg_Unique_Percentage$cluster <- TF_reg_Unique$cluster

  #Subset the regulatory clusters and order by percentage (Lowest to Highest):
  MDS <- subset(TF_reg_Unique_Percentage, cluster == "MDS")%>%
    arrange(Percentage_in_Cluster)
  TPDE_TMDS <- subset(TF_reg_Unique_Percentage, cluster == "TPDE_TMDS")%>%
    arrange(Percentage_in_Cluster)
  TPDE <- subset(TF_reg_Unique_Percentage, cluster == "TPDE")%>%
    arrange(Percentage_in_Cluster)
  TMDE <- subset(TF_reg_Unique_Percentage, cluster == "TMDE") %>%
    arrange(Percentage_in_Cluster)
  TPDS_TMDE <- subset(TF_reg_Unique_Percentage, cluster == "TPDS_TMDE")%>%
    arrange(Percentage_in_Cluster)
  TPDS <- subset(TF_reg_Unique_Percentage, cluster == "TPDS")%>%
    arrange(Percentage_in_Cluster)
  TMDS <- subset(TF_reg_Unique_Percentage, cluster == "TMDS")%>%
    arrange(Percentage_in_Cluster)
  MDE <- subset(TF_reg_Unique_Percentage, cluster == "MDE")%>%
    arrange(Percentage_in_Cluster)
  MDE_TMDS <- subset(TF_reg_Unique_Percentage, cluster == "MDE_TMDS")%>%
    arrange(Percentage_in_Cluster)
  MDS_TMDE <- subset(TF_reg_Unique_Percentage, cluster == "MDS_TMDE")%>%
    arrange(Percentage_in_Cluster)

  sircleTFEman(MDS, colorCol, fillNAVal, paste("MDS", titleStr, sep=""),  plotBar=plotBar, plotEman=plotEman)
  sircleTFEman(TPDE_TMDS, colorCol, fillNAVal, paste("TPDE_TMDS", titleStr, sep=""),  plotBar=plotBar, plotEman=plotEman)
  sircleTFEman(TPDE, colorCol, fillNAVal, paste("TPDE", titleStr, sep=""),  plotBar=plotBar, plotEman=plotEman)
  sircleTFEman(TMDE, colorCol, fillNAVal, paste("TMDE", titleStr, sep=""),  plotBar=plotBar, plotEman=plotEman)
  sircleTFEman(TPDS_TMDE, colorCol, fillNAVal, paste("TPDS_TMDE", titleStr, sep=""),  plotBar=plotBar, plotEman=plotEman)
  sircleTFEman(TPDS, colorCol, fillNAVal, paste("TPDS", titleStr, sep=""),  plotBar=plotBar, plotEman=plotEman)
  sircleTFEman(TMDS, colorCol, fillNAVal, paste("TMDS", titleStr, sep=""),  plotBar=plotBar, plotEman=plotEman)
  sircleTFEman(MDE, colorCol, fillNAVal, paste("MDE", titleStr, sep=""),  plotBar=plotBar, plotEman=plotEman)
  sircleTFEman(MDE_TMDS, colorCol, fillNAVal, paste("MDE_TMDS", titleStr, sep=""),  plotBar=plotBar, plotEman=plotEman)
  sircleTFEman(MDS_TMDE, colorCol, fillNAVal, paste("MDS_TMDE", titleStr, sep=""), plotBar=plotBar, plotEman=plotEman)
  sircleTFEman(TF_reg_Unique_Percentage, "cluster", fillNAVal=fillNAVal, titleStr=paste("All Clusters", titleStr, sep=""), continuous=F, plotBar=plotBar, plotEman=plotEman)

}


#' sircleTFEman
#'
#' Creates the eman plot for the TF vis
#'
#' @param TfDf Column name with the gene identifier
#' @param colorLabel Column with the value that we will filter on (usually the dorothea column)
#' @param fillNAVal \emph{Optional:}  Fill value for None columns \strong{Default: NA}
#' @param continuous \emph{Optional:} Whether or not the data are continous \strong{Default: T}
#' @param titleStr \emph{Optional:} String for the title \strong{Default: ""}
#' @param plotBar \emph{Optional:} Whether or not to plot the bar charts \strong{Default: T}
#' @param plotEman \emph{Optional:} Whether or not to plot the connection plot \strong{Default: T}
#' @return scimo an instance of the scimo package
#' @export
sircleTFEman <- function(TfDf, colorLabel, fillNAVal=NA, titleStr="", continuous=T, plotBar=T, plotEman=T) {

  TfDf$color[TfDf$color =="None"]<- fillNAVal  # Fill any None values
  TfDf$color <- as.numeric(TfDf$color) # Set as numeric

  if (plotBar) {
    # First plot the bar chart
    pdf(file=paste("TF-Bar_", titleStr,".pdf", sep=""))
    barplot(TfDf$Percentage_in_Cluster, main=titleStr, horiz=TRUE,
            names.arg=TfDf$motif, xlab="Amount of genes in cluster driven by the TF [%]", col="#69b3a2", las=1)
    dev.off()
  }

  if (plotEman) {
    # Now build our eman plot
    fromTf <- c()
    toTf <- c()
    motifsDone <- c()
    sharedGeneCount <- c()
    # Basically we get the shared gene count since % shared might not nake sense (different TFs will
    # have different amounts of genes per TF meaning it would require two values per line which is a bit much)
    for (motif1 in TfDf$motif) {
      motif1Df <- subset(TfDf, motif == motif1)
      for (motif2 in TfDf$motif) {
        if (!(motif2 %in% motifsDone)) {
          motif2Df <- subset(TfDf, motif == motif2)
          if (motif1 != motif2) {
            # Check the count of shared genes
            genesMotif1 <-  unlist(strsplit(motif1Df$genes, ","))
            genesMotif2 <- unlist(strsplit(motif2Df$genes, ","))
            countOverlap <- length(intersect(genesMotif1, genesMotif2))
            if (countOverlap > 0) {
              sharedGeneCount <- rbind(sharedGeneCount, countOverlap)
              fromTf <- rbind(fromTf, motif1)
              toTf <- rbind(toTf, motif2)
            }
          }
        }
      }
      # Keep track of which motifs we've done all pairs for i.e. do we don't get duplicates
      motifsDone <- rbind(motifsDone, motif1)
    }

    # Check if they actually have enough for a graph
    fromTFs <- unlist(fromTf)
    if (length(fromTFs) < 2) {
      print("You didn't have enough TFs (less than 2) so we are returning without plotting!")
      return()
    }
    if (length(unlist(toTf)) < 2) {
      print("You didn't have enough TFs (less than 2) so we are returning without plotting!")
      return()
    }
    TF <- data.frame(from = unlist(fromTf),
                     to = unlist(toTf))

    vertexMetadata <- data.frame(name = unlist(TfDf$motif),
                                 color = unlist(TfDf$color),
                                 cluster = unlist(TfDf$cluster),
                                 size = unlist(TfDf$Percentage_in_Cluster))
    # Create from the TF DF and then vertexMetdata

    TF_graph <- graph_from_data_frame(TF, directed = TRUE, vertices=vertexMetadata)
    # Add in some extra columns
    TF_graph$overlap <- unlist(sharedGeneCount)

    # Make the graph:
    if (continuous) {
      tfplt <- ggraph(TF_graph,
                      layout ="stress")+ #layout are the node-positions and (possible) attributes
        geom_edge_link(alpha=0.8,
                       width=3,#thickness of the edge
                       color ="grey",
                       aes(label=TF_graph$overlap),#overlap between TF driven genes
                       angle_calc = 'along',#position of the label on the edge
                       label_dodge = unit(2.5, 'mm'),#position below/above the line
                       label_colour = "black")+
        geom_node_point(alpha=0.9,
                        shape=21,
                        color="black",
                        aes_(size=vertex_attr(TF_graph, "size"), # This is how we can access the vertex attributes
                             fill=vertex_attr(TF_graph, "color"))) + #colour according to adjusted p-value
        scale_size(range=c(6, 20),
                   name="% of cluster genes regulated by TF")+
        scale_fill_continuous(low="red",
                              high="blue",
                              name = colorLabel,
                              guide=guide_colorbar(reverse=TRUE)) +
        geom_node_text(aes(label = name),
                       fontface = "bold",
                       nudge_x = 0,
                       nudge_y = 0)+
        theme(panel.background =element_rect(fill="white")) +
        theme(legend.position = "right",
              legend.key = element_rect(fill = "white"),
              legend.key.size = unit(0.5, "cm"),
              legend.key.width = unit(0.5,"cm"))

      ggsave(file=paste("TF-plot_", titleStr,".pdf", sep=""), plot=tfplt)
    } else {
      tfplt <- ggraph(TF_graph,
                      layout ="stress")+ #layout are the node-positions and (possible) attributes
        geom_edge_link(alpha=0.8,
                       width=3,#thickness of the edge
                       color ="grey",
                       aes(label=TF_graph$overlap),#overlap between TF driven genes
                       angle_calc = 'along',#position of the label on the edge
                       label_dodge = unit(2.5, 'mm'),#position below/above the line
                       label_colour = "black")+
        geom_node_point(alpha=0.9,
                        shape=21,
                        color="black",
                        aes_(size=vertex_attr(TF_graph, "size"), # This is how we can access the vertex attributes
                             fill=vertex_attr(TF_graph, "cluster"))) + #colour according to adjusted p-value
        scale_size(range=c(6, 20),
                   name="% of cluster genes regulated by TF")+
        scale_fill_discrete(name = "cluster") +
        geom_node_text(aes(label = name),
                       fontface = "bold",
                       nudge_x = 0,
                       nudge_y = 0)+
        theme(panel.background =element_rect(fill="white")) +
        theme(legend.position = "right",
              legend.key = element_rect(fill = "white"),
              legend.key.size = unit(0.5, "cm"),
              legend.key.width = unit(0.5,"cm"))

      ggsave(file=paste("TF-plot_", titleStr,".pdf", sep="" ), plot=tfplt)
    }
  }
}