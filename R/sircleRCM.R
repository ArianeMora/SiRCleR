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

#' sircleRCM_MRP
#'
#' Computes the regulatory clustering model (RCM) using the SiRCle regulatory rules for DNA-methylation, RNAseq and Proteomics data layers (MRP).
#'
#' @param rnaFile Filename for your RNAseq data (results from DeSeq2 and also your normalised expression counts)
#' @param methFile  Filename for your DNA methylation data (results from differential methylation analysis)
#' @param proteinFile Filename/path of you Protein data (results from DeSeq2 and also your normalised expression counts)
#' @param geneID Column name of geneId this MUST BE THE SAME in each of your protein, RNAseq and DNAmethylation files (we join on this)
#' @param rnaValueCol Column name of RNA value in rnaFile 
#' @param rnaPadjCol Column name of RNA p adjusted value in rnaFile 
#' @param methValueCol Column name of Methylation difference value in methFile 
#' @param methPadjCol Column name of Methylation p adjusted value in methFile 
#' @param proteinValueCol  Column name of protein log fold change in proteinFile
#' @param proteinPadjCol Column name of protein p adjusted value in proteinFile
#' @param rnaPadjCutoff  \emph{Optional: }Padjusted cutoff for RNAseq data \strong{Default=0.05}
#' @param rnaLogFCCutoff \emph{Optional: } LogFoldchange cutoff for RNAseq data \strong{Default=0.5}
#' @param proteinPadjCutoff \emph{Optional: } Padjusted cutoff for Protein data \strong{Default=0.05}
#' @param proteinValueCutoff \emph{Optional: } LogFoldchange cutoff for Protein data \strong{Default=0.3}
#' @param methPadjCutoff \emph{Optional: } Padjusted cutoff for DNA methylation data \strong{Default=0.05}
#' @param methDiffCutoff \emph{Optional: } DNA Methylation difference cutoff for DNA methylation \strong{Default=10}
#' @param backgroundMethod \emph{Optional: } Background method (NEED Description for each one here) \strong{Default="P|(M&R)"}
#' @param outputFileName \emph{Optional: } Output filename \strong{Default=SiRCle_RCM.csv}
#' @return rcm an instance of the rcm package
#' @export
#'
sircleRCM_MRP <- function(methFile, rnaFile, protFile, geneID, rnaValueCol="Log2FC", rnaPadjCol="padj", methValueCol="Diff", methPadjCol="padj", proteinValueCol="Log2FC", proteinPadjCol="padj", rna_padj_cutoff= 0.05, prot_padj_cutoff = 0.05, meth_padj_cutoff = 0.05,rna_FC_cutoff= 1, prot_FC_cutoff = 0.5, meth_Diff_cutoff = 0.1, backgroundMethod, OutputFileName = "Sircle_RCM.csv"){
  #Import the data:
  proteinDF <- as.data.frame(protFile)%>%
    rename("geneID"=paste(geneID),
           "ValueCol"=paste(proteinValueCol),
           "PadjCol"=paste(proteinPadjCol))
  rnaDF<- as.data.frame(rnaFile)%>%
    rename("geneID"=paste(geneID),
           "ValueCol"=paste(rnaValueCol),
           "PadjCol"=paste(rnaPadjCol))
  methDF<- as.data.frame(methFile)%>%
    rename("geneID"=paste(geneID),
           "ValueCol"=paste(methValueCol),
           "PadjCol"=paste(methPadjCol))
  
  #First check for duplicates in "geneID" and drop if there are any
  if(length(proteinDF[duplicated(proteinDF$geneID), "geneID"]) > 0){
    doublons <- as.character(proteinDF[duplicated(proteinDF$geneID), "geneID"])#number of duplications
    proteinDF <-proteinDF[!duplicated(proteinDF$geneID),]#remove duplications
    warning("Proteomics dataset contained duplicates based on geneID! Dropping duplicate IDs and kept only the first entry. You had ", length(doublons), " duplicates.")
    warning("Note that you should do this before running SiRCle.")
  }
  if(length(rnaDF[duplicated(rnaDF$geneID), "geneID"]) > 0){
    doublons <- as.character(rnaDF[duplicated(rnaDF$geneID), "geneID"])#number of duplications
    rnaDF <-rnaDF[!duplicated(rnaDF$geneID),]#remove duplications
    warning("RNAseq dataset contained duplicates based on geneID! Dropping duplicate IDs and kept only the first entry. You had ", length(doublons), " duplicates.")
    warning("Note that you should do this before running SiRCle.")
  }
  if(length(methDF[duplicated(methDF$geneID), "geneID"]) > 0){
    doublons <- as.character(methDF[duplicated(methDF$geneID), "geneID"])#number of duplications
    methDF <-methDF[!duplicated(methDF$geneID),]#remove duplications
    warning("DNA-methylation dataset contained duplicates based on geneID! Dropping duplicate IDs and kept only the first entry. You had ", length(doublons), " duplicates.")
    warning("Note that you should do this before running SiRCle.")
  }
  
  #Tag genes that are detected in each data layer
  proteinDF$Detected <- "TRUE"
  rnaDF$Detected <- "TRUE"
  methDF$Detected <- "TRUE"
  
  #Assign to Group based on individual Cutoff ("UP", "DOWN", "No Change")
  proteinDF <- proteinDF%>%
    mutate(Cutoff = case_when(proteinDF$PadjCol < prot_padj_cutoff & proteinDF$ValueCol > prot_FC_cutoff ~ 'UP',
                              proteinDF$PadjCol < prot_padj_cutoff & proteinDF$ValueCol < -prot_FC_cutoff ~ 'DOWN',
                              TRUE ~ 'No Change')) %>%
    mutate(Cutoff_Specific = case_when(Cutoff == "UP" ~ 'UP',
                                       Cutoff == "DOWN" ~ 'DOWN',
                                       Cutoff == "No Change" & proteinDF$PadjCol < prot_padj_cutoff & proteinDF$ValueCol > 0 ~ 'Significant Positive',
                                       Cutoff == "No Change" & proteinDF$PadjCol < prot_padj_cutoff & proteinDF$ValueCol < 0 ~ 'Significant Negative',
                                       Cutoff == "No Change" & proteinDF$PadjCol > prot_padj_cutoff ~ 'Not Significant',
                                       TRUE ~ 'FALSE'))
  
  rnaDF <- rnaDF%>%
    mutate(Cutoff = case_when(rnaDF$PadjCol < rna_padj_cutoff & rnaDF$ValueCol > rna_FC_cutoff ~ 'UP',
                              rnaDF$PadjCol < rna_padj_cutoff & rnaDF$ValueCol < -rna_FC_cutoff ~ 'DOWN',
                              TRUE ~ 'No Change'))
  
  methDF <- methDF%>%
    mutate(Cutoff = case_when(methDF$PadjCol < meth_padj_cutoff & methDF$ValueCol > meth_Diff_cutoff ~ 'Hypermethylation',
                              methDF$PadjCol < meth_padj_cutoff & methDF$ValueCol < -meth_Diff_cutoff ~ 'Hypomethylation',
                              TRUE ~ 'No Change'))
  
  #Merge the dataframes together: Merge the supplied RNAseq, DNA methylation and proteomics dataframes together. 
  ##Add prefix to column names to distinguish the different data types after merge
  colnames(proteinDF) <- paste0("proteinDF_", colnames(proteinDF))
  proteinDF <- proteinDF%>%
    rename("geneID" = "proteinDF_geneID")
  
  colnames(rnaDF) <- paste0("rnaDF_", colnames(rnaDF))
  rnaDF <- rnaDF%>%
    rename("geneID"="rnaDF_geneID")
  
  colnames(methDF) <- paste0("methDF_", colnames(methDF))
  methDF <- methDF%>%
    rename("geneID"="methDF_geneID")
  
  ##Merge
  MergeDF <- merge(methDF, rnaDF, by="geneID", all=TRUE)
  MergeDF <- merge(MergeDF, proteinDF, by="geneID", all=TRUE)
  
  ##Mark the undetected genes in each data layer
  MergeDF<-MergeDF %>% 
    mutate_at(c("proteinDF_Detected","rnaDF_Detected", "methDF_Detected"), ~replace_na(.,"FALSE"))%>% 
    mutate_at(c("proteinDF_Cutoff","rnaDF_Cutoff", "methDF_Cutoff"), ~replace_na(.,"No Change"))%>%
    mutate_at(c("proteinDF_Cutoff_Specific"), ~replace_na(.,"Not Detected"))
  
  #Apply Background filter (label genes that will be removed based on choosen background)
  if(backgroundMethod == "P|(M&R)"){# P|(M&R) = Protein AND (DNA methylation OR RNA)
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(methDF_Detected=="TRUE" & rnaDF_Detected=="TRUE" & proteinDF_Detected=="TRUE" ~ 'TRUE', # RNA & Methylation & Protein
                                   methDF_Detected=="TRUE" & rnaDF_Detected=="TRUE" & proteinDF_Detected=="FALSE" ~ 'TRUE', # RNA & Methylation ~Protein
                                   methDF_Detected=="TRUE" & rnaDF_Detected=="FALSE" & proteinDF_Detected=="TRUE" ~ 'TRUE', # ~ methylation RNA Protein
                                   methDF_Detected=="FALSE" & rnaDF_Detected=="TRUE" & proteinDF_Detected=="TRUE" ~ 'TRUE', # Methylation ~ RNA & protein
                                   methDF_Detected=="FALSE" & rnaDF_Detected=="FALSE" & proteinDF_Detected=="TRUE" ~ 'TRUE', # Just protein
                                   TRUE ~ 'FALSE'))
  }
  else if(backgroundMethod == "P|M|R"){ # Protein OR methylation OR RNA
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(methDF_Detected=="FALSE" & rnaDF_Detected=="FALSE" & proteinDF_Detected=="FALSE" ~ 'FALSE', # i.e. only one we don't want is Not detected in all
                                   TRUE ~ 'TRUE'))
  }
  else if(backgroundMethod == "P|R"){# Protein OR RNA
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(methDF_Detected=="FALSE" & rnaDF_Detected=="FALSE" & proteinDF_Detected=="FALSE" ~ 'FALSE', # i.e. only one we don't want is Not detected in all
                                   methDF_Detected=="TRUE" & rnaDF_Detected=="FALSE" & proteinDF_Detected=="FALSE" ~ 'FALSE',  # i.e. only one we don't want is NS in protein and RNA
                                   TRUE ~ 'TRUE'))
  }
  else if(backgroundMethod == "P&R"){# Protein AND RNA
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(methDF_Detected=="TRUE" & rnaDF_Detected=="TRUE" & proteinDF_Detected=="TRUE" ~ 'TRUE', # RNA & Methylation & Protein
                                   methDF_Detected=="FALSE" & rnaDF_Detected=="TRUE" & proteinDF_Detected=="TRUE" ~ 'TRUE', # Methylation ~ RNA & protein
                                   TRUE ~ 'FALSE'))
  }
  else if(backgroundMethod == "P&M&R"){# Protein AND Methylation AND RNA
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(methDF_Detected=="TRUE" & rnaDF_Detected=="TRUE" & proteinDF_Detected=="TRUE" ~ 'TRUE', # RNA & Methylation & Protein
                                   TRUE ~ 'FALSE'))
  }
  else if(backgroundMethod == "(P&M)|(P&R)|(M&R)"){ # At least two are significant
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(methDF_Detected=="TRUE" & rnaDF_Detected=="TRUE" & proteinDF_Detected=="TRUE" ~ 'TRUE', # RNA & Methylation & Protein
                                   methDF_Detected=="TRUE" & rnaDF_Detected=="FALSE" & proteinDF_Detected=="TRUE" ~ 'TRUE', # P&M
                                   methDF_Detected=="FALSE" & rnaDF_Detected=="TRUE" & proteinDF_Detected=="TRUE" ~ 'TRUE', # P&R
                                   methDF_Detected=="TRUE" & rnaDF_Detected=="TRUE" & proteinDF_Detected=="FALSE" ~ 'TRUE', # M&R
                                   TRUE ~ 'FALSE'))
  }
  else if(backgroundMethod == "(P&M)|(P&R)"){ # Protein and one other
    MergeDF <- MergeDF%>%
      mutate(BG_Method = case_when(methDF_Detected=="TRUE" & rnaDF_Detected=="TRUE" & proteinDF_Detected=="TRUE" ~ 'TRUE', # RNA & Methylation & Protein
                                   methDF_Detected=="TRUE" & rnaDF_Detected=="FALSE" & proteinDF_Detected=="TRUE" ~ 'TRUE', # P&M
                                   methDF_Detected=="FALSE" & rnaDF_Detected=="TRUE" & proteinDF_Detected=="TRUE" ~ 'TRUE', # P&R
                                   TRUE ~ 'FALSE'))
  }
  else if(backgroundMethod == "*"){ # Use all genes as the background
    MergeDF$BG_Method <- "TRUE"
  }
  else{
    stop("Please use one of the following backgroundMethods: P|(M&R), P|M|R, P|R, P&R, P&M&R, (P&M)|(P&R)|(M&R), (P&M)|(P&R), *")#error message
  }
  
  #Assign SiRCle cluster names to the genes
  MergeDF <- MergeDF%>%
    mutate(RG1_All = case_when(BG_Method =="FALSE"~ 'Background = FALSE',
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="DOWN" ~ 'Hypermethylation + RNA DOWN + Protein DOWN',#State 1
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="DOWN" ~ 'Hypermethylation + RNA No change + Protein DOWN',#State 2
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="DOWN" ~ 'Hypermethylation + RNA UP + Protein DOWN',#State 3
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'Hypermethylation + RNA DOWN + Protein Undetected',#State 7
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'Hypermethylation + RNA No change + Protein Undetected',#State 8
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'Hypermethylation + RNA UP + Protein Undetected',#State 9
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'Hypermethylation + RNA DOWN + Protein not-significant',#State 4
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'Hypermethylation + RNA No change + Protein not-significant',#State 5
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'Hypermethylation + RNA UP + Protein not-significant',#State 6
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'Hypermethylation + RNA DOWN + Protein significant-negative',#State 10
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'Hypermethylation + RNA No change + Protein significant-negative',#State 11
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'Hypermethylation + RNA UP + Protein significant-negative',#State 12
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'Hypermethylation + RNA DOWN + Protein significant-positive',#State 10
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'Hypermethylation + RNA No change + Protein significant-positive',#State 11
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'Hypermethylation + RNA UP + Protein significant-positive',#State 12
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="UP" ~ 'Hypermethylation + RNA DOWN + Protein UP',#State 10
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="UP" ~ 'Hypermethylation + RNA No change + Protein UP',#State 11
                               methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="UP" ~ 'Hypermethylation + RNA UP + Protein UP',#State 12
                               
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="DOWN" ~ 'Hypomethylation + RNA DOWN + Protein DOWN',#State 13
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="DOWN" ~ 'Hypomethylation + RNA No change + Protein DOWN',#State 14
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="DOWN" ~ 'Hypomethylation + RNA UP + Protein DOWN',#State 15
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'Hypomethylation + RNA DOWN + Protein Undetected',#State 22
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'Hypomethylation + RNA No change + Protein Undetected',#State 23
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'Hypomethylation + RNA UP + Protein Undetected',#State 24
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'Hypomethylation + RNA DOWN + Protein not-significant',#State 16
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'Hypomethylation + RNA No change + Protein not-significant',#State 17
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'Hypomethylation + RNA UP + Protein not-significant',#State 18
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'Hypomethylation + RNA DOWN + Protein significant-negative',#State 13
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'Hypomethylation + RNA No change + Protein significant-negative',#State 14
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'Hypomethylation + RNA UP + Protein significant-negative',#State 15
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'Hypomethylation + RNA DOWN + Protein significant-positive',#State 13
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'Hypomethylation + RNA No change + Protein significant-positive',#State 14
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'Hypomethylation + RNA UP + Protein significant-positive',#State 15
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="UP" ~ 'Hypomethylation + RNA DOWN + Protein UP',#State 19
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="UP" ~ 'Hypomethylation + RNA No change + Protein UP',#State 20
                               methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="UP" ~ 'Hypomethylation + RNA UP + Protein UP',#State 21
                               
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="DOWN" ~ 'Methylation No Change + RNA DOWN + Protein DOWN',#State 25
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="DOWN" ~ 'Methylation No Change + RNA No change + Protein DOWN',#State 26
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="DOWN" ~ 'Methylation No Change + RNA UP + Protein DOWN',#State 27
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'Methylation No Change + RNA DOWN + Protein Undetected',#State 31
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'Methylation No Change + RNA No change + Protein Undetected',#State 32
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'Methylation No Change + RNA UP + Protein Undetected',#State 33
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'Methylation No Change + RNA DOWN + Protein not-significant',#State 28
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'Methylation No Change + RNA No change + Protein not-significant',#State 29
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'Methylation No Change + RNA UP + Protein not-significant',#State 30
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'Methylation No Change + RNA DOWN + Protein significant-negative',#State 25
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'Methylation No Change + RNA No change + Protein significant-negative',#State 26
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'Methylation No Change + RNA UP + Protein significant-negative',#State 27
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'Methylation No Change + RNA DOWN + Protein significant-positive',#State 34
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'Methylation No Change + RNA No change + Protein significant-positive',#State 35
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'Methylation No Change + RNA UP + Protein significant-positive',#State 36
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="UP" ~ 'Methylation No Change + RNA DOWN + Protein UP',#State 34
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="UP" ~ 'Methylation No Change + RNA No change + Protein UP',#State 35
                               methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="UP" ~ 'Methylation No Change + RNA UP + Protein UP',#State 36
                               TRUE ~ 'NA'))%>%
    mutate(RG2_Changes = case_when(BG_Method =="FALSE"~ 'Background = FALSE',
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="DOWN" ~ 'MDS',#State 1
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TMDS',#State 2
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TPDE_TMDS',#State 3
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'MDS_TMDE',#State 7
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'None',#State 8
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'TPDE_TMDS',#State 9
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'MDS_TMDE',#State 4
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'None',#State 5
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'TPDE_TMDS',#State 6
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'MDS_TMDE',#State 10
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 11
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'TPDE_TMDS',#State 12
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'MDS_TMDE',#State 10
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 11
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'TPDE_TMDS',#State 12
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="UP" ~ 'MDS_TMDE',#State 10
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="UP" ~ 'TMDE',#State 11
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="UP" ~ 'TPDE',#State 12
                                   
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TPDS',#State 13
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TMDS',#State 14
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="DOWN" ~ 'MDE_TMDS',#State 15
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'TPDS_TMDE',#State 22
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'None',#State 23
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'MDE_TMDS',#State 24
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'TPDS_TMDE',#State 16
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'None',#State 17
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'MDE_TMDS',#State 18
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'TPDS_TMDE',#State 13
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 14
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'MDE_TMDS',#State 15
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'TPDS_TMDE',#State 13
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 14
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'MDE_TMDS',#State 15
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="UP" ~ 'TPDS_TMDE',#State 19
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="UP" ~ 'TMDE',#State 20
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="UP" ~ 'MDE',#State 21
                                   
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TPDS',#State 25
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TMDS',#State 26
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TPDE_TMDS',#State 27
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'TPDS_TMDE',#State 31
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'None',#State 32
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'TPDE_TMDS',#State 33
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'TPDS_TMDE',#State 28
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'None',#State 29
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'TPDE_TMDS',#State 30
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'TPDS_TMDE',#State 25
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 26
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'TPDE_TMDS',#State 27
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'TPDS_TMDE',#State 34
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 35
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'TPDE_TMDS',#State 36
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="UP" ~ 'TPDS_TMDE',#State 34
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="UP" ~ 'TMDE',#State 35
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="UP" ~ 'TPDE',#State 36
                                   TRUE ~ 'NA'))%>%
    mutate(RG3_Protein = case_when(BG_Method =="FALSE"~ 'Background = FALSE',
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="DOWN" ~ 'MDS',#State 1
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TMDS',#State 2
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TMDS',#State 3
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'MDS',#State 7
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'None',#State 8
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'TPDE',#State 9
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'None',#State 4
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'None',#State 5
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'None',#State 6
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'MDS',#State 10
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 11
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'TMDS',#State 12
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'TMDE',#State 10
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 11
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'TPDE',#State 12
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="UP" ~ 'TMDE',#State 10
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="UP" ~ 'TMDE',#State 11
                                   methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="UP" ~ 'TPDE',#State 12
                                   
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TPDS',#State 13
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TMDS',#State 14
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TMDS',#State 15
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'TPDS',#State 22
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'None',#State 23
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'MDE',#State 24
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'None',#State 16
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'None',#State 17
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'None',#State 18
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'TPDS',#State 13
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 14
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'TMDS',#State 15
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'TMDE',#State 13
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 14
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'MDE',#State 15
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="UP" ~ 'TMDE',#State 19
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="UP" ~ 'TMDE',#State 20
                                   methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="UP" ~ 'MDE',#State 21
                                   
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TPDS',#State 25
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TMDS',#State 26
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TMDS',#State 27
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'TPDS',#State 31
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'None',#State 32
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'TPDE',#State 33
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'None',#State 28
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'None',#State 29
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'None',#State 30
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'TPDS',#State 25
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 26
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'TMDS',#State 27
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'TMDE',#State 34
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 35
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'TPDE',#State 36
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="UP" ~ 'TMDE',#State 34
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="UP" ~ 'TMDE',#State 35
                                   methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="UP" ~ 'TPDE',#State 36
                                   TRUE ~ 'NA'))%>%
    mutate(RG4_Detection = case_when(BG_Method =="FALSE"~ 'Background = FALSE',
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="DOWN" ~ 'MDS',#State 1
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TMDS',#State 2
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TPDE_TMDS',#State 3
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'MDS',#State 7
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'None',#State 8
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'TPDE',#State 9
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'MDS_TMDE',#State 4 #####
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'None',#State 5
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'TPDE_TMDS',#State 6 #####
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'MDS',#State 10
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 11
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'TPDE_TMDS',#State 12
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'MDS_TMDE',#State 10
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 11
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'TPDE',#State 12
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="UP" ~ 'MDS_TMDE',#State 10
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="UP" ~ 'TMDE',#State 11
                                     methDF_Cutoff=="Hypermethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="UP" ~ 'TPDE',#State 12
                                     
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TPDS',#State 13
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TMDS',#State 14
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="DOWN" ~ 'MDE_TMDS',#State 15
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'TPDS_TMDE',#State 22
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'None',#State 23
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'MDE',#State 24
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'TPDS_TMDE',#State 16 ######
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'None',#State 17
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'MDE_TMDS',#State 18  ######
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'TPDS',#State 13
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 14
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'MDE_TMDS',#State 15
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'TPDS_TMDE',#State 13
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 14
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'MDE',#State 15
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="UP" ~ 'TPDS_TMDE',#State 19
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="UP" ~ 'TMDE',#State 20
                                     methDF_Cutoff=="Hypomethylation" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="UP" ~ 'MDE',#State 21
                                     
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TPDS',#State 25
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TMDS',#State 26
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="DOWN" ~ 'TPDE_TMDS',#State 27
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'TPDS',#State 31
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'None',#State 32
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Detected" ~ 'TPDE',#State 33
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'TPDS_TMDE',#State 28 ######
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'None',#State 29
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Not Significant" ~ 'TPDE_TMDS',#State 30 ######
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'TPDS',#State 25
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'None',#State 26
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Negative" ~ 'TPDE_TMDS',#State 27
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'TPDS_TMDE',#State 34
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'None',#State 35
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="Significant Positive" ~ 'TPDE',#State 36
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="DOWN" & proteinDF_Cutoff_Specific=="UP" ~ 'TPDS_TMDE',#State 34
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="No Change" & proteinDF_Cutoff_Specific=="UP" ~ 'TMDE',#State 35
                                     methDF_Cutoff=="No Change" & rnaDF_Cutoff=="UP" & proteinDF_Cutoff_Specific=="UP" ~ 'TPDE',#State 36
                                     TRUE ~ 'NA'))
  
  #Safe the DF and return the groupings
  ##RCM DF (Merged InputDF filtered for background with assigned SiRcle cluster names)
  MergeDF_Select1 <- MergeDF[, c("geneID", "methDF_Detected","methDF_ValueCol","methDF_PadjCol","methDF_Cutoff", "rnaDF_Detected","rnaDF_ValueCol","rnaDF_PadjCol","rnaDF_Cutoff", "proteinDF_Detected", "proteinDF_ValueCol","proteinDF_PadjCol","proteinDF_Cutoff", "proteinDF_Cutoff_Specific", "BG_Method", "RG1_All", "RG2_Changes", "RG3_Protein", "RG4_Detection")]
  
  proteinValueCol_Unique<-paste("proteinDF_",proteinValueCol)
  proteinPadjCol_Unique <-paste("proteinDF_",proteinPadjCol)
  rnaValueCol_Unique<-paste("rnaDF_",rnaValueCol)
  rnaPadjCol_Unique <-paste("rnaDF_",rnaPadjCol)
  methValueCol_Unique<-paste("methDF_",methValueCol)
  methPadjCol_Unique <-paste("methDF_",methPadjCol)
  
  MergeDF_Select2<- subset(MergeDF, select=-c(methDF_Detected,methDF_Cutoff, rnaDF_Detected,rnaDF_Cutoff, proteinDF_Detected,proteinDF_Cutoff, proteinDF_Cutoff_Specific, BG_Method, RG1_All, RG2_Changes, RG3_Protein, RG4_Detection))%>%
    rename(!!proteinValueCol_Unique :="proteinDF_ValueCol",#This syntax is needed since paste(geneID)="geneID" is not working in dyplr
           !!proteinPadjCol_Unique :="proteinDF_PadjCol",
           !!rnaValueCol_Unique :="rnaDF_ValueCol",
           !!rnaPadjCol_Unique :="rnaDF_PadjCol",
           !!methValueCol_Unique :="methDF_ValueCol",
           !!methPadjCol_Unique :="methDF_PadjCol")
  
  MergeDF_Rearrange <- merge(MergeDF_Select1, MergeDF_Select2, by="geneID")
  
  write.csv(MergeDF_Rearrange, OutputFileName, row.names = FALSE)
  return(MergeDF_Rearrange)
  
  ##Summary SiRCle clusters (number of genes assigned to each SiRCle cluster in each grouping)
  ClusterSummary_RG1 <- MergeDF_Rearrange[,c("geneID", "RG1_All")]%>%
    count(RG1_All, name="Number of Genes")%>%
    rename("SiRCle cluster Name"= "RG1_All")
  ClusterSummary_RG1$`Regulation Grouping` <- "RG1_All"
  
  ClusterSummary_RG2 <- MergeDF_Rearrange[,c("geneID", "RG2_Changes")]%>%
    count(RG2_Changes, name="Number of Genes")%>%
    rename("SiRCle cluster Name"= "RG2_Changes")
  ClusterSummary_RG2$`Regulation Grouping` <- "RG2_Changes"
  
  ClusterSummary_RG3 <- MergeDF_Rearrange[,c("geneID", "RG3_Protein")]%>%
    count(RG3_Protein, name="Number of Genes")%>%
    rename("SiRCle cluster Name"= "RG3_Protein")
  ClusterSummary_RG3$`Regulation Grouping` <- "RG3_Protein"
  
  ClusterSummary_RG4 <- MergeDF_Rearrange[,c("geneID","RG4_Detection")]%>%
    count(RG4_Detection, name="Number of Genes")%>%
    rename("SiRCle cluster Name"= "RG4_Detection")
  ClusterSummary_RG4$`Regulation Grouping` <- "RG4_Detection"
  
  ClusterSummary <- rbind(ClusterSummary_RG1, ClusterSummary_RG2,ClusterSummary_RG3 , ClusterSummary_RG4)
  ClusterSummary <- ClusterSummary[,c(3,1,2)]
  
  write.csv(ClusterSummary, paste("Summary_",OutputFileName), row.names = FALSE)
}

