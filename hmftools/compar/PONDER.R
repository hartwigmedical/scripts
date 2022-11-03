#!/usr/bin/Rscript
# PONDER: To think in detail, in this case: to extract data from COMPAR
# Author: Teoman Deger
# Version 0.3, 03-11-2022
# USE: Rscript PONDER.R COMPARoutput.csv
# -----------------------./
args = commandArgs(trailingOnly=TRUE)

#sanity checks
if (length(args)==0) {
  stop("Give a COMPAR output-file as first argument, add EXPORT as second argument to save output", call.=FALSE)
} else if (length(args) < 2) {
  args[2] = "placeholder"
}

input <- read.csv(args[1])
name <- substr(args[1],1,nchar(args[1])-4)

analysisComparOutput <- function(x, export = FALSE){
  #Remove all VALUE-type 'mismatches' from DF
  name <- deparse(substitute(x))
  noValue      <- x[x$MismatchType != "VALUE",]
  #from data-groups (NEW_ONLY, REF_ONLY, VALUE, SHARED)
  comparGroups <- noValue[!duplicated(noValue$MismatchType),colnames(noValue) == "MismatchType"]
  #make export-tables, SNPs_Somatic
  groups_SNP <- c("A>C", "A>G", "A>T", "C>G", "C>T", "G>T")
  countsSNP  <- as.data.frame(matrix(nrow=length(comparGroups), ncol=length(groups_SNP)))
  colnames(countsSNP) <- groups_SNP
  rownames(countsSNP) <- comparGroups
  #make export-tables, SNPs_Germline
  countsSNP_GL  <- as.data.frame(matrix(nrow=length(comparGroups), ncol=length(groups_SNP)))
  colnames(countsSNP_GL) <- groups_SNP
  rownames(countsSNP_GL) <- comparGroups
  #make export-tables, StrucVars
  groups_SV <- c("SV_BND", "SV_DEL", "SV_DUP", "SV_INV", "SV_SGL", "SV_INS", "SV_INF")
  countsSV  <- as.data.frame(matrix(nrow=length(comparGroups), ncol=length(groups_SV)))
  colnames(countsSV) <- groups_SV
  rownames(countsSV) <- comparGroups
  #make export-tables, GLVars
  groups_GLV <- c("GLV_BND", "GLV_DEL", "GLV_DUP", "GLV_INV", "GLV_SGL", "GLV_INS", "GLV_INF")
  countsGLV  <- as.data.frame(matrix(nrow=length(comparGroups), ncol=length(groups_GLV)))
  colnames(countsGLV) <- groups_GLV
  rownames(countsGLV) <- comparGroups
  #make export-tables, countsINDELMNP_Somatic
  groups_MNPINDEL <- c("MNP-2", "MNP-3", "MNP-4+","INDELs")
  countsMNPINDEL  <- as.data.frame(matrix(nrow=length(comparGroups), ncol=length(groups_MNPINDEL)))
  colnames(countsMNPINDEL) <- groups_MNPINDEL
  rownames(countsMNPINDEL) <- comparGroups
  #make export-tables, countsINDELMNP_Germline
  countsMNPINDEL_GL  <- as.data.frame(matrix(nrow=length(comparGroups), ncol=length(groups_MNPINDEL)))
  colnames(countsMNPINDEL_GL) <- groups_MNPINDEL
  rownames(countsMNPINDEL_GL) <- comparGroups
  #make export-tables, countsCN
  groups_CN <- c("CN(0-2)", "CN(2-4)", "CN(4+)")
  countsCN  <- as.data.frame(matrix(nrow=length(comparGroups), ncol=length(groups_CN)))
  colnames(countsCN) <- groups_CN
  rownames(countsCN) <- comparGroups
  #make export-tables, totals
  groups_total <-c("SomatVars","GermlineVars", "StrucVars", "GL-SVs", "CNs", "Driver")
  countsTotal  <- as.data.frame(matrix(nrow=length(comparGroups), ncol=length(groups_total)))
  colnames(countsTotal) <- groups_total
  rownames(countsTotal) <- comparGroups
  #for each data-group do:
  for (i in comparGroups){
    #extract data from group only
    #i = "REF_ONLY"
    selection  <- noValue[noValue$MismatchType == i,]
    #extract Drivers only
    selectionDrivers <- nrow(selection[selection$Category == "DRIVER",])
    #extract Somatic Variants only
    selectionSomVars <- selection[selection$Category == "SOMATIC_VARIANT",]
    #extract SNPs only
    SNPgrep      <- selectionSomVars[grep("SNP",selectionSomVars$Key),colnames(selectionSomVars) == "Key"]
    #group per SNP-type
    SNPtypes <- data.frame(do.call('rbind', strsplit(as.character(SNPgrep), " ", fixed=TRUE)))
    #combine rev-comp SNP-types
    A_Csnp <- sum(SNPtypes[,2] == "A>C", SNPtypes[,2] == "T>G")
    A_Gsnp <- sum(SNPtypes[,2] == "A>G", SNPtypes[,2] == "T>C")
    A_Tsnp <- sum(SNPtypes[,2] == "A>T", SNPtypes[,2] == "T>A")
    C_Gsnp <- sum(SNPtypes[,2] == "C>G", SNPtypes[,2] == "G>C")
    C_Tsnp <- sum(SNPtypes[,2] == "C>T", SNPtypes[,2] == "G>A")
    G_Tsnp <- sum(SNPtypes[,2] == "G>T", SNPtypes[,2] == "C>A")
    countsSNP[i,] <- c(A_Csnp, A_Gsnp, A_Tsnp, C_Gsnp, C_Tsnp, G_Tsnp)
    #extract MNPs only
    MNPgrep  <- selectionSomVars[grep("MNP",selectionSomVars$Key),colnames(selectionSomVars) == "Key"]
    #extract MNP-info only
    if (length(MNPgrep) >=  1){
      MNPtypes <- as.character(data.frame(do.call('rbind', strsplit(as.character(MNPgrep), " ", fixed=TRUE)))[,2])
      #get lengths in categories (2,3,more)
      MNP_2    <- sum((nchar(MNPtypes)-1)/2 == 2)
      MNP_3    <- sum((nchar(MNPtypes)-1)/2 == 3)
      MNP_more <- sum((nchar(MNPtypes)-1)/2 >3)
    } else {
      MNP_2    <- 0
      MNP_3    <- 0
      MNP_more <- 0
    }
    #extract number of indels
    INDELs            <- length(selectionSomVars[grep("INDEL",selectionSomVars$Key),colnames(selectionSomVars) == "Key"])
    countsMNPINDEL[i,]<- c(MNP_2, MNP_3, MNP_more,INDELs)
    #extract Somatic Variants only
    selectionGLVars <- selection[selection$Category == "GERMLINE_VARIANT",]
    if (nrow(selectionGLVars) >= 1){
      #extract SNPs only
      SNPgrep_GL      <- selectionGLVars[grep("SNP",selectionGLVars$Key),colnames(selectionGLVars) == "Key"]
      #group per SNP-type
      SNPtypes_GL <- data.frame(do.call('rbind', strsplit(as.character(SNPgrep_GL), " ", fixed=TRUE)))
      #combine rev-comp SNP-types
      A_Csnp_GL <- sum(SNPtypes_GL[,2] == "A>C", SNPtypes_GL[,2] == "T>G")
      A_Gsnp_GL <- sum(SNPtypes_GL[,2] == "A>G", SNPtypes_GL[,2] == "T>C")
      A_Tsnp_GL <- sum(SNPtypes_GL[,2] == "A>T", SNPtypes_GL[,2] == "T>A")
      C_Gsnp_GL <- sum(SNPtypes_GL[,2] == "C>G", SNPtypes_GL[,2] == "G>C")
      C_Tsnp_GL <- sum(SNPtypes_GL[,2] == "C>T", SNPtypes_GL[,2] == "G>A")
      G_Tsnp_GL <- sum(SNPtypes_GL[,2] == "G>T", SNPtypes_GL[,2] == "C>A")
      countsSNP_GL[i,] <- c(A_Csnp_GL, A_Gsnp_GL, A_Tsnp_GL, C_Gsnp_GL, C_Tsnp_GL, G_Tsnp_GL)
    } else {
      countsSNP_GL[i,] <- c(0, 0, 0, 0, 0, 0)
    }
    #extract MNPs only
    MNPgrep_GL  <- selectionGLVars[grep("MNP",selectionGLVars$Key),colnames(selectionGLVars) == "Key"]
    #extract MNP-info only
    if (length(MNPgrep_GL) >=  1){
      MNPtypes_GL <- as.character(data.frame(do.call('rbind', strsplit(as.character(MNPgrep_GL), " ", fixed=TRUE)))[,2])
      #get lengths in categories (2,3,more)
      MNP_2_GL    <- sum((nchar(MNPtypes_GL)-1)/2 == 2)
      MNP_3_GL    <- sum((nchar(MNPtypes_GL)-1)/2 == 3)
      MNP_more_GL <- sum((nchar(MNPtypes_GL)-1)/2 >3)
    } else {
      MNP_2_GL    <- 0
      MNP_3_GL    <- 0
      MNP_more_GL <- 0
    }
    #extract number of indels
    INDELs_GL            <- length(selectionGLVars[grep("INDEL",selectionGLVars$Key),colnames(selectionGLVars) == "Key"])
    countsMNPINDEL_GL[i,]<- c(MNP_2_GL, MNP_3_GL, MNP_more_GL,INDELs_GL)
    #extract SV-events only
    selection_SV     <- selection[selection$Category == "DISRUPTION",]
    selection_SV$Key <- as.character(selection_SV$Key)
    #get counts per event type
    SV_BND       <- sum(do.call('rbind', strsplit(do.call('rbind', strsplit(selection_SV$Key, "_"))[,2], " "))[,1] == "BND")
    SV_DEL       <- sum(do.call('rbind', strsplit(do.call('rbind', strsplit(selection_SV$Key, "_"))[,2], " "))[,1] == "DEL")
    SV_DUP       <- sum(do.call('rbind', strsplit(do.call('rbind', strsplit(selection_SV$Key, "_"))[,2], " "))[,1] == "DUP")
    SV_INV       <- sum(do.call('rbind', strsplit(do.call('rbind', strsplit(selection_SV$Key, "_"))[,2], " "))[,1] == "INV")
    SV_SGL       <- sum(do.call('rbind', strsplit(do.call('rbind', strsplit(selection_SV$Key, "_"))[,2], " "))[,1] == "SGL")
    SV_INS       <- sum(do.call('rbind', strsplit(do.call('rbind', strsplit(selection_SV$Key, "_"))[,2], " "))[,1] == "INS")
    SV_INF       <- sum(do.call('rbind', strsplit(do.call('rbind', strsplit(selection_SV$Key, "_"))[,2], " "))[,1] == "INF")
    countsSV[i,] <- c(SV_BND, SV_DEL, SV_DUP, SV_INV, SV_SGL, SV_INS, SV_INF)
    #extract GLV-events only
    selection_GLV     <- selection[selection$Category == "GERMLINE_SV",]
    selection_GLV$Key <- as.character(selection_GLV$Key)
    #get counts per event type
    GLV_BND       <- sum(substr(do.call('rbind', strsplit(selection_GLV$Key, ":"))[,2],1,3) == "BND")
    GLV_DEL       <- sum(substr(do.call('rbind', strsplit(selection_GLV$Key, ":"))[,2],1,3) == "DEL")
    GLV_DUP       <- sum(substr(do.call('rbind', strsplit(selection_GLV$Key, ":"))[,2],1,3) == "DUP")
    GLV_INV       <- sum(substr(do.call('rbind', strsplit(selection_GLV$Key, ":"))[,2],1,3) == "INV")
    GLV_SGL       <- sum(substr(do.call('rbind', strsplit(selection_GLV$Key, ":"))[,2],1,3) == "SGL")
    GLV_INS       <- sum(substr(do.call('rbind', strsplit(selection_GLV$Key, ":"))[,2],1,3) == "INS")
    GLV_INF       <- sum(substr(do.call('rbind', strsplit(selection_GLV$Key, ":"))[,2],1,3) == "INF")
    countsGLV[i,] <- c(GLV_BND, GLV_DEL, GLV_DUP, GLV_INV, GLV_SGL, GLV_INS, GLV_INF)
    #extract copynumbers only
    selectionCN  <- selectionSomVars <- selection[selection$Category == "COPY_NUMBER",]
    selectionCN$AllValues <- as.character(selectionCN$AllValues)
    CNcounts     <- as.numeric(do.call('rbind', strsplit(do.call('rbind', strsplit(selectionCN$AllValues, ";"))[,1], "="))[,2])
    CN_loss      <- sum(CNcounts <= 2)
    CN_gain      <- sum(CNcounts > 2 & CNcounts <= 4)
    CN_amp       <- sum(CNcounts > 4)
    countsCN[i,] <- c(CN_loss, CN_gain, CN_amp)
    #extract totals
    SomVars_Total   <- sum(as.integer(rowSums(countsSNP[i,])), MNP_2, MNP_3, MNP_more, INDELs)
    GLVars_Total    <- sum(as.integer(rowSums(countsSNP_GL[i,])), MNP_2_GL, MNP_3_GL, MNP_more_GL, INDELs_GL)
    SV_total        <- sum(as.integer(rowSums(countsSV[i,])))
    GLV_total       <- sum(as.integer(rowSums(countsGLV[i,])))
    CN_total        <- sum(as.integer(rowSums(countsCN[i,])))
    countsTotal[i,] <- c(SomVars_Total, GLVars_Total,  SV_total, GLV_total, CN_total, selectionDrivers)
  }
  #export all
  outputList <- list(name,countsTotal,"Somatic",countsSNP,countsMNPINDEL,"Germline",countsSNP_GL,countsMNPINDEL_GL,countsCN,countsSV,countsGLV)
  print(outputList)#, file = paste(deparse(substitute(x)),"_",i))
  if (export == "EXPORT"){
    capture.output(outputList, file = paste0(name,"_COMPARed_Output"))
  }
}
analysisComparOutput(input, args[2])