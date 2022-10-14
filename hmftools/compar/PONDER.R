#!/usr/bin/Rscript
# PONDER: To think in detail, in this case: to extract data from COMPAR
# Author: Teoman Deger
# Version 0.1, 14-10-2022
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
  COMPARname   <- deparse(substitute(x))
  noValue      <- x[x$MismatchType != "VALUE",]
  #from data-groups (NEW_ONLY, REF_ONLY, VALUE, SHARED)
  comparGroups <- noValue[!duplicated(noValue$MismatchType),colnames(noValue) == "MismatchType"]
  #make export-tables, SNPs
  groups_SNP <- c("A>C", "A>G", "A>T", "C>G", "C>T", "G>T", "SNP_total")
  countsSNP <- as.data.frame(matrix(nrow=length(comparGroups), ncol=length(groups_SNP)))
  colnames(countsSNP) <- groups_SNP
  rownames(countsSNP) <- comparGroups
  #make export-tables, SV
  groups_SV <- c("SV_BND", "SV_DEL", "SV_DUP", "SV_INV", "SV_SGL", "SV_INS", "SV_total")
  countsSV <- as.data.frame(matrix(nrow=length(comparGroups), ncol=length(groups_SV)))
  colnames(countsSV) <- groups_SV
  rownames(countsSV) <- comparGroups
  #make export-tables, countsINDELMNP
  groups_MNPINDEL <- c("MNP-2", "MNP-3", "MNP-4+","MNP_total", "INDELs")
  countsMNPINDEL <- as.data.frame(matrix(nrow=length(comparGroups), ncol=length(groups_MNPINDEL)))
  colnames(countsMNPINDEL) <- groups_MNPINDEL
  rownames(countsMNPINDEL) <- comparGroups
  #make export-tables, countsCN
  groups_CN <- c("CN (0-2)", "CN (2-4)", "CN (4+)", "Total")
  countsCN <- as.data.frame(matrix(nrow=length(comparGroups), ncol=length(groups_CN)))
  colnames(countsCN) <- groups_CN
  rownames(countsCN) <- comparGroups
  #for each data-group do:
  for (i in comparGroups){
    #extract data from group only
    selection  <- noValue[noValue$MismatchType == i,]
    #extract Somatic Variants only
    selectionSomVars <- selection[selection$Category == "SOMATIC_VARIANT",]
    #extract SNPs only
    SNPgrep      <- selectionSomVars[grep("SNP",selectionSomVars$Key),colnames(selectionSomVars) == "Key"]
    #group per SNP-type
    SNPtypes <- data.frame(do.call('rbind', strsplit(as.character(SNPgrep), " ", fixed=TRUE)))
    #combine rev-comp SNP-types
    SNPtypes[,2][SNPtypes[,2] == "T>G"] <- "A>C"
    SNPtypes[,2][SNPtypes[,2] == "T>C"] <- "A>G"
    SNPtypes[,2][SNPtypes[,2] == "T>A"] <- "A>T"
    SNPtypes[,2][SNPtypes[,2] == "G>C"] <- "C>G"
    SNPtypes[,2][SNPtypes[,2] == "G>A"] <- "C>T"
    SNPtypes[,2][SNPtypes[,2] == "C>A"] <- "G>T"
    A_Csnp <- sum(SNPtypes[,2] == "A>C")
    A_Gsnp <- sum(SNPtypes[,2] == "A>G")
    A_Tsnp <- sum(SNPtypes[,2] == "A>T")
    C_Gsnp <- sum(SNPtypes[,2] == "C>G")
    C_Tsnp <- sum(SNPtypes[,2] == "C>T")
    G_Tsnp <- sum(SNPtypes[,2] == "G>T")
    snp_Total <- sum(A_Csnp, A_Gsnp, A_Tsnp, C_Gsnp, C_Tsnp, G_Tsnp)
    countsSNP[i,]<- c(A_Csnp, A_Gsnp, A_Tsnp, C_Gsnp, C_Tsnp, G_Tsnp, snp_Total)
    #extract MNPs only
    MNPgrep <- selectionSomVars[grep("MNP",selectionSomVars$Key),colnames(selectionSomVars) == "Key"]
    #extract MNP-info only
    MNPtypes <- as.character(data.frame(do.call('rbind', strsplit(as.character(MNPgrep), " ", fixed=TRUE)))[,2])
    #get lengths in categories (2,3,more)
    MNP_2    <- sum((nchar(MNPtypes)-1)/2 == 2)
    MNP_3    <- sum((nchar(MNPtypes)-1)/2 == 3)
    MNP_more <- sum((nchar(MNPtypes)-1)/2 >3)
    MNP_total <- sum(MNP_2, MNP_3, MNP_more)
    #extract number of indels
    INDELs            <- length(selectionSomVars[grep("INDEL",selectionSomVars$Key),colnames(selectionSomVars) == "Key"])
    countsMNPINDEL[i,]<- c(MNP_2, MNP_3, MNP_more, MNP_total, INDELs)
    #extract SV-events only
    selection_SV <- selection[selection$Category == "DISRUPTION",]
    selection_SV$Key <- as.character(selection_SV$Key)
    #get counts per event type
    SV_BND   <- sum(do.call('rbind', strsplit(do.call('rbind', strsplit(selection_SV$Key, "_"))[,2], " "))[,1] == "BND")
    SV_DEL   <- sum(do.call('rbind', strsplit(do.call('rbind', strsplit(selection_SV$Key, "_"))[,2], " "))[,1] == "DEL")
    SV_DUP   <- sum(do.call('rbind', strsplit(do.call('rbind', strsplit(selection_SV$Key, "_"))[,2], " "))[,1] == "DUP")
    SV_INV   <- sum(do.call('rbind', strsplit(do.call('rbind', strsplit(selection_SV$Key, "_"))[,2], " "))[,1] == "INV")
    SV_SGL   <- sum(do.call('rbind', strsplit(do.call('rbind', strsplit(selection_SV$Key, "_"))[,2], " "))[,1] == "SGL")
    SV_INS   <- sum(do.call('rbind', strsplit(do.call('rbind', strsplit(selection_SV$Key, "_"))[,2], " "))[,1] == "INS")
    SV_total <- sum(SV_BND, SV_DEL, SV_DUP, SV_INV, SV_SGL, SV_INS)
    countsSV[i,] <- c(SV_BND, SV_DEL, SV_DUP, SV_INV, SV_SGL, SV_INS, SV_total)
    #extract copynumbers only
    selectionCN  <- selectionSomVars <- selection[selection$Category == "COPY_NUMBER",]
    selectionCN$AllValues <- as.character(selectionCN$AllValues)
    CNcounts     <- as.numeric(do.call('rbind', strsplit(do.call('rbind', strsplit(selectionCN$AllValues, ";"))[,1], "="))[,2])
    CN_loss  <- sum(CNcounts <= 2)
    CN_gain  <- sum(CNcounts > 2 & CNcounts <= 4)
    CN_amp   <- sum(CNcounts > 4)
    CN_total <- sum(CN_loss, CN_gain, CN_amp)
    countsCN[i,] <- c(CN_loss, CN_gain, CN_amp, CN_total)
    }
  #export all
  outputList <- list(name,countsSNP,countsMNPINDEL, countsCN, countsSV)
  print(outputList)#, file = paste(deparse(substitute(x)),"_",i))
  if (args[2] == "EXPORT"){
    capture.output(outputList, file = paste0(name,"_COMPARed_Output"))
  }
}
analysisComparOutput(input, args[2])