#!/usr/bin/Rscript
# Author: Teoman Deger
# -----------------------./
args = commandArgs(trailingOnly=TRUE)


#sanity checks
if (length(args)!=4) {
  stop("plotter requires 4 files, not all were found", call.=FALSE)
} 

setwd(args[4])
TypeName<-args[3]

myFiles<-sort(list.files()[grepl(list.files(), pattern = "000") & grepl(list.files(), pattern = ".vcf$")])


for(i in 1:length(myFiles)){
  example <-read.table(myFiles[i])
  colnames(example) <-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","REF","TUM")
  formatS <- strsplit(as.character(example$FORMAT[1]),":")[[1]]
  colV <- which(formatS == "AF")
  colN <- sapply(strsplit(as.character(example$TUM[1]),":"), length)
  vafs <- as.numeric(matrix(unlist(strsplit(as.character(example$TUM),":")), ncol=colN, byrow=TRUE)[,colV])
  if (i == 1){
    png(paste0(args[i],"_",TypeName,"_Vaf.png"), height = 768, width = 1024)
    hist(vafs, breaks=500, main = paste0(args[i],"_only"))
    dev.off()
  } else if (i == 2){
    png(paste0(args[i],"_",TypeName,"_Vaf.png"), height = 768, width = 1024)
    hist(vafs, breaks=500, main = paste0(args[i],"_only"))
    dev.off()
  } else if (i == 3){
    example2 <-read.table(myFiles[4])
    colnames(example2) <-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","REF","TUM")
    colN <- sapply(strsplit(as.character(example2$TUM[1]),":"), length)
    vafs2 <- as.numeric(matrix(unlist(strsplit(as.character(example2$TUM), ":")), ncol=colN, byrow=TRUE)[,4])
    minimat <-cbind(vafs,vafs2)
    vafs_avg <- rowMeans(minimat)
    rm(minimat)
    png(paste0(args[1],"_",args[2],"_shared_",TypeName,"_Vaf.png"), height = 768, width = 1024)
    hist(vafs_avg, breaks=500, main = paste0(args[1],"&",args[2],"_shared"))
    dev.off()
  } else if (i == 4){
    example2 <-read.table(myFiles[3])
    colnames(example2) <-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","REF","TUM")
    colN <- sapply(strsplit(as.character(example2$TUM[1]),":"), length)
    vafs2 <- as.numeric(matrix(unlist(strsplit(as.character(example2$TUM), ":")), ncol=colN, byrow=TRUE)[,4])
    png(paste0(args[1],"_",args[2],"_shared_",TypeName,"Vaf_Scatter.png"), height = 768, width = 1024)
    plot(vafs2,vafs, pch=20, cex=0.1, xlab=paste0("Vafs ",args[2]), ylab=paste0("Vafs ",args[1]), main=paste0(TypeName,"Scatterplot shareds: ",args[1]," and ",args[2]))
    dev.off()
  }
  #hist(vafs, breaks=500, main = myFiles[i])
}
