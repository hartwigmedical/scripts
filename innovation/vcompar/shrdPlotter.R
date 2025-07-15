#!/usr/bin/Rscript
# Author: Teoman Deger
# -----------------------./
args = commandArgs(trailingOnly=TRUE)

#sanity checks
if (length(args)!=2) {
  stop("plotter requires 2 argument, not all were found", call.=FALSE)
}

TypeName<-args[1]
setwd(args[2])

myFiles<-sort(list.files()[grepl(list.files(), pattern = "000") & grepl(list.files(), pattern = ".vcf$")])

#plot vafs (1/2)only
only_12_1 <- read.table("only_12_1.vcf.gz")
colnames(only_12_1)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","REF","TUM")
colN <- sapply(strsplit(as.character(only_12_1$TUM[1]),":"), length)
only_12_1_vafs <- as.numeric(matrix(unlist(strsplit(only_12_1$TUM, ":")), ncol=colN, byrow=TRUE)[,4])
#
only_12_2 <- read.table("only_12_2.vcf.gz")
colnames(only_12_2)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","REF","TUM")
colN <- sapply(strsplit(as.character(only_12_2$TUM[1]),":"), length)
only_12_2_vafs <- as.numeric(matrix(unlist(strsplit(only_12_2$TUM, ":")), ncol=colN, byrow=TRUE)[,4])
#
minimat <-cbind(only_12_1_vafs, only_12_2_vafs)
vafs_avg <- rowMeans(minimat)
rm(minimat)
png(paste0(TypeName,"_12_only.png"), height = 768, width = 1024)
hist(vafs_avg, breaks=500, main ="VAFs - Uniques, old batch")
invisible(dev.off())

#plot vafs (3/4)only
only_34_3 <- read.table("only_34_3.vcf.gz")
colnames(only_34_3)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","REF","TUM")
colN <- sapply(strsplit(as.character(only_34_3$TUM[1]),":"), length)
only_34_3_vafs <- as.numeric(matrix(unlist(strsplit(only_34_3$TUM, ":")), ncol=colN, byrow=TRUE)[,4])
#
only_34_4 <- read.table("only_34_4.vcf.gz")
colnames(only_34_4)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","REF","TUM")
colN <- sapply(strsplit(as.character(only_34_4$TUM[1]),":"), length)
only_34_4_vafs <- as.numeric(matrix(unlist(strsplit(only_34_4$TUM, ":")), ncol=colN, byrow=TRUE)[,4])
#
minimat <-cbind(only_34_3_vafs, only_34_4_vafs)
vafs_avg <- rowMeans(minimat)
rm(minimat)
png(paste0(TypeName,"_34_only.png"), height = 768, width = 1024)
hist(vafs_avg, breaks=500, main ="VAFs - Uniques, new batch")
invisible(dev.off())

#shared vafs
shared_all_1 <- read.table("shared_all_1.vcf.gz")
shared_all_2 <- read.table("shared_all_2.vcf.gz")
shared_all_3 <- read.table("shared_all_3.vcf.gz")
shared_all_4 <- read.table("shared_all_4.vcf.gz")
colnames(shared_all_1)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","REF","TUM")
colnames(shared_all_2)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","REF","TUM")
colnames(shared_all_3)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","REF","TUM")
colnames(shared_all_4)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","REF","TUM")
colN1 <- sapply(strsplit(as.character(shared_all_1$TUM[1]),":"), length)
colN2 <- sapply(strsplit(as.character(shared_all_2$TUM[1]),":"), length)
colN3 <- sapply(strsplit(as.character(shared_all_3$TUM[1]),":"), length)
colN4 <- sapply(strsplit(as.character(shared_all_4$TUM[1]),":"), length)
shared_vafs_1 <- as.numeric(matrix(unlist(strsplit(shared_all_1$TUM, ":")), ncol=colN1, byrow=TRUE)[,4])
shared_vafs_2 <- as.numeric(matrix(unlist(strsplit(shared_all_2$TUM, ":")), ncol=colN2, byrow=TRUE)[,4])
shared_vafs_3 <- as.numeric(matrix(unlist(strsplit(shared_all_3$TUM, ":")), ncol=colN3, byrow=TRUE)[,4])
shared_vafs_4 <- as.numeric(matrix(unlist(strsplit(shared_all_4$TUM, ":")), ncol=colN4, byrow=TRUE)[,4])
minimat<-cbind(shared_vafs_1,shared_vafs_2,shared_vafs_3,shared_vafs_4)
shared_vafs_old<-rowMeans(minimat[,c(1:2)])
shared_vafs_new<-rowMeans(minimat[,c(3:4)])
shared_vafs_all<-rowMeans(minimat[,c(1:4)])
rm(minimat)
png(paste0(TypeName,"_Vaf_shared_all.png"), height = 768, width = 1024)
hist(shared_vafs_all, breaks=500, main ="VAFs - Shared, all")
invisible(dev.off())
png(paste0(TypeName,"_Scatter_shared_all.png"), height = 768, width = 1024)
plot(shared_vafs_old,shared_vafs_new, pch=20, cex=0.1, xlab="Old batch", ylab="New batch", main=paste0(TypeName," Scatterplot shareds"))
invisible(dev.off())


