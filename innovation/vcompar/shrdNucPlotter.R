#!/usr/bin/Rscript
# Author: Teoman Deger
# -----------------------./
args = commandArgs(trailingOnly=TRUE)

#sanity checks
if (length(args)!=7) {
  stop("plotter requires 7 arguments, not all were found", call.=FALSE)
}

name1<-args[1]
name2<-args[2]
name3<-args[3]
name4<-args[4]

TypeName<-args[5]
refGen<-args[6]
setwd(args[7])

if (refGen == "hg19"){
  genomeVer <- "BSgenome.Hsapiens.UCSC.hg19"
} else if (refGen == "hg38"){
  genomeVer <- "BSgenome.Hsapiens.UCSC.hg38"
} else {
  print("provide hg19 or hg38")
  q()
}

suppressMessages(library(genomeVer, character.only = TRUE))
suppressMessages(library(MutationalPatterns))

read12_1 <- read_vcfs_as_granges("only_12_1.vcf.gz", paste0("Unique_old_",name1), genomeVer)
read12_2 <- read_vcfs_as_granges("only_12_2.vcf.gz", paste0("Unique_old_",name2), genomeVer)
read34_3 <- read_vcfs_as_granges("only_34_3.vcf.gz", paste0("Unique_new_",name3), genomeVer)
read34_4 <- read_vcfs_as_granges("only_34_4.vcf.gz", paste0("Unique_new_",name4), genomeVer)
readShr1 <- read_vcfs_as_granges("shared_all_1.vcf.gz", paste0("Shared_all_",name1), genomeVer)
readShr2 <- read_vcfs_as_granges("shared_all_2.vcf.gz", paste0("Shared_all_",name2), genomeVer)
readShr3 <- read_vcfs_as_granges("shared_all_3.vcf.gz", paste0("Shared_all_",name3), genomeVer)
readShr4 <- read_vcfs_as_granges("shared_all_4.vcf.gz", paste0("Shared_all_",name4), genomeVer)

TriNucleotidePlotter <- function(ref,new,ref_genome,firstName,secondName){
  refmm <- mut_matrix(ref, ref_genome)
  newmm <- mut_matrix(new, ref_genome)
  plot_compare_profiles(refmm[,1], newmm[,1], profile_name=c(firstName,secondName))
}

out<-TriNucleotidePlotter(read12_1, read12_2, genomeVer, name1, name2)
png(paste0("MutPat_OldUniques_",name1,"_",name2,TypeName,".png"), height = 768, width = 1024)
plot(out)
dev.off()

out<-TriNucleotidePlotter(read34_3, read34_4, genomeVer, name3, name4)
png(paste0("MutPat_NewUniques_",name3,"_",name4,TypeName,".png"), height = 768, width = 1024)
plot(out)
dev.off()

out<-TriNucleotidePlotter(readShr1, readShr2, genomeVer, name1, name2)
png(paste0("MutPat_AllShared_",name1,"_",name2 ,TypeName,".png"), height = 768, width = 1024)
plot(out)
dev.off()

out<-TriNucleotidePlotter(readShr3, readShr4, genomeVer, name3, name4)
png(paste0("MutPat_AllShared_",name3,"_",name4 ,TypeName,".png"), height = 768, width = 1024)
plot(out)
dev.off()
