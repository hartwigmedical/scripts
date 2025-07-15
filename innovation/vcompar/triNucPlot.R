#!/usr/bin/Rscript
# Author: Teoman Deger
# -----------------------./
args = commandArgs(trailingOnly=TRUE)

#sanity checks
if (length(args)!=5) {
  stop("plotter requires 5 arguments, not all were found", call.=FALSE)
}

name1<-args[1]
name2<-args[2]
TypeName<-args[3]
refGen<-args[4]
setwd(args[5])

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

myFiles<-sort(list.files()[grepl(list.files(), pattern = "000") & grepl(list.files(), pattern = ".vcf$")])

readIn0 <- read_vcfs_as_granges(myFiles[1], paste0(name1,"_only"), genomeVer, predefined_dbs_mbs = TRUE)
readIn1 <- read_vcfs_as_granges(myFiles[2], paste0(name2,"_only"), genomeVer, predefined_dbs_mbs = TRUE)
readIn2 <- read_vcfs_as_granges(myFiles[3], paste0(name1,"_",name2,"_shared_1"), genomeVer, predefined_dbs_mbs = TRUE)
readIn3 <- read_vcfs_as_granges(myFiles[4], paste0(name1,"_",name2,"_shared_2"), genomeVer, predefined_dbs_mbs = TRUE)


TriNucleotidePlotter <- function(ref,new,ref_genome,firstName,secondName){
  refmm <- mut_matrix(ref, ref_genome)
  newmm <- mut_matrix(new, ref_genome)
  plot_compare_profiles(refmm[,1], newmm[,1], profile_name=c(firstName,secondName))
}

out<-TriNucleotidePlotter(readIn0, readIn1, genomeVer, name1, name2)
png(paste0("MutPat_Uniques_",name1,"_",name2,TypeName,".png"), height = 768, width = 1024)
plot(out)
dev.off()

out<-TriNucleotidePlotter(readIn2, readIn3, genomeVer, name1, name2)
png(paste0("MutPat_Shared_",name1,"_",name2 ,TypeName,".png"), height = 768, width = 1024)
plot(out)
dev.off()



#readMm <- mut_matrix(readIn, genomeVer, extension = 1)
#out <- plot_96_profile(readMm, colors = NA, ymax = 0.2, condensed = FALSE)
#png(paste0(args[i],"_",TypeName,"_MutPat.png"), height = 768, width = 1024)
#plot(out)
#dev.off()

