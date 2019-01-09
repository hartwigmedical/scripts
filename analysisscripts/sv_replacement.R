#library(purple)
library(RMySQL)
library(tidyverse)
library(Biostrings)
library(StructuralVariantAnnotation)
library(testthat)
tmpwd=getwd()
setwd("../gridss/")
source("libgridss.R")
setwd(tmpwd)
remove(tmpwd)
source("libSvAnalyser.R")

library(BSgenome.Hsapiens.UCSC.hg19)
library(stringdist)
library(ggplot2)

dbobj = dbConnect(MySQL(), dbname='hmfpatients')
query = paste("SELECT * ",
              "FROM structuralVariant ",
              "WHERE abs(endPosition - startPosition - length(insertSequence)) < 5 AND type='DEL'",
              sep="")
dbdf = dbGetQuery(dbobj, query)
gr = GRanges(seqnames=dbdf$startChromosome, ranges=IRanges(start=dbdf$startPosition, end=dbdf$endPosition), insSeq=dbdf$insertSequence)
seqlevelsStyle(gr) = "UCSC"
gr$refSeq = getSeq(BSgenome.Hsapiens.UCSC.hg19, names=gr, as.character=TRUE)
gr$revSeq = as.character(reverseComplement(DNAStringSet(gr$refSeq)))
gr$fwdEditDistance=stringdist(gr$insSeq, gr$refSeq, method="lv")
gr$invEditDistance=stringdist(gr$insSeq, gr$revSeq, method="lv")

ggplot(as.data.frame(gr)) +
  aes(fwdEditDistance / nchar(insSeq), invEditDistance / nchar(insSeq)) +
  geom_point()

as.data.frame(gr) %>%
  filter(fwdEditDistance / nchar(insSeq) > 0.5 & invEditDistance / nchar(insSeq) < 0.2)
