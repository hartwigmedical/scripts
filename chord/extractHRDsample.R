#!/usr/bin/Rscript

library('devtools')
args <- commandArgs(TRUE)

load_all(paste(args[1], '/mutSigExtractor', sep = ""))
setwd(args[2])
sample = extractSigsForHrdClassifier(args[4], args[4], args[5], args[3], "gridss")
write.table(sample, file = paste(args[3], "_sample_mut_signatures.txt", sep = ""))