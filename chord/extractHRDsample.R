#!/usr/bin/Rscript

library('devtools')
args <- commandArgs(TRUE)

load_all(args[1])
setwd(args[2])
sample = extractSigsForHrdClassifier(args[4], args[4], args[5], args[3])
write.table(sample, file = paste(args[3], "_sample_mut_signatures.txt", sep = ""))