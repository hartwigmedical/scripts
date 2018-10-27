#!/usr/bin/Rscript

library('devtools')
load_all('/data/common/tools/chord_v1.0/mutSigExtractor')
args <- commandArgs(TRUE)
setwd(args[1])
sample = extractSigsForHrdClassifier(args[3], args[3], args[4], args[2])
write.table(sample, file = paste(args[2], "_sample_mut_signatures.txt", sep = ""))