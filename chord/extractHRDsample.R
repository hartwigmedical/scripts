#!/usr/bin/Rscript

library('devtools')
load_all('/data/common/tools/chord_v0.1/mutSigExtractor')
args <- commandArgs(TRUE)
setwd(args[4])
sample = extractSigsForHrdClassifier(args[1], args[1], args[2], args[3])
write.table(sample, "_sample_mut_signatures.txt", sep = "")