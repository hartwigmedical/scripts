#!/usr/bin/Rscript

library('devtools')
load_all('/data/experiments/hrd_classifier_evaluation/tools/mutSigExtractor')
args <- commandArgs(TRUE)
setwd(args[4])
sample = extractSigsForHrdClassifier(args[1], args[1], args[2], args[3])
write.table(sample, "sample_mut_signatures.txt")