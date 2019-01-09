#!/usr/bin/Rscript

library('devtools')
args <- commandArgs(TRUE)

cat(args[1] , "/n")
cat(args[2] , "/n")
cat(args[3] , "/n")
cat(args[4] , "/n")
cat(args[5] , "/n")

load_all(paste(args[1], '/mutSigExtractor', sep = ""))
setwd(args[2])
sample = extractSigsForHrdClassifier(args[4], args[4], args[5], args[3], "gridss")
write.table(sample, file = paste(args[3], "_sample_mut_signatures.txt", sep = ""))