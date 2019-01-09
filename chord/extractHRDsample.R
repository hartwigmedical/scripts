#!/usr/bin/Rscript

library('devtools')
args <- commandArgs(TRUE)

load_all(paste(args[1], '/mutSigExtractor', sep = ""))
setwd(args[2])

cat("[INFP] script extractHRDsample", "\n")
cat("[INFO] chord_dir: ", args[1], "\n")
cat("[INFO] working_dir: ", args[2], "\n")
cat("[INFO] sample: ", args[3], "\n")
cat("[INFO] somatic vcf: ", args[4], "\n")
cat("[INFO] structural vcf: ", args[5], "\n")

sample = extractSigsForHrdClassifier(args[4], args[4], args[5], args[3], "gridss")
write.table(sample, file = paste(args[3], "_sample_mut_signatures.txt", sep = ""))