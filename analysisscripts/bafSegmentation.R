# Parse the arguments
args <- commandArgs(trailing=T)
bafFile <- args[1]
column <- args[2]
pcfFile   <- args[3]
pcfFile1  <- paste(args[3], "1", sep="")

library(copynumber)
baf <- read.table(bafFile, header=TRUE)
baf <- baf[,c("Chromosome","Position",column)]
baf.seg<-pcf(baf,verbose=FALSE,gamma=100,kmin=1,save.res = TRUE, file.names = c(pcfFile1, pcfFile))