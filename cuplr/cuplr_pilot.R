#!/usr/bin/Rscript
options(stringsAsFactors=F) # to avoid invalid factor level warning

args <- commandArgs(TRUE)

   toolDir <- args[1]
workingDir <- args[2]
sampleName <- args[3]
    somVcf <- args[4] # .purple.somatic.vcf.gz
     svVcf <- args[5] # .purple.sv.vcf.gz
 purplePur <- args[6] # .purple.purity.tsv
 purpleCnv <- args[7] # .purple.cnv.somatic.tsv
   linxFus <- args[8] # .linx.fusion.tsv
   linxIns <- args[9] # .linx.viral_inserts.tsv
   linxVis <- args[10] # .linx.vis_sv_data.tsv
   germVcf <- args[11] # .germline.vcf.gz
 genomeVsn <- args[12] # HG19 or HG38
featOutTxt <- paste0(workingDir,'/', sampleName,'_cuplr_features.txt')
predOutTxt <- paste0(workingDir,'/', sampleName,'_cuplr_prediction.txt')
highOutTxt <- paste0(workingDir,'/', sampleName,'_cuplr_highest.csv')

inputFiles <- c(somVcf, svVcf, purplePur, purpleCnv, linxFus, linxIns, linxVis, germVcf)
outputFiles <- c(featOutTxt, predOutTxt, highOutTxt)

cat("[INFO] Starting CUPLR for", sampleName, "\n")
setwd(workingDir)

suppressPackageStartupMessages(library('devtools'))
suppressPackageStartupMessages(library('randomForest'))
suppressPackageStartupMessages(library('rfFC'))
suppressPackageStartupMessages(library("seqminer"))
suppressPackageStartupMessages(library("GenomeInfoDb"))

suppressPackageStartupMessages(load_all(paste0(toolDir, '/mutSigExtractor')))
suppressPackageStartupMessages(load_all(paste0(toolDir, '/cuplr/commonUtils')))
suppressPackageStartupMessages(load_all(paste0(toolDir, '/cuplr/featureExtractor')))
suppressPackageStartupMessages(load_all(paste0(toolDir, '/cuplr/geneDriverAnnotator')))
suppressPackageStartupMessages(load_all(paste0(toolDir, '/cuplr/cuplr')))

cat("[INFO] Package NamespaceVersions:\n")
for (pkgName in c("mutSigExtractor", "cuplr")){
  pkgVsn <- getNamespaceVersion(pkgName)[["version"]]
  cat("[INFO]  ", pkgName, "has version", pkgVsn, "\n")
}

## Convert genome name to BSGenome name
if (genomeVsn == "HG19") {
  suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
  refGenome <- BSgenome.Hsapiens.UCSC.hg19
} else if (genomeVsn == "HG38") {
  suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  refGenome <- BSgenome.Hsapiens.UCSC.hg38
} else {
  stop("Unsupported ref genome version: ", genomeVsn," (should be HG19 or HG38)\n")
}

cat("[INFO] Settings:\n")
cat("[INFO]   ToolDir =", toolDir, "\n")
cat("[INFO]   WorkingDir =", workingDir, "\n")
cat("[INFO]   SampleName =", sampleName, "\n")
cat("[INFO]   GenomeVersion =", genomeVsn, "\n")
cat("[INFO] Input Files:\n")
for (inputFile in inputFiles){
  cat("[INFO]  ", basename(inputFile), "\n")
}
cat("[INFO] Output Files:\n")
for (outputFile in outputFiles){
  cat("[INFO]  ", basename(outputFile), "\n")
}

cat("[INFO] Performing signature extraction...\n")
features <- extractFeaturesCuplr(
   germ.vcf.path=germVcf,
   som.vcf.path=somVcf,
   sv.vcf.path=svVcf,
   purple.cnv.path=purpleCnv,
   purple.purity.path=purplePur,
   linx.fusion.path=linxFus,
   linx.viral.inserts.path=linxIns,
   linx.vis.sv.data.path=linxVis,
   out.dir=workingDir, 
   sample.name=sampleName,
   write.features=FALSE, # with TRUE would write features.txt.gz to disk
   return.features=TRUE, # needed because we do not write features.txt.gz to disk
   verbose=2, # level 2 verbosity passes through errors from external functions
   colname.translations=list( # translations needed because column names changed
     purple.cnv=c(
       chrom='chromosome',
       start='start',
       end='end',
       total_cn='copyNumber',
       major_cn='majorAlleleCopyNumber',
       minor_cn='minorAlleleCopyNumber'
     )
   )
)

cat("[INFO] Reading features...\n")
features <- readFeaturesCuplr(features)

cat("[INFO] Loading random forest model...\n")
cuplr <- readRDS(CUPLR)

cat("[INFO] Performing prediction...\n")
prediction <- predict.randomForestEnsemble(cuplr, features, type='prob')

highestScore <- prediction[apply(prediction,1,which.max)]
highestType <- colnames(prediction)[apply(prediction,1,which.max)]

## Output
cat("[INFO] Writing output files\n")
write.table(features, file=featOutTxt, sep="\t")
write.table(prediction, file=predOutTxt, sep="\t", quote=FALSE, row.names=FALSE)
writeLines(paste0(highestType, ",", highestScore), highOutTxt)
cat("[INFO] Output in:", workingDir, "\n")

cat("[INFO] Finished CUPLR for", sampleName, "\n")

