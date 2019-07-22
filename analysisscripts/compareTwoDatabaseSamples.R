#!/usr/bin/Rscript
library('getopt')

SQLGROUP = 'RAnalysis'
  LABEL1 = 'TEST'
  LABEL2 = 'PROD'

spec = matrix(c(
  'help',       'h', 0, "logical",
  'configCsv',  'c', 2, "character",
  'outputPath', 'o', 2, "character",
  'outputName', 'n', 2, "character",
  'dbName1',    't', 1, "character",
  'dbName2',    'p', 2, "character",
  'sample1',    'y', 1, "character",
  'sample2',    'z', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

if ( length(commandArgs(TRUE)) < 1 || !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  cat( "  example1: --dbName1 pilot --sample1 CPCTXXXXXXXX\n" )
  cat( "  example2: --dbName1 reference_validation_sets --sample1 COLO829v001T --dbName2 reference_validation_sets --sample2 COLO829v002T\n" )
  cat( "  example3: --configCsv <config.csv>\n" )
  cat( "    sample1,dbName1,sample2,dbName2,tag\n" )
  cat( "    COLO829v001T,reference_validation_sets,COLO829v002T,reference_validation_sets,COLO_v1_vs_v2\n" )
  q(status=1)
}

if ( is.null(opt$outputPath) ){ 
    opt$outputPath = '.' 
}
if ( is.null(opt$configCsv) ){
  for (reqParam in c( 'dbName1', 'sample1' )){
      if ( is.null(opt[[ reqParam ]]) ) {
          cat( "[EXIT] Missing required parameter:", reqParam, " (-h for help)\n" )
          q(status=1)
      }
  }
}

## function: Parses config for multisample comparison
parseConfig <- function(configCsv){
    cat( paste0("[INFO] Parsing config file",configCsv,"\n"))
    comparisonTable <- read.csv(configCsv, stringsAsFactors=FALSE, comment.char = "#")
    return(comparisonTable)
}

## function: Gathers Somatic and Structural variants for one sample
getSampleDataFromDb <- function(dbName, sampleId, label){
    cat( paste0("[INFO]   Querying database ", dbName, " for sample ", sampleId," (label=",label,")\n"))
    dbConn = dbConnect(MySQL(), dbname = dbName, groups = SQLGROUP)
    cohort = data.frame(sampleId=sampleId)
    somVar = purple::query_somatic_variants(dbConn, cohort) %>% mutate(source = label) 
    somVar = somVar %>% mutate(patientId = sampleId, sampleId = paste0(sampleId, "-", source))
     svVar = purple::query_structural_variants(dbConn, cohort) %>% mutate(source = label)
     svVar = svVar %>% mutate(patientId = sampleId, sampleId = paste0(sampleId, "-", source))
    output = list("somaticVariants" = somVar, "structuralVariants" = svVar)
    respon = dbDisconnect(dbConn)
    return(output)
}

## function: Takes variants from two samples and compares the two
compareTwoSamples <- function(testData, prodData, outName, pdfPath){
    #cat( paste0("[INFO] Combining (", dbName, ") for sample (", sampleId,")\n"))
    somVarTest = testData$somaticVariants
    strVarTest = testData$structuralVariants
    somVarProd = prodData$somaticVariants
    strVarProd = prodData$structuralVariants
    
    ## Combine TEST and PROD
    cat(paste0("[INFO]   Combining data and processing ", outName, "\n"))
    somVars = bind_rows(somVarProd, somVarTest)
    strVars = bind_rows(strVarProd, strVarTest)

    ## Extract patient data
    sampleIds = somVars %>% distinct(sampleId) %>% pull()
    sample1 = sampleIds[1]
    patientSomVars = scope(somVars[somVars$sample %in% sampleIds,])
    patientStrVars = scope(strVars[strVars$sample %in% sampleIds,])

    ## Signatures
    patientMutSig = purple::mutational_signature_by_scope(patientSomVars)
    patientIndSig = purple::indel_signature_by_scope(patientSomVars)
    patientStrSig = purple::sv_signature_by_scope(patientStrVars)
    mutSig = list(); indSig = list(); strSig = list()
    mutSig[[sample1]] <- patientMutSig
    indSig[[sample1]] <- patientIndSig
    strSig[[sample1]] <- patientStrSig

    ## Plotting
    cat(paste0("[INFO]   Preparing plots\n"))
    somPloidyPlots = somatic_ploidy_plots(patientSomVars)
    strPloidyPlots = structural_ploidy_plots(patientStrVars)
    p7 <- plot_indel_signature(patientIndSig)
    p8 <- plot_sv_signature(patientStrSig)
    p9 = p8 # Note: if identival samples the following line won't work so we dummy p9 with p8
    p9 <- plot_cosmic_signature(patientMutSig)

    cat(paste0("[INFO]   Creating plots\n"))
    pdf(file = pdfPath, height = 14, width = 20)
    multiplot(
      somPloidyPlots[[1]] + ggtitle(sampleIds[1]), 
      somPloidyPlots[[2]] + ggtitle(sampleIds[2]), 
      somPloidyPlots[[3]],
      strPloidyPlots[[1]] + ggtitle(sampleIds[1]), 
      strPloidyPlots[[2]] + ggtitle(sampleIds[2]), 
      strPloidyPlots[[3]], 
      p7, p8, p9,
      cols = 3
    )
    dev.off()
    cat(paste0("[INFO]   Output written (", pdfPath, ")\n"))
}

## Load requirements
cat( paste0("[INFO] Loading all required packages\n"))
suppressPackageStartupMessages(library(RMySQL))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(MutationalPatterns))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
options(warn=-1)
  suppressPackageStartupMessages(library(purple))
options(warn=0)

## Set and reset variables based on potential user input
outputPath = opt$outputPath
 configCsv = opt$configCsv

if ( !is.null(configCsv) ){ 
  comparisonTable = parseConfig(configCsv) 
} else { 
  dbName1 = opt$dbName1
  sample1 = opt$sample1
  dbName2 = opt$dbName2
  sample2 = opt$sample2
  if ( is.null(dbName2) ){ dbName2 = 'hmfpatients' }
  if ( is.null(sample2) ){ sample2 = sample1 }
  comparisonTable = data.frame( 
    'sample1' = sample1, 'dbName1' = dbName1, 
    'sample2' = sample2, 'dbName2' = dbName2,
    'tag' = paste0(sample1,'-',dbName1,'_vs_',sample2,'-',dbName2),
    stringsAsFactors=FALSE
  ) 
}

## Process all combinations
for (i in 1:nrow(comparisonTable)){
  pairTag = comparisonTable[ i, 'tag' ]
  cat( paste0("[INFO] At comparison ", i, ": ", pairTag,"\n"))
  sample1 = comparisonTable[ i, 'sample1' ]
  sample2 = comparisonTable[ i, 'sample2' ]
  dbName1 = comparisonTable[ i, 'dbName1' ]
  dbName2 = comparisonTable[ i, 'dbName2' ]
  pairTag = comparisonTable[ i, 'tag' ]
    data1 = getSampleDataFromDb(dbName = dbName1, sampleId = sample1, label = LABEL1)
    data2 = getSampleDataFromDb(dbName = dbName2, sampleId = sample2, label = LABEL2)
  pdfPath = paste0(outputPath, '/sampleComparison_', pairTag, '.pdf')
  
  compareTwoSamples(data1, data2, pairTag, pdfPath)
}
