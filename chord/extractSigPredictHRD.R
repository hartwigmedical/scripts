#!/usr/bin/Rscript
args <- commandArgs(TRUE)

chordToolDir <- args[1]
  workingDir <- args[2]
  sampleName <- args[3]
   snvIndVcf <- args[4]
       svVcf <- args[5]
   sigOutTxt <- paste0( workingDir, '/', sampleName, '_chord_signatures.txt')
   prdOutTxt <- paste0( workingDir, '/', sampleName, '_chord_prediction.txt')

cat("[INFO] START CHORD signature extraction and HRD prediction", "\n")
setwd(workingDir)

suppressPackageStartupMessages(library('devtools'))
suppressPackageStartupMessages(library('randomForest'))
suppressPackageStartupMessages(load_all(paste0(chordToolDir, '/mutSigExtractor-1.03')))
suppressPackageStartupMessages(load_all(paste0(chordToolDir, '/CHORD-60.02')))

cat("[INFO] Settings:\n")
cat("[INFO]   Chord dir:", chordToolDir, "\n")
cat("[INFO]   Working dir:", workingDir, "\n")
cat("[INFO]   Sample name:", sampleName, "\n")
cat("[INFO]   Somatic SNV/IND vcf:", snvIndVcf, "\n")
cat("[INFO]   Somatic SV vcf:", svVcf, "\n")
cat("[INFO]   Signature out file:", sigOutTxt, "\n")
cat("[INFO]   Prediction out file:", prdOutTxt, "\n")

## Extract signature counts and predict HRD
cat("[INFO] Performing chord signature extraction\n")
signatures <- extractSigsChord(
  vcf.snv = snvIndVcf,
  vcf.indel = snvIndVcf,
  vcf.sv = svVcf,
  sample.name = sampleName,
  sv.caller = "gridss"
)
cat("[INFO] Performing chord HRD prediction\n")
prediction <- chordPredict(
  signatures, 
  rf.model = CHORD, 
  hrd.cutoff = 0.5
)

## Output
write.table(signatures, file=sigOutTxt, sep="\t")
write.table(prediction, file=prdOutTxt, sep="\t", quote=FALSE, row.names=FALSE)

cat("[INFO] FINISHED CHORD signature extraction and HRD prediction\n")

