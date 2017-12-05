library(RMySQL)
library(plyr)
LOAD_FROM_FILE = TRUE
DATA_FILE = "~/hmf/purple.RData"

leftJoin<-function(left, right) {
  jon = merge(x = left, y = right, by="row.names", all.x=TRUE)
  return (sqlColumnToRownames(jon, "Row.names"))
}

sampleIdRowNames<-function(df, .drop = TRUE) {
  row.names(df) <- df$sampleId
  if (.drop) {
    df$sampleId <- NULL
  }
  return (df)
}

applySamples<-function(cohort, sampleFunction) {
  result = ddply(cohort, 1, sampleFunction, .drop = FALSE, .progress = "text", .parallel = FALSE)
  return (sampleIdRowNames(result))
}

####### PURITY #######
query_purity<-function(dbConnect) {
  query = paste(
    "SELECT p.sampleId, p.purity AS purplePurity, round(p.ploidy,2) as ploidy, gender",
    " FROM purity p",
    "WHERE qcStatus = 'PASS'",
    "AND status <> 'NO_TUMOR'",
    sep = " ")
  return (sampleIdRowNames(dbGetQuery(dbConnect, query), .drop = FALSE))
}

####### Sample Data #######
query_sample_data<-function(dbConnect) {
  query = paste(
    "SELECT sampleId, tumorPercentage AS pathologyPurity",
    " FROM sample c",
    sep = " ")
  return (sampleIdRowNames(dbGetQuery(dbConnect, query)))
}

attach_sample_data<-function(dbConnect, samples) {
  cancerTypes = query_sample_data(dbConnect)
  return (leftJoin(samples, cancerTypes))
}

####### Clinical Data #######
query_clinical_data<-function(dbConnect) {
  query = paste(
    "SELECT c.sampleId, c.cancertype",
    " FROM clinical c",
    sep = " ")
  return (sampleIdRowNames(dbGetQuery(dbConnect, query)))
}

attach_clinical_data<-function(dbConnect, samples) {
  cancerTypes = query_clinical_data(dbConnect)
  return (leftJoin(samples, cancerTypes))
}

####### Microsatellite Instability #######
query_msi_sample<-function(dbConnect, sample) {
  query = paste(
    "SELECT count(*) as msiScore ",
    " FROM somaticVariant ",
    "WHERE filter = 'PASS'",
    "   AND length(alt) <> length(ref)",
    "   AND length(alt) <= 50",
    "   AND length(ref) <= 50",
    "   AND sampleId = '", sample$sampleId, "'",
    " GROUP BY sampleId",
    sep = "")
  score = dbGetQuery(dbConnect, query)
  return (score)
}

attach_msi<-function(dbConnect, cohort) {
  cat("Attaching microsatellite instability...")
  score = applySamples(cohort, function(x) {query_msi_sample(dbConnect, x)})
  score$msiStatus <-ifelse(score$msiScore > 3095, "UNSTABLE", "STABLE")
  cat("Attaching microsatellite instability complete!")
  return (leftJoin(cohort, score))
}

####### Structural Variants #######
query_structural_variants<-function(dbConnect, sample) {
  query = paste(
    "SELECT sampleId, ",
    "       SUM(1) AS svTotal, ",
    "       SUM(IF (type ='DEL', 1, 0)) AS svDel, ",
    "       SUM(IF (type ='DUP', 1, 0)) AS svDup, ",
    "       SUM(IF (type ='INS', 1, 0)) AS svIns, ",
    "       SUM(IF (type ='INV', 1, 0)) AS svInv, ",
    "       SUM(IF (type ='BND', 1, 0)) AS svBnd ",
    " FROM structuralVariant ",
    "WHERE sampleId = '", sample$sampleId, "'",
    " GROUP BY sampleId",
    sep = "")
  return (dbGetQuery(dbConnect, query))
}

attach_structural_variants<-function(dbConnect, cohort) {
  cat("Attaching somatics....")
  somatics = applySamples(cohort, function(x) {query_structural_variants(dbConnect, x)})
  cat("Attaching somatics complete!")
  return (leftJoin(cohort, somatics))
}

####### QC Score #######
query_qc_score<-function(dbConnect, sample) {
  query = paste(
    "SELECT sampleId, count(*) as unsupportedSegments ",
    " FROM copyNumber ",
    "WHERE sampleId = '", sample$sampleId, "'",
    "  AND segmentStartSupport = 'NONE'",
    " GROUP BY sampleId",
    sep = "")
  return (dbGetQuery(dbConnect, query))
}

attach_qc_score<-function(dbConnect, cohort) {
  cat("Attaching qc score")
  qcScore = applySamples(cohort, function(x) {query_qc_score(dbConnect, x)})
  result = leftJoin(cohort, qcScore) 
  result$qcScore <- round(result$unsupportedSegments / result$ploidy,0)
  result$unsupportedSegments <- NULL
  cat("Attaching qc score!")
  return (result)
}

####### Somatics #######
query_somatic_sample<-function(dbConnect, sample) {
  query = paste(
    "SELECT sampleId, ",
    "       SUM(1) AS somaticVariants, ",
    "       SUM(IF (clonality ='CLONAL', 1, 0)) AS clonal, ",
    "       SUM(IF (clonality ='SUBCLONAL', 1, 0)) AS subclonal, ",
    "       SUM(IF (clonality ='INCONSISTENT', 1, 0)) AS inconsistentClonality, ",
    "       SUM(IF (length(ref) = 1 AND length(alt) = 1, 1, 0)) as snv, ",
    "       SUM(IF (length(ref) <> length(alt), 1, 0)) as indels, ",
    "       SUM(IF (effect like '%missense%', 1, 0)) AS mutationalLoad ",
    " FROM somaticVariant ",
    "WHERE filter = 'PASS' ",
    "   AND sampleId = '", sample$sampleId, "'",
    " GROUP BY sampleId",
    sep = "")
  return (dbGetQuery(dbConnect, query))
}

attach_somatics<-function(dbConnect, cohort) {
  cat("Attaching somatics....")
  somatics = applySamples(cohort, function(x) {query_somatic_sample(dbConnect, x)})
  cat("Attaching somatics complete!")
  return (leftJoin(cohort, somatics))
}

####### START OF EXECUTION #######
if (LOAD_FROM_FILE) {
  load(DATA_FILE)
} else {
  
  #Pilot
  dbConnect = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
  allSamples = query_purity(dbConnect)
  
  cohort = allSamples
  #cohort = allSamples[1:7, , drop = FALSE]
  cohort = attach_qc_score(dbConnect, cohort)
  cohort = attach_sample_data(dbConnect, cohort)
  cohort = attach_structural_variants(dbConnect, cohort)
  cohort = attach_msi(dbConnect, cohort)
  cohort = attach_somatics(dbConnect, cohort)
  dbDisconnect(dbConnect)
  
  #Production
  prodDB = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
  cohort = attach_clinical_data(prodDB, cohort)
  dbDisconnect(prodDB)
  
  save(allSamples, cohort, file = DATA_FILE)
  dbDisconnect(dbConnect)
}

### Missing Data
#samplesMissingCancerTypes = samples[is.na(samples$cancerType),]
#samplesMissingStructuralVariants = samples[samples$svCount == 0,]

### Experiment
#breastSamples = rownames(samples[samples$cancerType == "Breast",])

jon = allSamples[1:5, 1, drop = FALSE]
