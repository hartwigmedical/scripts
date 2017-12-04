library(RMySQL)
library(plyr)
LOG_QUERIES = FALSE
LOAD_FROM_FILE = FALSE
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
    "SELECT p.sampleId, p.purity",
    " FROM purity p",
    "WHERE qcStatus = 'PASS'",
    "AND status <> 'NO_TUMOR'",
    sep = " ")
  return (sampleIdRowNames(dbGetQuery(dbConnect, query), .drop = FALSE))
}

####### Cancer Types #######
query_cancer_type<-function(dbConnect) {
  query = paste(
    "SELECT s.sampleId, c.cancertype",
    " FROM clinical c, sample s",
    "WHERE s.sampleId = c.sampleId",
    sep = " ")
  return (sampleIdRowNames(dbGetQuery(dbConnect, query)))
}

attach_cancerTypes<-function(dbConnect, samples) {
  cancerTypes = query_cancer_type(dbConnect)
  return (leftJoin(samples, cancerTypes))
}


####### Structual Variants #######
query_sv_count<-function(dbConnect) {
  query = paste(
    "SELECT sampleId, count(*) as svCount",
    " FROM structuralVariant",
    "GROUP BY sampleId",
    sep = " ")
  return (sampleIdRowNames(dbGetQuery(dbConnect, query)))
}

attach_structural_variants<-function(dbConnect, samples) {
    svCounts = query_sv_count(dbConnect)
    return (leftJoin(samples, svCounts))
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
  if (LOG_QUERIES) {
    cat(paste(query, "\n", sep = ""))
  }

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

####### Somatics #######
query_somatic_sample<-function(dbConnect, sample) {
  query = paste(
    "SELECT sampleId, sum(if (clonality ='CLONAL', 1, 0)) as clonal, sum(if (clonality ='SUBCLONAL', 1, 0)) as subclonal",
    " FROM somaticVariant ",
    "WHERE filter = 'PASS' ",
    "   AND sampleId = '", sample$sampleId, "'",
    " GROUP BY sampleId",
    sep = "")
  if (LOG_QUERIES) {
    cat(paste(query, "\n", sep = ""))
  }
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
  dbConnect = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
  allSamples = query_purity(dbConnect)

  samples = allSamples[1:7, , drop = FALSE]
  samples = attach_cancerTypes(dbConnect, samples)
  samples = attach_structural_variants(dbConnect, samples)
  samples = attach_msi(dbConnect, samples)
  samples = attach_somatics(dbConnect, samples)

  save(allSamples, samples, file = DATA_FILE)
  dbDisconnect(dbConnect)
}

### Missing Data
#samplesMissingCancerTypes = samples[is.na(samples$cancerType),]
#samplesMissingStructuralVariants = samples[samples$svCount == 0,]

### Experiment
#breastSamples = rownames(samples[samples$cancerType == "Breast",])

jon = allSamples[1:5, 1, drop = FALSE]
