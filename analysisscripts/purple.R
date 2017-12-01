library(RMySQL)
LOG_QUERIES = TRUE
LOAD_FROM_FILE = TRUE
DATA_FILE = "~/hmf/purple.RData"

leftJoin<-function(left, right) {
  jon = merge(x = left, y = right, by="row.names", all.x=TRUE)
  return (sqlColumnToRownames(jon, "Row.names"))
}

query_purity<-function(dbConnect) {
  query = paste(
    "SELECT p.sampleId, p.purity",
    " FROM purity p",
    "WHERE qcStatus = 'PASS'",
    "AND status <> 'NO_TUMOR'",
    sep = " ")
  return (sqlColumnToRownames(dbGetQuery(dbConnect, query), "sampleId"))
}

####### Cancer Types #######
attach_cancerTypes<-function(dbConnect, samples) {
  cancerTypes = query_cancer_type(dbConnect)
  return (leftJoin(samples, cancerTypes))
}

query_cancer_type<-function(dbConnect) {
  query = paste(
    "SELECT s.sampleId, c.cancertype",
    " FROM clinical c, sample s",
    "WHERE s.sampleId = c.sampleId",
    sep = " ")
  return (sqlColumnToRownames(dbGetQuery(dbConnect, query), "sampleId"))
}

####### Structual Variants #######
attach_structural_variants<-function(dbConnect, samples) {
  svCounts = query_sv_count(dbConnect)
  return (leftJoin(samples, svCounts))
}

query_sv_count<-function(dbConnect) {
  query = paste(
    "SELECT sampleId, count(*) as svCount",
    " FROM structuralVariant",
    "GROUP BY sampleId",
    sep = " ")
  return (sqlColumnToRownames(dbGetQuery(dbConnect, query), "sampleId"))
}

####### Microsatellite Instability #######
attach_msi<-function(dbConnect, samples) {
  samples$msiScore <- query_msi_cohort(dbConnect, samples)
  samples$msiStatus <- ifelse(samples$msiScore > 3095, "UNSTABLE", "STABLE")
  samples$msiStatus <- ifelse(samples$msiScore == 0, NA, samples$msiStatus)
  return (samples)
}

query_msi_cohort<-function(dbConnect, cohort) {
  return (sapply(rownames(cohort),function(x) {query_msi_sample(dbConnect, x)}))
}

query_msi_sample<-function(dbConnect, sample) {
  query = paste(
    "SELECT count(*) as msiScore ",
    " FROM somaticVariant ",
    "WHERE filter = 'PASS'",
    "   AND length(alt) <> length(ref)",
    "   AND length(alt) <= 50",
    "   AND length(ref) <= 50",
    "   AND sampleId = '", sample, "'",
    " GROUP BY sampleId",
    sep = "")
  if (LOG_QUERIES) {
    cat(paste(query, "\n", sep = ""))
  }
  
  score = dbGetQuery(dbConnect, query)
  return (score)
}


####### START OF EXECUTION #######
if (LOAD_FROM_FILE) {
  load(DATA_FILE)
} else {
  dbConnect = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
  allSamples = query_purity(dbConnect)
  allSamples = attach_cancerTypes(dbConnect, allSamples)
  
  samples = allSamples[1:5, ]
  samples = attach_structural_variants(dbConnect, samples)
  samples = attach_msi(dbConnect, samples)

  save(allSamples, samples, file = DATA_FILE)
  dbDisconnect(dbConnect)
}

### Missing Data
#samplesMissingCancerTypes = samples[is.na(samples$cancerType),]
#samplesMissingStructuralVariants = samples[samples$svCount == 0,]

### Experiment
#breastSamples = rownames(samples[samples$cancerType == "Breast",])

