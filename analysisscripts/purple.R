library(RMySQL)

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

query_cancer_type<-function(dbConnect) {
  query = paste(
    "SELECT s.sampleId, c.cancertype",
    " FROM clinical c, sample s",
    "WHERE s.sampleId = c.sampleId",
    sep = " ")
  return (sqlColumnToRownames(dbGetQuery(dbConnect, query), "sampleId")) 
}

query_sv_count<-function(dbConnect) {
  query = paste(
    "SELECT sampleId, count(*) as svCount",
    " FROM structuralVariant",
    "GROUP BY sampleId",
    sep = " ")
  return (sqlColumnToRownames(dbGetQuery(dbConnect, query), "sampleId")) 
}

load_data_from_database<-function() {
  # Query DB
  dbConnect = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
  purities = query_purity(dbConnect)
  cancerTypes = query_cancer_type(dbConnect)
  svCounts = query_sv_count(dbConnect)
  dbDisconnect(dbConnect)
  
  #Join
  samples <- leftJoin(purities, cancerTypes)
  samples <- leftJoin(samples, svCounts)
  return (samples)
}

##### START OF EXECUTION
dataFile = "~/hmf/purple.RData"
fromFile = TRUE
if (fromFile) {
  load(dataFile)
} else {
  samples = load_data_from_database();
  save(samples, file = dataFile)
}

### Clean Data
samplesMissingCancerTypes = samples[is.na(samples$cancerType),]
samplesMissingStructuralVariants = samples[samples$svCount == 0,]
samples = samples[!is.na(samples$cancerType),]


### Experiment
breastSamples = rownames(samples[samples$cancerType == "Breast",])


