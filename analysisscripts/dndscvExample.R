library(MASS)
library(dndscv)

sample_mutations<-function(dbConnect, cohort)
{
  query = paste(
    "select sv.sampleId, chromosome as chr, position as pos, ref,alt from somaticVariant sv,sample s,patient p ",
    "where sv.sampleId = s.sampleId and s.patientId = p.id and p.primaryTumorLocation='",
    cohort,
    "' and length(alt) = length(ref) and filter = 'PASS' limit 500000",
    sep = "")
  
  raw_data = dbGetQuery(dbConnect, query)
}

dbConnect = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
mutations<-sample_mutations(dbConnect,'Melanoma')
dbDisconnect(dbConnect)
dndscv(mutations)