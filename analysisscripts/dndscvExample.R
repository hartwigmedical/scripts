### ssh -L 3307:localhost:3306 peter@ext-hmf-datastore
library(MASS)
library(dndscv)
library(MutationalPatterns)
library(RMySQL)
library(data.table)
library("NMF")
library(devtools)

sample_mutations_all_snv<-function(dbConnect)
{
  query = paste(
    "select s.sampleId, chromosome as chr, position as pos, ref,alt  from somaticVariant s, purity p where s.sampleId = p.sampleId ",
    "and length(alt) = length(ref) and filter = 'PASS' and gene <> '' ",
    "and status <> 'NO_TUMOR' and qcstatus = 'PASS' ",
    "and right(s.sampleId,1) ='T'",
    sep = "")
  print(query)
  raw_data = dbGetQuery(dbConnect, query)
}

sample_mutations<-function(dbConnect, cohort)
{
  query = paste(
    "select sv.sampleId, chromosome as chr, position as pos, ref,alt from somaticVariant sv,sample s,patient p ",
    "where sv.sampleId = s.sampleId and s.patientId = p.id and upper(primaryTumorLocation) like '%",
    cohort,
    "%' and length(alt) = length(ref) and filter = 'PASS' and gene <> '' limit 20000000 and right(s.sampleId,1) ='T'",
    sep = "")
  print(query)
  raw_data = dbGetQuery(dbConnect, query)
}

DNDSCV_DATA_FILE="~/hmf/dndscv.RData"

####### LOAD FROM DB ################
dbConnect = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
#mutations<-sample_mutations(dbConnect,'BREAST')
### TO DO - ENSURE WE ONLY GET 1 copy of each tumor!!!!
mutations<-sample_mutations_all_snv(dbConnect)
dbDisconnect(dbConnect)
save(mutations, cohort,  file = DNDSCV_DATA_FILE)

###########  LOAD FROM FILE

load(DNDSCV_DATA_FILE)

########### RUN dndsCV #############

output<-dndscv(mutations)
output

nrow(mutations)
