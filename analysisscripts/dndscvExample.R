library(MASS)
library(dndscv)
library(MutationalPatterns)
library(RMySQL)
library(data.table)
library("NMF")

sample_mutations_all_snv<-function(dbConnect, cohort)
{
  query = paste(
    "select sv.sampleId, chromosome as chr, position as pos, ref,alt from somaticVariant sv",
    "where length(alt) = length(ref) and filter = 'PASS' and gene <> '' limit 20000000",
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
    "%' and length(alt) = length(ref) and filter = 'PASS' and gene <> '' limit 10000000",
    sep = "")
  print(query)
  raw_data = dbGetQuery(dbConnect, query)
}

dbConnect = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
mutations<-sample_mutations(dbConnect,'BREAST')
dbDisconnect(dbConnect)
output<-dndscv(mutations)
output

nrow(mutations)
