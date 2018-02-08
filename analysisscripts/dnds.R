library(RMySQL)
library(dndscv)
library(IRanges)
detach("package:purple", unload=TRUE); library(purple);

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
cat("Querying purple")
rawCohort = purple::query_purity(dbProd)
cohort = rawCohort

# PatientIds
cat("Mapping samples to patients")
patientIdLookups = query_patient_id_lookup(dbProd)
patientIds = purple::apply_to_cohort(cohort, function(x) {purple::sample_to_patient_id(x$sampleId, patientIdLookups)})
cohort$patientId <- patientIds$V1

#Clinical Data
cat("Querying clinical data")
clinicalData = purple::query_clinical_data(dbProd)
cohort = left_join(cohort, clinicalData[, c("sampleId", "cancerType")])

# Cohort
cohort = purple::highest_purity_patients(cohort)
save(cohort, file="/Users/jon/hmf/RData/dnds/dndsCohort.RData")

# Somatics 
cat("Querying somatics")
rawSomatics = purple::query_somatic_variants(dbProd, cohort, filterEmptyGenes = TRUE)
save(rawSomatics, file="/Users/jon/hmf/RData/dnds/dndsSomatics.RData")
somatics = rawSomatics[rawSomatics$type == "INDEL" | rawSomatics$type == "SNP", c("sampleId", "chromosome", "position", "ref", "alt")]
colnames(somatics) <- c("sampleId", "chr", "pos", "ref", "alt")

# Takes AGES to annotate with cancer type... don't bother!
#somatics$cancerType <- "UNKNOWN"
#for (sampleId in unique(somatics$sampleId)) {
#  matchedCancerType = cohort[match(sampleId, cohort$sampleId), c("cancerType")]
#  cat("SampleId:", sampleId, ", cancerType:", matchedCancerType, "\n")
#  somatics[somatics$sampleId == sampleId, ]$cancerType <- matchedCancerType
#}

# Clean up DB Connection
dbDisconnect(dbProd)
rm(dbProd)
rm(patientIdLookups)
rm(patientIds)
rm(clinicalData)

cancerTypes = unique(cohort$cancerType)
cancerTypes = cancerTypes[!is.na(cancerTypes)]

dndsResults = list()
for (cancerType in cancerTypes[!is.na(cancerTypes)]) {
  cat("Processing", cancerType)
  cancerTypeSampleIds = cohort[!is.na(cohort$cancerType) & cohort$cancerType == cancerType, c("sampleId")]
  input = somatics[somatics$sampleId %in% cancerTypeSampleIds, c("sampleId", "chr", "pos", "ref", "alt")]
  output = dndscv(input)
  dndsResults[[cancerType]] <- output$sel_cv
  save(dndsResults, file="/Users/jon/hmf/RData/dnds/dndsSelCV.RData")
}

output = dndscv(somatics)

output$sel_cv
dndsResults[["All"]] <- output$sel_cv
names(dndsResults)

panCancer <-  output$sel_cv


pcawgRaw = read.csv("/Users/jon/hmf/pcawg/PCAWG_counts.txt", sep = '\t')
