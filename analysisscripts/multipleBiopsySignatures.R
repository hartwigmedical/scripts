detach("package:purple", unload=TRUE)
library(purple)
library(RMySQL)
library(data.table)
library(IRanges)

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")

## Get cohort
cohort = purple::query_purity(dbProd)
patientIdLookups = query_patient_id_lookup(dbProd)
patientIds = purple::apply_to_cohort(cohort, function(x) {purple::sample_to_patient_id(x$sampleId, patientIdLookups)})
cohort$patientId <- patientIds$V1
save(cohort, file="~/hmf/multipleBiopysyCohort.RData")

## Get multiple biopsies
multipleBiopsyCohort = purple::multiple_biopsy(cohort)
multipleBiopsyPatientsId = unique(multipleBiopsyCohort$patientId)

## Get somatic variants
allVariants = purple::query_variant_trinucleotides(dbProd, multipleBiopsyCohort)
save(cohort, allVariants, file="~/hmf/multipleBiopysyVariants.RData")

### Mutational Signatures
multipleBiopysySignatures = list()
for (patientId in multipleBiopsyPatientsId) {
  cat("Processing", patientId, "\n")
  
  patientSampleIds = unique(multipleBiopsyCohort[multipleBiopsyCohort$patientId == patientId, c("sampleId")])
  patientVariants = allVariants[allVariants$sample %in% patientSampleIds,]
  patientSignature = purple::signature_matrix_by_scope(patientVariants)
  
  multipleBiopysySignatures[[patientId]] <- patientSignature
}
save(multipleBiopysySignatures, file="~/hmf/multipleBiopysySignatures.RData")




### Clean up
dbDisconnect(dbProd)
rm(dbProd)

#### TEMP



