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
multipleBiopsyCohort = purple::multiple_biopsy(cohort)
multipleBiopsyPatientsId = unique(multipleBiopsyCohort$patientId)
save(cohort, multipleBiopsyCohort, multipleBiopsyPatientsId, file="~/hmf/multipleBiopysyCohort.RData")

## Get somatic variants
multipleBiopsySomaticVariants = purple::query_variant_trinucleotides(dbProd, multipleBiopsyCohort)
save(multipleBiopsySomaticVariants, file="~/hmf/multipleBiopsySomaticVariants.RData")

## Get indels
multipleBiopsyIndels = purple::query_indel_signature(dbProd, multipleBiopsyCohort)
save(multipleBiopsyIndels, file="~/hmf/multipleBiopsyIndels.RData")

## Get structural variants
multipleBiopsyStructuralVariants = purple::query_structural_variants(dbProd, multipleBiopsyCohort)
save(multipleBiopsyStructuralVariants, file="~/hmf/multipleBiopsyStructuralVariants.RData")

### Clean up
dbDisconnect(dbProd)
rm(dbProd)

#### LOAD DATA
load(file="~/hmf/multipleBiopysyCohort.RData")
load(file="~/hmf/multipleBiopsySomaticVariants.RData")
load(file="~/hmf/multipleBiopsyStructuralVariants.RData")
load(file="~/hmf/multipleBiopsyIndels.RData")

### Mutational Signatures
multipleBiopsySomaticSignature = list()
multipleBiopsySVSignature = list()
multipleBiopsyIndelSignature = list()
for (patientId in multipleBiopsyPatientsId) {
  cat("Processing", patientId, "\n")
  
  patientSampleIds = unique(multipleBiopsyCohort[multipleBiopsyCohort$patientId == patientId, c("sampleId")])
  
  patientSomaticVariants = multipleBiopsySomaticVariants[multipleBiopsySomaticVariants$sample %in% patientSampleIds,]
  patientSomaticSignature = purple::signature_matrix_by_scope(patientSomaticVariants)
  multipleBiopsySomaticSignature[[patientId]] <- patientSomaticSignature
  
  patientStructualVariants = multipleBiopsyStructuralVariants[multipleBiopsyStructuralVariants$sample %in% patientSampleIds,]
  patientSVSignature = purple::sv_signature_by_scope(patientStructualVariants)
  multipleBiopsySVSignature[[patientId]] <- patientSVSignature
  
  patientIndels = multipleBiopsyIndels[multipleBiopsyIndels$sampleId %in% patientSampleIds,]
  patientIndelSignature = purple::indel_signature_by_scope(patientIndels)
  multipleBiopsyIndelSignature[[patientId]] <- patientIndelSignature
}

save(multipleBiopsySomaticSignature, multipleBiopsySVSignature, multipleBiopsyIndelSignature, file="~/hmf/multipleBiopsySignatures.RData")