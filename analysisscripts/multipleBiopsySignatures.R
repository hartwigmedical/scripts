detach("package:purple", unload=TRUE)
library(purple)
library(RMySQL)
library(data.table)
library(IRanges)

load(file="~/hmf/cohort.RData")
multipleBiopsyCohort = purple::multiple_biopsy(cohort)
multipleBiopsyPatientsId = unique(multipleBiopsyCohort$patientId)
rm(cohort)

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")


multipleBiopysySignatures = list()
for (patientId in multipleBiopsyPatientsId) {
  cat("Processing", patientId, "\n")
  
  samples = multipleBiopsyCohort[multipleBiopsyCohort$patientId == patientId, ]
  patientVariants = purple::query_variant_trinucleotides(dbProd, samples)
  patientSignature = purple::signature_matrix_by_scope(patientVariants)
  
  multipleBiopysySignatures[[patientId]] <- patientSignature
}

save(multipleBiopysySignatures, file="~hmf/multipleBiopysySignatures.RData")

### Clean up
dbDisconnect(dbProd)
rm(dbProd)

