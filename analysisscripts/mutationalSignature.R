detach("package:purple", unload=TRUE);
library(purple);
library(RMySQL)
library(data.table)
library(MutationalPatterns)
library(ggplot2)


dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")

######### Create cohort #########
sampleCohort = purple::query_purity(dbProd)
clinical = purple::query_clinical_data(dbProd)
sampleCohort = merge(sampleCohort, clinical[, c("sampleId", "cancerType")], all.x = T)
patientIdLookups = query_patient_id_lookup(dbProd)
sampleCohort$patientId <- purple::apply_to_cohort(sampleCohort, function(x) {purple::sample_to_patient_id(x$sampleId, patientIdLookups)})$V1
patientCohort = purple::highest_purity_patients(sampleCohort)
lungCohort = patientCohort[!is.na(patientCohort$cancerType) & patientCohort$cancerType == "Lung", ]
rm(patientIdLookups)
rm(clinical)

####### Querying Variants ########
lungVariants = purple::query_somatic_variants(dbProd, lungCohort)


dbDisconnect(dbProd)
rm(dbProd)


## Do entire cohort!
cohortCosmicSignature = purple::mutational_signature_by_clonality(lungVariants)
plot_cosmic_signature(cohortCosmicSignature)

## do individaul samples
for (sampleId in lungCohort$sampleId) {
  cat("Processing", sampleId, "\n")
  sampleCosmicSignature = mutational_signature_by_clonality(lungVariants[lungVariants$sampleId == sampleId, ])
}

