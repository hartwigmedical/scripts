detach("package:purple", unload=TRUE)
library(purple)
library(RMySQL)

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")

# Query Cohort Data
purity = query_purity(dbProd)
sample = query_sample_data(dbProd)
clinical = query_clinical_data(dbProd)
clinical$ageAtBiopsyDate = year(clinical$biopsyDate) - clinical$birthYear
cohort = purity[, c("sampleId","gender", "purity", "ploidy")]

cat("Querying microsatellite instability", "\n")
msi = purple::apply_to_cohort(cohort, function(x) {purple::query_msi_sample(dbProd, x$sampleId)})

cat("Querying structural variants", "\n")
svTypes = purple::apply_to_cohort(cohort, function(x) {purple::query_structural_variant_overview(dbProd, x$sampleId)})

cat("Querying unsupported copy number segments", "\n")
unsupportedSegments = purple::apply_to_cohort(cohort, function(x) {purple::query_unsupported_segments(dbProd, x$sampleId)})

cat("Querying somatic variants", "\n")
somatics = purple::apply_to_cohort(cohort, function(x) {purple::query_somatic_overview(dbProd, x$sampleId)})

cat("Querying patientIds", "\n")
patientIdLookups = query_patient_id_lookup(dbProd)
patientIds = purple::apply_to_cohort(cohort, function(x) {purple::sample_to_patient_id(x$sampleId, patientIdLookups)})
colnames(patientIds) <- c("sampleId", "patientId")

cat("Querying duplicated autosomes", "\n")
wholeGenomeDuplications = purple::apply_to_cohort(cohort, function(x) {purple::query_whole_genome_duplication(dbProd, x$sampleId)})

# Save
save(clinical, purity, sample, msi, svTypes, unsupportedSegments, somatics, patientIds, wholeGenomeDuplications, file="~/hmf/purpleCohortQueries.RData")

# Load
load("~/hmf/purpleCohortQueries.RData")

# Attached queried data
cohort = purity[purity$sampleId != "CPCT02010063T", c("sampleId","gender", "version", "purity", "ploidy")]
cohort = left_join(cohort, clinical[, c("sampleId", "cancerType", "cancerSubtype", "biopsySite", "biopsyLocation", "treatment", "treatmentType", "ageAtBiopsyDate")])
cohort = left_join(cohort, sample[, c("sampleId", "pathologyPurity")])
cohort = left_join(cohort, msi)
cohort = left_join(cohort, svTypes)
cohort = left_join(cohort, somatics)
cohort = left_join(cohort, patientIds)
cohort = left_join(cohort, unsupportedSegments)
cohort = left_join(cohort, wholeGenomeDuplications)


# Enrich
cohort$qcScore <- round(cohort$unsupportedSegments / cohort$ploidy,0)
cohort$unsupportedSegments <- NULL
cohort$duplicatedAutosomes <- ifelse(is.na(cohort$duplicatedAutosomes), 0, cohort$duplicatedAutosomes)
cohort$WGD <- ifelse(is.na(cohort$WGD), FALSE, cohort$WGD)

# Clean
cohort["_CLONAL"] <- NULL
cohort["_SUBCLONAL"] <- NULL
cohort["_UNKNOWN"] <- NULL


write.csv(cohort, file = '~/hmf/cohortOverview.csv')
