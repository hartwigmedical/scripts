detach("package:purple", unload=TRUE);
library(purple);
library(RMySQL)
library(data.table)
library(MutationalPatterns)
library(ggplot2)

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
purity = purple::query_purity(dbProd)
clinical = purple::query_clinical_data(dbProd)[, c("sampleId", "cancerType")]
sampleCohort = merge(purity, clinical, by=c("sampleId"), all.x = TRUE)

patientIdLookups = query_patient_id_lookup(dbProd)
patientIds = purple::apply_to_cohort(sampleCohort, function(x) {purple::sample_to_patient_id(x$sampleId, patientIdLookups)})
colnames(patientIds) <- c("sampleId", "patientId")

sampleCohort = merge(sampleCohort, patientIds, by=c("sampleId"), all.x = TRUE)
patientCohort = highest_purity_patients(sampleCohort)
rm(purity)
rm(patientIds)
rm(patientIdLookups)
rm(clinical)

patientVariants = purple::query_somatic_variants(dbProd, patientCohort)
save(sampleCohort, patientCohort, patientVariants, file = "~/hmf/RData/patientVariants.RData")
load(file = "~/hmf/RData/patientVariants.RData")

dbDisconnect(dbProd)
rm(dbProd)

snps = patientVariants[patientVariants$type == 'SNP', ]
snps = purple::add_scope_to_variants(snps, by = c("chromosome", "position", "alt", "ref", "type"))

######### FIND COMMON VARIANTS BETWEEN TWO SAMPLES
firstSampleId = "CPCT02370006T"
secondSampleId = "CPCT02160030T"
cat("Processing", firstSampleId, "with", secondSampleId, "\n")
combinedVariants = snps[snps$sampleId == firstSampleId | snps$sampleId == secondSampleId,]
sharedVariants = combinedVariants[, .N,  by=c("chromosome", "position", "alt", "ref", "type")]
sharedVariantCount = sum(sharedVariants[sharedVariants$N > 1, ]$N)

######### FIND REMAINING DUPLICATES
snps = data.table(patientVariants[patientVariants$type == 'SNP', ])

snps$TwoOccurences <- F
snps[snps[, .I[.N == 2], by=c("chromosome", "position", "alt", "ref", "type")]$V1, ]$TwoOccurences <- T

variantsWithExactly2Occurences = snps[snps$TwoOccurences == T, ] 
duplicatedVariants = variantsWithExactly2Occurences[, .(paste(sort(.SD$sampleId), collapse = "|")), by=c("chromosome", "position", "alt", "ref", "type")]
duplicatedVariants = duplicatedVariants[, .N, by=c("V1")]
foo <- data.frame(do.call('rbind', strsplit(as.character(duplicatedVariants$V1),'|',fixed=TRUE)))
foo$N <- duplicatedVariants$N
duplicatedVariants = foo
names(duplicatedVariants) <- c("left", "right", "N")

significantDuplications = duplicatedVariants[duplicatedVariants$N > 400, ]
leftApparences = significantDuplications[, .N, by=left]
names(leftApparences) <- c("sample", "N")
rightApparences = significantDuplications[, .N, by=right]
names(rightApparences) <- c("sample", "N")
