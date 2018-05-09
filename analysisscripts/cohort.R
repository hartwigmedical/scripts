detach("package:purple", unload=TRUE)
library(purple)
library(RMySQL)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(multidplyr)

### DATABASE
dbProd = dbConnect(MySQL(), dbname='hmfpatients_20180418', groups="RAnalysis")
dbDisconnect(dbProd)
rm(dbProd)

cat("Querying gene deletes")
geneDeletes = purple::query_gene_deletes(dbProd)
save(geneDeletes, file = '~/hmf/RData/reference/geneDeletes.RData')
#load('~/hmf/RData/reference/geneDeletes.RData')

cat("Querying purity")
highestPurityCohort = purple::query_highest_purity_cohort(dbProd, geneDeletes)
save(highestPurityCohort, file = "~/hmf/RData/reference/highestPurityCohort.RData")

cat("Copy Numbers")
highestPurityCopyNumbers = purple::query_copy_number(dbProd, highestPurityCohort)
save(highestPurityCopyNumbers, file = "~/hmf/RData/reference/highestPurityCopyNumbers.RData")

cat("Querying canonical transcripts")
canonicalTranscripts = purple::query_canonical_transcript(dbProd)
save(canonicalTranscripts, file = "~/hmf/RData/reference/canonicalTranscripts.RData")

cat("Querying gene copy number deletes")
geneCopyNumberDeletes = purple::query_gene_copy_number_deletes(dbProd, highestPurityCohort)
geneCopyNumberDeletes = left_join(geneCopyNumberDeletes, highestPurityCohort %>% select(sampleId, cancerType = primaryTumorLocation), by = "sampleId")
save(geneCopyNumberDeletes, file = "~/hmf/RData/reference/geneCopyNumberDeletes.RData")

cat("Querying gene copy number amplifications")
geneCopyNumberAmplifications = purple::query_gene_copy_number_amplifications(dbProd, highestPurityCohort)
geneCopyNumberAmplifications = left_join(geneCopyNumberAmplifications, highestPurityCohort %>% select(sampleId, cancerType = primaryTumorLocation), by = "sampleId")
save(geneCopyNumberAmplifications, file = "~/hmf/RData/reference/geneCopyNumberAmplifications.RData")

cat("Querying somatics")
highestPuritySomatics_p1 = purple::query_somatic_variants(dbProd, highestPurityCohort[1:1000, ])
save(highestPuritySomatics_p1, file = "~/hmf/RData/reference/highestPuritySomatics_p1.RData")
highestPuritySomatics_p2 = purple::query_somatic_variants(dbProd, highestPurityCohort[1001:nrow(highestPurityCohort), ])
save(highestPuritySomatics_p2, file = "~/hmf/RData/reference/highestPuritySomatics_p2.RData")

cat("Find any missing somatics")
somaticSamples = c(unique(highestPuritySomatics_p1$sampleId), unique(highestPuritySomatics_p2$sampleId))
missingSomaticSamples = highestPurityCohort %>% filter(!sampleId %in% somaticSamples)
#highestPuritySomatics_p3 = purple::query_somatic_variants(dbProd, missingSomaticSamples)
#save(highestPuritySomatics_p3, file = "~/hmf/RData/reference/highestPuritySomatics_p3.RData")

cat("Determing exonic somatics")
load(file = "~/hmf/RData/reference/HmfRefCDS.RData")
exonic_somatics <- function(somatics, gr_genes) {
  gr_muts = GRanges(somatics$chromosome, IRanges(somatics$position,somatics$position))
  ol = as.matrix(findOverlaps(gr_muts, gr_genes, type="any", select="all"))
  return (somatics[unique(ol[, 1]), ])
}
exonic_p1 = exonic_somatics(highestPuritySomatics_p1, gr_genes)
exonic_p2 = exonic_somatics(highestPuritySomatics_p2, gr_genes)
#exonic_p3 = exonic_somatics(highestPuritySomatics_p3, gr_genes)
highestPurityExonicSomatics = rbind(exonic_p1, exonic_p2)
save(highestPurityExonicSomatics, file = "~/hmf/RData/reference/highestPurityExonicSomatics.RData")
rm(exonic_p1, exonic_p2)

cat("Somatic cohort level stats")
#load(file = "~/hmf/RData/reference/highestPuritySomatics_p1.RData")
#load(file = "~/hmf/RData/reference/highestPuritySomatics_p2.RData")
#load(file = "~/hmf/RData/reference/highestPuritySomatics_p3.RData")
somatics_summary_p1 = cohort_somatic_summary(highestPuritySomatics_p1)
somatics_summary_p2 = cohort_somatic_summary(highestPuritySomatics_p2)
#somatics_summary_p3 = cohort_somatic_summary(highestPuritySomatics_p3)
highestPuritySomaticSummary = rbind(somatics_summary_p1, somatics_summary_p2)
save(highestPuritySomaticSummary, file = "~/hmf/RData/reference/highestPuritySomaticSummary.RData")
rm(somatics_summary_p1, somatics_summary_p2, somatics_summary_p3)

cat("Structual Variant Overview")
highestPurityStructuralVariantSummary = query_structural_variant_summary(dbProd, highestPurityCohort)
save(highestPurityStructuralVariantSummary, file = "~/hmf/RData/reference/highestPurityStructuralVariantSummary.RData")

cat("Tert promoters")
tertPromoters = purple::query_tert_promoters(dbProd, highestPurityCohort)
save(tertPromoters, file = "~/hmf/RData/reference/tertPromoters.RData")

cat("Fusions")
fusions = purple::query_fusions(dbProd, highestPurityCohort)
save(fusions, file = "~/hmf/RData/reference/fusions.RData")

cat("Whole Genome Duplication")
wgd = purple::query_whole_genome_duplication(dbProd, highestPurityCohort)
save(wgd, file = "~/hmf/RData/reference/wholeGenomeDuplication.RData")

clinicalData = purple::query_clinical_data(dbProd)
save(clinicalData, file = "~/hmf/RData/reference/clinicalData.RData")

sampleData = purple::query_sample_data(dbProd)
save(sampleData, file = "~/hmf/RData/reference/sampleData.RData")

#### COMBINE
load(file = "~/hmf/RData/reference/highestPurityCohort.RData")
load(file = "~/hmf/RData/reference/wholeGenomeDuplication.RData")
load(file = "~/hmf/RData/reference/clinicalData.RData")
load(file = "~/hmf/RData/reference/sampleData.RData")
load(file = "~/hmf/RData/reference/highestPuritySomaticSummary.RData")
load(file = "~/hmf/RData/reference/highestPurityStructuralVariantSummary.RData")

clinicalSummary = clinicalData %>% select(sampleId, cancerSubtype, biopsyDate, biopsySite, biopsyType, biopsyLocation, treatment, treatmentType, birthYear)
cohortSummary = left_join(highestPurityCohort, clinicalSummary, by = "sampleId") %>%
  left_join(sampleData, by = "sampleId") %>%
  left_join(highestPuritySomaticSummary, by = "sampleId") %>%
  left_join(highestPurityStructuralVariantSummary, by = "sampleId") %>%
  left_join(wgd, by = "sampleId") %>%  mutate(duplicatedAutosomes = ifelse(is.na(duplicatedAutosomes), 0, duplicatedAutosomes), WGD = ifelse(is.na(WGD), F, WGD))

save(cohortSummary, file = "~/hmf/RData/processed/cohortSummary.RData")  

############## MULTIPLE BIOPSY COHORT
multipleBiopsyCohort = purple::query_multiple_biopsy_cohort(dbProd)
save(multipleBiopsyCohort, file = "~/hmf/RData/reference/multipleBiopsyCohort.RData")

multipleBiopsyScope = multipleBiopsyCohort %>% 
  left_join(clinicalData %>% select(sampleId, sampleArrivalDate), by = "sampleId") %>%
  arrange(patientId, sampleArrivalDate) %>%
  select(patientId, sampleId) %>% 
  group_by(patientId) %>% 
  mutate(scope = paste0("Sample",row_number()) ) %>%
  ungroup()
save(multipleBiopsyScope, file = "~/hmf/RData/reference/multipleBiopsyScope.RData")

multipleBiopsyStructuralVariants = query_structural_variants(dbProd, multipleBiopsyCohort)
save(multipleBiopsyStructuralVariants, file = "~/hmf/RData/reference/multipleBiopsyStructuralVariants.RData")

multipleBiopsySomatics = purple::query_somatic_variants(dbProd, multipleBiopsyCohort)
save(multipleBiopsySomatics, file = "~/hmf/RData/reference/multipleBiopsySomatics.RData")

multipleBiopsyStructuralVariantsWithScope = left_join(multipleBiopsyStructuralVariants, multipleBiopsyScope, by = "sampleId") %>%
  group_by(patientId, startChromosome, endChromosome, startPosition,endPosition,startOrientation,endOrientation,type) %>%
  mutate(scope = ifelse(n_distinct(sampleId) > 1, "Shared", scope))
save(multipleBiopsyStructuralVariantsWithScope, file = "~/hmf/RData/reference/multipleBiopsyStructuralVariantsWithScope.RData")

multipleBiopsySomaticsWithScope = multipleBiopsySomatics %>% 
  filter(filter == 'PASS') %>%
  left_join(multipleBiopsyScope, by = "sampleId") %>%
  partition(patientId, chromosome, position, ref, alt, type) %>%
  mutate(scope = ifelse(n_distinct(sampleId) > 1, "Shared", scope)) %>%
  collect() %>% 
  as_tibble() %>%
  ungroup()
save(multipleBiopsySomaticsWithScope, file = "~/hmf/RData/reference/multipleBiopsySomaticsWithScope.Rdata")

multipleBiopsyMSI = purple::cohort_msi(multipleBiopsySomaticsWithScope %>% ungroup())
save(multipleBiopsyMSI, file = "~/hmf/RData/reference/multipleBiopsyMSI.RData")

load(file = "~/hmf/RData/reference/multipleBiopsyMSI.RData")
load(file = "~/hmf/RData/reference/multipleBiopsyScope.RData")
load(file = "~/hmf/RData/reference/multipleBiopsyCohort.RData")
load(file = "~/hmf/RData/reference/multipleBiopsySomaticsWithScope.Rdata")
load(file = "~/hmf/RData/reference/multipleBiopsyStructuralVariantsWithScope.RData")


multipleBiopsyStructuralVariantSummary = multipleBiopsyStructuralVariantsWithScope %>%
  ungroup() %>% 
  distinct(patientId, scope, startChromosome, endChromosome, startPosition, endPosition, startOrientation, endOrientation, type) %>% 
  group_by(patientId, type, scope) %>% summarise(count = n()) %>% 
  unite(type, scope, type, sep = "_") %>% 
  spread(type, count)

multipleBiopsySomaticVariantSummary = multipleBiopsySomaticsWithScope %>%
  ungroup() %>% 
  distinct(patientId, scope, chromosome, position, ref, alt, type, clonality) %>% 
  group_by(patientId, scope, type, clonality) %>%
  summarise(count = n()) %>%
  unite(type, scope, type, clonality, sep = "_") %>% spread(type, count)

multipleBiopsyPatientMsi = multipleBiopsyMSI %>% left_join(multipleBiopsyScope, by = "sampleId") %>%
  gather(type, value, msiScore, msiStatus) %>%
  unite(type, scope, type) %>% 
  select(-sampleId) %>%
  spread(type, value)

multipleBiopsySampleIds = multipleBiopsyScope %>% spread(scope, sampleId)

multipleBiopsyPurity = multipleBiopsyCohort %>% select(patientId, sampleId, purity, ploidy) %>% 
  left_join(multipleBiopsyScope %>% select(sampleId, scope), by = "sampleId") %>% 
  gather(type, value, purity, ploidy) %>%  
  unite(type, scope, type) %>%
  select(-sampleId) %>%
  spread(type, value)

multipleBiopsySummary = multipleBiopsyCohort %>% select(patientId, sampleId, gender, primaryTumorLocation) %>% 
  left_join(multipleBiopsyScope %>% select(sampleId, scope), by = "sampleId") %>% 
  filter(scope == "Sample1") %>% select(patientId, gender, primaryTumorLocation) %>%
  left_join(multipleBiopsySampleIds, by = "patientId") %>%
  left_join(multipleBiopsyPurity, by = "patientId") %>%
  left_join(multipleBiopsyPatientMsi, by = "patientId") %>%
  left_join(multipleBiopsyStructuralVariantSummary, by = "patientId") %>%
  left_join(multipleBiopsySomaticVariantSummary, by = "patientId")
save(multipleBiopsySummary, file = "~/hmf/RData/processed/multipleBiopsySummary.RData")


#### VISUALISATION
load(file = "~/hmf/RData/reference/highestPurityCohort.RData")
cohortByPrimaryTumorLocation = highestPurityCohort %>% group_by(primaryTumorLocation) %>% summarise(N = n())
save(cohortByPrimaryTumorLocation, file = '~/hmf/RData/reference/cohortByPrimaryTumorLocation.RData')
primaryTumorLocations = unique(highestPurityCohort$primaryTumorLocation)
primaryTumorLocations= primaryTumorLocations[!is.na(primaryTumorLocations)]

cosmicSignatureColours = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
                           "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a",
                           "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
                           "#dea185","#a0729d","#8a392f")
primaryTumorLocationColours = setNames(cosmicSignatureColours[1:length(primaryTumorLocations)], primaryTumorLocations)
save(primaryTumorLocationColours, file = "~/hmf/RData/reference/primaryTumorLocationColours.RData")


