detach("package:purple", unload=TRUE)
library(purple)
library(RMySQL)
library(GenomicRanges)
library(dplyr)
library(tidyr)

### DATABASE
dbProd = dbConnect(MySQL(), dbname='hmfpatients_20180418', groups="RAnalysis")
dbDisconnect(dbProd)
rm(dbProd)


cat("Querying canonical transcripts")
canonicalTranscripts = purple::query_canonical_transcript(dbProd)
save(canonicalTranscripts, file = "~/hmf/RData/input/canonicalTranscripts.RData")

cat("Querying purple")
highestPurityCohort = purple::query_highest_purity_cohort(dbProd)
save(highestPurityCohort, file = "~/hmf/RData/input/highestPurityCohort.RData")

cat("Querying somatics")
highestPuritySomatics_p1 = purple::query_somatic_variants(dbProd, highestPurityCohort[1:1000, ])
save(highestPuritySomatics_p1, file = "~/hmf/RData/input/highestPuritySomatics_p1.RData")
highestPuritySomatics_p2 = purple::query_somatic_variants(dbProd, highestPurityCohort[1001:nrow(highestPurityCohort), ])
save(highestPuritySomatics_p2, file = "~/hmf/RData/input/highestPuritySomatics_p2.RData")

cat("Find any missing somatics")
somaticSamples = c(unique(highestPuritySomatics_p1$sampleId), unique(highestPuritySomatics_p2$sampleId))
missingSomaticSamples = highestPurityCohort %>% filter(!sampleId %in% somaticSamples)
highestPuritySomatics_p3 = purple::query_somatic_variants(dbProd, missingSomaticSamples)
save(highestPuritySomatics_p3, file = "~/hmf/RData/input/highestPuritySomatics_p3.RData")

cat("Determing exonic somatics")
load(file = "~/hmf/RData/input/HmfRefCDS.RData")
exonicSomatics <- function(somatics, gr_genes) {
  gr_muts = GRanges(somatics$chromosome, IRanges(somatics$position,somatics$position))
  ol = as.matrix(findOverlaps(gr_muts, gr_genes, type="any", select="all"))
  return (somatics[unique(ol[, 1]), ])
}
exonic_p1 = exonic_somatics(highestPuritySomatics_p1, gr_genes)
exonic_p2 = exonic_somatics(highestPuritySomatics_p2, gr_genes)
exonic_p3 = exonic_somatics(highestPuritySomatics_p3, gr_genes)
highestPurityExonicSomatics = rbind(rbind(exonic_p1, exonic_p2), exonic_p3)
save(highestPurityExonicSomatics, file = "~/hmf/RData/input/highestPurityExonicSomatics.RData")
rm(exonic_p1, exonic_p2, exonic_p3)

cat("Somatic cohort level stats")
somatics_summary_p1 = cohort_somatic_summary(highestPuritySomatics_p1)
somatics_summary_p2 = cohort_somatic_summary(highestPuritySomatics_p2)
somatics_summary_p3 = cohort_somatic_summary(highestPuritySomatics_p3)
highestPuritySomaticSummary = rbind(rbind(somatics_summary_p1, somatics_summary_p2), somatics_summary_p3)
save(highestPuritySomaticSummary, file = "~/hmf/RData/input/highestPuritySomaticSummary.RData")
rm(somatics_summary_p1, somatics_summary_p2, somatics_summary_p3)

cat("Structual Variant Overview")
highestPurityStructuralVariantSummary = query_structural_variant_summary(dbProd, highestPurityCohort)
save(highestPurityStructuralVariantSummary, file = "~/hmf/RData/input/highestPurityStructuralVariantSummary.RData")

cat("Copy Numbers")
highestPurityCopyNumbers = purple::query_copy_number(dbProd, highestPurityCohort)
save(highestPurityCopyNumbers, file = "~/hmf/RData/input/highestPurityCopyNumbers.RData")

tertPromoters = purple::query_tert_promoters(dbProd, highestPurityCohort)
save(tertPromoters, file = "~/hmf/RData/input/tertPromoters.RData")


#### COMBINE
load(file = "~/hmf/RData/input/highestPurityCohort.RData")
load(file = "~/hmf/RData/input/highestPuritySomaticSummary.RData")
load(file = "~/hmf/RData/input/highestPurityStructuralVariantSummary.RData")

cohortSummary = left_join(highestPurityCohort, highestPuritySomaticSummary)


#### VISUALISATION
load(file = "~/hmf/RData/highestPurityCohort.RData")
cohortByPrimaryTumorLocation = highestPurityCohort %>% group_by(primaryTumorLocation) %>% summarise(N = n())
save(cohortByPrimaryTumorLocation, file = '~/hmf/RData/cohortByPrimaryTumorLocation.RData')
primaryTumorLocations = unique(highestPurityCohort$primaryTumorLocation)
primaryTumorLocations= primaryTumorLocations[!is.na(primaryTumorLocations)]

cosmicSignatureColours = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
                           "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a",
                           "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
                           "#dea185","#a0729d","#8a392f")
primaryTumorLocationColours = setNames(cosmicSignatureColours[1:length(primaryTumorLocations)], primaryTumorLocations)
save(primaryTumorLocationColours, file = "~/hmf/RData/primaryTumorLocationColours.RData")



