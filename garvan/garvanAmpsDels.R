detach("package:purple", unload=TRUE)
library(purple)
library(dplyr)

outputDir = "~/garvan/RData/"
outputDir = "~/Documents/LKCGP_projects/RData/"

referenceDir = paste0(outputDir, "reference/")
processedDir = paste0(outputDir, "processed/")


####### LOAD DATA FROM FILE #######
load(paste0(processedDir, "highestPurityCohortSummary.RData"))
load(paste0(referenceDir, "cohortGeneCopyNumberDeletes.RData"))
load(paste0(referenceDir, "cohortGeneCopyNumberAmplifications.RData"))
load(paste0(referenceDir, "canonicalTranscripts.RData"))
genes = canonicalTranscripts %>% select(gene, chromosome, start = geneStart, end = geneEnd)

####### EXECUTE ALGORITHM #######
geneCopyNumberDeletes = cohortGeneCopyNumberDeletes %>% 
  left_join(highestPurityCohortSummary %>% select(sampleId, cancerType), by = "sampleId") %>% filter(
  !is.na(cancerType),
  germlineHetRegions == 0,
  germlineHomRegions == 0,
  !(minRegionStartSupport == "TELOMERE" & minRegionEndSupport == "CENTROMERE"),
  !(minRegionStartSupport == "CENTROMERE" & minRegionEndSupport == "TELOMERE"))


geneCopyNumberAmplifications = cohortGeneCopyNumberAmplifications %>% 
  left_join(highestPurityCohortSummary %>% select(sampleId, cancerType), by = "sampleId") %>% 
  filter(germlineHetRegions == 0, germlineHomRegions == 0, !is.na(cancerType))

deletionsOutput = copy_number_drivers(genes, geneCopyNumberDeletes, adjacent = "arm", absCandidateScore = 2, relativeCandidateScore = 0.85); date()
amplificationOutput = copy_number_drivers(genes, geneCopyNumberAmplifications, adjacent = "arm", absCandidateScore = 5, relativeCandidateScore = 0.75); date()

geneCopyNumberAmplificationSummary = amplificationOutput$summary
geneCopyNumberDeletionsSummary = deletionsOutput$summary

save(geneCopyNumberAmplificationSummary, file = paste0(processedDir,"geneCopyNumberAmplificationSummary.RData"))
save(geneCopyNumberDeletionsSummary, file = paste0(processedDir,"geneCopyNumberDeletionsSummary.RData"))

