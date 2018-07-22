library(dplyr)
library(tidyr)

#outputDir = "~/garvan/RData/"
outputDir = "~/Documents/LKCGP_projects/RData/"
referenceDir = paste0(outputDir, "reference/")
processedDir = paste0(outputDir, "processed/")

load(paste0(referenceDir,"cohort.RData"))
load(paste0(referenceDir,"cohortSomatics.RData"))
load(paste0(referenceDir,"cohortStructuralVariantSummary.RData"))
load(paste0(referenceDir,"cohortWGD.RData"))
load(paste0(referenceDir,"samples.RData"))

load(paste0(processedDir, "cohortIndelSummary.RData"))
load(paste0(processedDir, "cohortSNPSummary.RData"))
load(paste0(processedDir, "cohortSomaticsSummary.RData"))

clinicalSummary = samples %>% 
  mutate(cancer_type = ifelse(cancer_type == 'Other ', 'Other', cancer_type)) %>%
  transmute(sampleId = sample_id, patientId = patient_id, cancerType = cancer_type, ageAtSample = age_at_sample)

cohortSummary = cohort %>%
  select(sampleId, gender, status, qcStatus, purity, ploidy) %>%
  left_join(clinicalSummary, by = "sampleId") %>%
  left_join(cohortSomaticsSummary, by = "sampleId") %>%
  left_join(cohortStructuralVariantSummary, by = "sampleId") %>%
  left_join(cohortWGD, by = "sampleId") %>%
  mutate(
    duplicatedAutosomes = ifelse(is.na(duplicatedAutosomes), 0, duplicatedAutosomes), 
    WGD = ifelse(is.na(WGD), F, WGD))
save(cohortSummary, file = paste0(processedDir, "cohortSummary.RData"))

highestPurityCohortSummary = cohortSummary %>% filter(!is.na(cancerType))
save(highestPurityCohortSummary, file = paste0(processedDir, "highestPurityCohortSummary.RData"))

hpcCancerTypeCounts = highestPurityCohortSummary %>%
  group_by(cancerType) %>%
  summarise(
    N = n(),
    medianMutationalLoad = median(TOTAL_SNV + TOTAL_INDEL)) %>%
  arrange(medianMutationalLoad)
save(hpcCancerTypeCounts, file = paste0(processedDir, "hpcCancerTypeCounts.RData"))

cancerTypes = sort(unique(highestPurityCohortSummary$cancerType))
cancerTypeColours = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
                      "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#ff4791","#01837a",
                      "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
                      "#dea185","#a0729d","#8a392f")
cancerTypeColours = setNames(cancerTypeColours[1:length(cancerTypes)], cancerTypes)
save(cancerTypeColours, file = paste0(processedDir, "cancerTypeColours.RData"))


simplifiedDrivers = c("Amp","Del","FragileDel","Fusion","Indel","Missense","Multihit","Nonsense","Promoter","Splice")
simplifiedDriverColours = c("#fb8072","#bc80bd","#bebada", "#fdb462","#80b1d3","#8dd3c7","#b3de69","#fccde5","#ffffb3","#d9d9d9")
simplifiedDriverColours = setNames(simplifiedDriverColours, simplifiedDrivers)
save(simplifiedDrivers, simplifiedDriverColours, file = paste0(processedDir, "simplifiedDrivers.RData"))
