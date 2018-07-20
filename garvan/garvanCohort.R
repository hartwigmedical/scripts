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


