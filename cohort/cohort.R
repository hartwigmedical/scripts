#### 1. cohort
#### 2. dndsRefresh
#### 3. dndsClassification
#### 4. dndDriverLikelihood
#### 5. genePanel
#### 6. amps and dels
#### 7. amp and del targets

library(dplyr)
library(tidyr)
library(RMySQL)

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis", host = "127.0.0.1")

hmfIds = dbGetQuery(dbProd, "SELECT sampleId, hmfId from sampleMapping") %>%
  mutate(patientId = substr(hmfId, 1, 9)) %>% select(sampleId, patientId)
purpleCohort = purple::query_cohort(dbProd) %>%  
  select(-patientId) %>% inner_join(hmfIds, by = "sampleId")
purpleClinical = purple::query_clinical_data(dbProd) %>% 
  filter(sampleId %in% purpleCohort$sampleId) %>%
  select(-patientId) %>% left_join(hmfIds, by = "sampleId")

dbDisconnect(dbProd)
rm(dbProd)

load(file = "~/hmf/RData/Reference/allClinicalData.RData")
paperClinical = allClinicalData %>% select(-patientId) %>%
  filter(sampleId %in% purpleCohort$sampleId, (cancerType != 'Other' | primaryTumorLocation == "Double primary")) %>%
  inner_join(hmfIds, by = "sampleId") %>%
  select(patientId, paperCancerType = cancerType)

paperClinical %>% group_by(patientId) %>% count() %>% filter(n() > 1)

sampleClinical = purpleClinical %>% ungroup() %>%
  select(sampleId, patientId, primaryTumorLocation) %>%
  left_join(paperClinical, by = "patientId") %>%
  distinct() %>%
  mutate(
    primaryTumorLocation = ifelse(primaryTumorLocation == "Unknown", NA, primaryTumorLocation),
    cancerType = coalesce(paperCancerType, primaryTumorLocation, "Unknown"),
    cancerType = ifelse(cancerType == "Bone/soft tissue", "Bone/Soft tissue", cancerType),
    cancerType = ifelse(cancerType == "Head and Neck", "Head and neck", cancerType),
    cancerType = ifelse(cancerType == 'Net', 'NET', cancerType),
    cancerType = ifelse(cancerType == 'Small intestine', 'Small Intestine', cancerType),
    cancerType = ifelse(cancerType == 'Nervous system', 'CNS', cancerType)
  ) %>%
  filter(cancerType != "Unknown") %>%
  group_by(patientId) %>% mutate(nSamples = n()) %>% group_by(patientId, cancerType) %>% mutate(nSamples = max(nSamples),  nCancerType = n())

sampleClinical %>% filter(nSamples != nCancerType)

patientClinical = sampleClinical %>%
  ungroup() %>%
  mutate(cancerType = ifelse(nSamples != nCancerType, "Other", cancerType)) %>%
  group_by(patientId, cancerType) %>% 
  summarise() %>%
  group_by(cancerType) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  mutate(cancerType = ifelse(n < 10, "Other", cancerType)) %>%
  select(patientId, cancerType)

cohort = purpleCohort %>% left_join(patientClinical, by = "patientId") %>% mutate(cancerType = ifelse(is.na(cancerType), "Unknown", cancerType))
highestPurityCohort = purple::highest_purity_cohort(cohort) %>% select(-patientId)

multipleBiopsyMapping = cohort %>% group_by(patientId) %>% 
  filter(n() > 1) %>% 
  select(patientId, sampleId) %>% 
  mutate(sample = paste0("sample", row_number())) %>%
  spread(sample, sampleId) %>%
  ungroup() %>%
  mutate(patientId = row_number())
  
cohort = cohort %>% mutate(hpc = sampleId %in% highestPurityCohort$sampleId) %>% select(-patientId)

save(highestPurityCohort, file = "~/hmf/analysis/cohort/reference/highestPurityCohort.RData")
save(cohort, highestPurityCohort, multipleBiopsyMapping, file = "~/hmf/analysis/cohort/reference/cohort.RData")

load(file = "~/hmf/analysis/cohort/highestPurityCohort.RData")
cancerTypeCounts = highestPurityCohort %>% group_by(cancerType) %>% count() %>% arrange(n)
sum(!is.na(multipleBiopsyMapping %>% select(-patientId, -Sample1)))

########### STRUCTURAL VARIANTS ########### 
sampleIdString = paste("'", highestPurityCohort$sampleId, "'", collapse = ",", sep = "")
svQuery = paste0("select * from structuralVariant where filter = 'PASS' AND sampleId in (",sampleIdString, ")")
hpcStructuralVariants = dbGetQuery(dbProd, svQuery)
save(hpcStructuralVariants, file = "~/hmf/analysis/cohort/hpcStructuralVariants.RData")


svQuery = paste0("select * from structuralVariant where filter = 'INFERRED' AND sampleId in (",sampleIdString, ")")
hpcInferredStructuralVariants = dbGetQuery(dbProd, svQuery)
save(hpcInferredStructuralVariants, file = "~/hmf/analysis/cohort/hpcInferredStructuralVariants.RData")

########### Somatics Counts ########### 
sampleIdString = paste("'", highestPurityCohort$sampleId, "'", collapse = ",", sep = "")
somaticCountsQuery = paste0("select sampleId, type, count(*) as total from somaticVariant where filter = 'PASS' and sampleId in (",sampleIdString, ") group by sampleId, type")
somaticCounts = dbGetQuery(dbProd, somaticCountsQuery)
somaticCounts = somaticCounts %>% mutate(type = paste0("TOTAL_", type)) %>% spread(type, total)
save(somaticCounts, file = "~/hmf/analysis/genepanel/somaticCounts.RData")


########### PRIOR SOMATICS ########### 
load(file = "~/hmf/RData/reference/HmfRefCDS.RData")
exonic_somatics <- function(somatics, gr_genes) {
  gr_muts = GRanges(somatics$chromosome, IRanges(somatics$position,somatics$position + nchar(somatics$ref) - 1))
  ol = as.matrix(findOverlaps(gr_muts, gr_genes, type="any", select="all"))
  return (somatics[unique(ol[, 1]), ])
}

load(file = "~/hmf/analysis/dnds/somatics1.RData")
somatics1 = somatics1 %>% filter(sampleId %in% highestPurityCohort$sampleId)
somatics1 = somatics1 %>% filter(!sampleId %in% somatic3samples)
somatic1samples = unique(somatics1$sampleId)
save(somatics1, file = "~/hmf/analysis/genepanel/somatics1.RData")
exonicSomatics1 = exonic_somatics(somatics1, gr_genes)
save(exonicSomatics1, file = "~/hmf/analysis/genepanel/exonicSomatics1.RData")

load(file = "~/hmf/analysis/dnds/somatics2.RData")
somatics2 = somatics2 %>% filter(sampleId %in% highestPurityCohort$sampleId)
somatics2 = somatics2 %>% filter(!sampleId %in% somatic3samples)
somatic2samples = unique(somatics2$sampleId)
save(somatics2, file = "~/hmf/analysis/genepanel/somatics2.RData")
exonicSomatics2 = exonic_somatics(somatics2, gr_genes)
save(exonicSomatics2, file = "~/hmf/analysis/genepanel/exonicSomatics2.RData")

########### SOMATIC VARIANTS ########### 


load(file = "~/hmf/RData/Reference/hpcExonicSomatics.RData")
samplesWithSomatics = unique(hpcExonicSomatics$sampleId)
samplesWithoutSomatics = highestPurityCohort %>% filter(!sampleId %in% samplesWithSomatics) %>% pull(sampleId)

sampleIdString = paste("'", samplesWithoutSomatics, "'", collapse = ",", sep = "")
somaticQuery = paste0("select * from somaticVariant where filter = 'PASS' AND sampleId in (",sampleIdString, ")")
somatics3 = dbGetQuery(dbProd, somaticQuery)
save(somatics3, file = "~/hmf/analysis/genepanel/somatics3.RData")

exonicSomatics3 = exonic_somatics(somatics3, gr_genes)
save(exonicSomatics3, file = "~/hmf/analysis/genepanel/exonicSomatics3.RData")

load(file = "~/hmf/analysis/genepanel/hpcExonicSomatics.RData")
missingExonicsCohort = highestPurityCohort %>% filter(!sampleId %in% hpcExonicSomatics$sampleId)
unwantedExonicsSomatics = hpcExonicSomatics %>% filter(!sampleId %in% highestPurityCohort$sampleId)

hpcExonicSomatics = bind_rows(bind_rows(exonicSomatics1, exonicSomatics2), exonicSomatics3)
save(hpcExonicSomatics, file = "~/hmf/analysis/genepanel/hpcExonicSomatics.RData")

########### COPY NUMBER ########### 
hpcGeneCopyNumberDeletes = purple::query_gene_copy_number_deletes(dbProd, highestPurityCohort)
hpcGeneCopyNumberDeletes = left_join(hpcGeneCopyNumberDeletes, highestPurityCohort %>% select(sampleId, cancerType), by = "sampleId")
save(hpcGeneCopyNumberDeletes, file = "~/hmf/analysis/genepanel/hpcGeneCopyNumberDeletes.RData")

hpcGeneCopyNumberAmplifications = purple::query_gene_copy_number_amplifications(dbProd, highestPurityCohort)
hpcGeneCopyNumberAmplifications = left_join(hpcGeneCopyNumberAmplifications, highestPurityCohort %>% select(sampleId, cancerType), by = "sampleId")
save(hpcGeneCopyNumberAmplifications, file = "~/hmf/analysis/genepanel/hpcGeneCopyNumberAmplifications.RData")

hpcCopyNumbers = purple::query_copy_number(dbProd, highestPurityCohort)
save(hpcCopyNumbers, file = "~/hmf/analysis/genepanel/hpcCopyNumbers.RData")

########### VERIFY PATIENTS ########### 
load(file = "~/hmf/analysis/genepanel/hpcStructuralVariants.RData")
load(file = "~/hmf/analysis/genepanel/highestPurityCohort.RData")

svVerification = hpcStructuralVariants %>% group_by(startChromosome, endChromosome, startOrientation, endOrientation, startPosition, endPosition) %>%
  mutate(n = n()) %>%
  filter(n > 2)

