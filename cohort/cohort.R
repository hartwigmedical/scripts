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

hmfIds = read.csv(file = "~/hmf/resources/idgenerator/amber_anonymized.csv") %>%
  select(sampleId = OriginalId, hmfId = AnonymousId) %>%
  mutate(patientId = substr(hmfId, 1, 9)) %>% select(sampleId, patientId)

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis", host = "127.0.0.1")
purpleCohort = purple::query_cohort(dbProd)
purpleClinical = purple::query_clinical_data(dbProd)
dbDisconnect(dbProd)
rm(dbProd)

patientCohort = purpleCohort %>% select(-patientId) %>% inner_join(hmfIds, by = "sampleId")
highestPurityCohort = purple::highest_purity_cohort(patientCohort) %>% select(-patientId)

load(file = "~/hmf/RData/Reference/allClinicalData.RData")
paperClinical = allClinicalData %>% select(sampleId, paperCancerType = cancerType)

highestPurityCohortClinical = purpleClinical %>% 
  filter(sampleId %in% highestPurityCohort$sampleId) %>% 
  select(sampleId, primaryTumorLocation) %>%
  left_join(paperClinical, by = "sampleId") %>%
  mutate(
    primaryTumorLocation = ifelse(primaryTumorLocation == "Unknown", NA, primaryTumorLocation),
    cancerType = coalesce(primaryTumorLocation, paperCancerType, "Unknown"),
    cancerType = ifelse(cancerType == "Bone/soft tissue", "Bone/Soft tissue", cancerType),
    cancerType = ifelse(cancerType == "Head and Neck", "Head and neck", cancerType),
    cancerType = ifelse(cancerType == 'Net', 'NET', cancerType),
    cancerType = ifelse(cancerType == 'Small intestine', 'Small Intestine', cancerType)
    ) %>% 
  group_by(cancerType) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  mutate(cancerType = ifelse(n < 10, "Other", cancerType)) %>%
  select(sampleId, cancerType)

cancerTypeCounts = highestPurityCohortClinical %>% group_by(cancerType) %>% count() %>% arrange(n)
cancerTypeCounts

highestPurityCohort = highestPurityCohort %>% left_join(highestPurityCohortClinical, by = "sampleId")
save(highestPurityCohort, file = "~/hmf/analysis/genepanel/highestPurityCohort.RData")
load(file = "~/hmf/analysis/genepanel/highestPurityCohort.RData")


########### STRUCTURAL VARIANTS ########### 
sampleIdString = paste("'", highestPurityCohort$sampleId, "'", collapse = ",", sep = "")
svQuery = paste0("select * from structuralVariant where filter = 'PASS' AND sampleId in (",sampleIdString, ")")
hpcStructuralVariants = dbGetQuery(dbProd, svQuery)
save(hpcStructuralVariants, file = "~/hmf/analysis/genepanel/hpcStructuralVariants.RData")


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

