
### REFERENCE
load(file = "~/hmf/RData/Reference/allPurity.RData")

dbProd = dbConnect(MySQL(), dbname='hmfpatients_20180418', groups="RAnalysis")
allFusions = purple::query_fusions(dbProd, allPurity)
save(allFusions, file = "~/hmf/RData/reference/allFusions.RData")
dbDisconnect(dbProd)
rm(dbProd, allPurity)

dbEnsemble = dbConnect(MySQL(), dbname='homo_sapiens_core_89_37', groups="RAnalysis")
allFusionCodingRegions = purple::query_coding_regions(dbEnsemble, unique(c(allFusions$`5pTranscript`, allFusions$`3pTranscript`)))
save(allFusionCodingRegions, file = "~/hmf/RData/reference/allFusionCodingRegions.RData")
dbDisconnect(dbEnsemble)
rm(dbEnsemble)


### PROCESSED
load(file = "~/hmf/RData/reference/allFusions.RData")
load(file = "~/hmf/RData/reference/allFusionCodingRegions.RData")

fusions = allFusions %>%
  mutate(`5pPosition` = ifelse(`5pIsStartEnd` == 1, startPosition, endPosition)) %>%
  left_join(allFusionCodingRegions %>% select(`5pTranscript` = transcriptId, `5pCodingStart` = coding_start, `5pCodingEnd` = coding_end), by = "5pTranscript")  %>%
  mutate(
    `5pPreCoding` = ifelse((`5pStrand` == 1 & `5pPosition` < `5pCodingStart`) | (`5pStrand` == -1 & `5pPosition` > `5pCodingEnd`), T, F),
    `5pPostCoding` = ifelse((`5pStrand` == 1 & `5pPosition` > `5pCodingEnd`) | (`5pStrand` == -1 & `5pPosition` < `5pCodingStart`), T, F)) %>%
  select(-`5pPosition`, -`5pCodingStart`, -`5pCodingEnd`) %>%
  mutate(`3pPosition` = ifelse(`3pIsStartEnd` == 1, startPosition, endPosition)) %>%
  left_join(allFusionCodingRegions %>% select(`3pTranscript` = transcriptId, `3pCodingStart` = coding_start, `3pCodingEnd` = coding_end), by = "3pTranscript") %>%
  mutate(
    `3pPreCoding` = ifelse((`3pStrand` == 1 & `3pPosition` < `3pCodingStart`) | (`3pStrand` == -1 & `3pPosition` > `3pCodingEnd`), T, F),
    `3pPostCoding` = ifelse((`3pStrand` == 1 & `3pPosition` > `3pCodingEnd`) | (`3pStrand` == -1 & `3pPosition` < `3pCodingStart`), T, F)) %>%
  select(-`3pPosition`, -`3pCodingStart`, -`3pCodingEnd`) %>%
  mutate(intragenic = `5pGene` == `3pGene`) %>%
  filter(`3pPostCoding` == F, is.na(`5pPostCoding`) | `5pPostCoding` == F, !(intragenic & `5pPhase` == -1)) %>%
  mutate(driver = "Fusion-Coding", driver = ifelse(is.na(`5pPreCoding`) | `5pPreCoding`, "Fusion-UTR", driver), driver = ifelse(intragenic, "Fusion-Intragenic", driver)) %>%
  select(-ends_with("IsStartEnd"),  -ends_with("Coding"), -filter, -intragenic)

load(file = "~/hmf/RData/reference/highestPurityCohort.RData")
hpcFusions = fusions %>% 
  filter(sampleId %in% highestPurityCohort$sampleId) %>%
  select(-starts_with("start"), -starts_with("end"))
save(hpcFusions, file = "~/hmf/RData/Processed/hpcFusions.RData")
  
load(file = "~/hmf/RData/reference/multipleBiopsyCohort.RData")
load(file = "~/hmf/RData/reference/multipleBiopsyScope.RData")
mbFusions = fusions %>% 
  filter(sampleId %in% multipleBiopsyCohort$sampleId) %>%
  left_join(multipleBiopsyScope, by = "sampleId") %>%
  group_by(patientId, startChromosome, endChromosome, startPosition,endPosition) %>%
  mutate(scope = ifelse(n_distinct(sampleId) > 1, "Shared", scope)) %>%
  ungroup() %>%
  select(-starts_with("start"), -starts_with("end"))
save(mbFusions, file = "~/hmf/RData/Processed/mbFusions.RData")

