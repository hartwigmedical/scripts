library(RMySQL)
library(dplyr)

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

transcriptString = paste("'", unique(c(allFusions$`3pTranscript`)), "'", collapse = ",", sep = "")
biotypeQuery = paste0("SELECT stable_id as '3pTranscript', biotype as '3pBiotype'  FROM transcript where stable_id in (",transcriptString, ")")
allFusionBioTypes = dbGetQuery(dbEnsemble, biotypeQuery)
save(allFusionBioTypes, file = "~/hmf/RData/reference/allFusionBioTypes.RData")
rm(biotypeQuery, transcriptString)

dbDisconnect(dbEnsemble)
rm(dbEnsemble)

### PROCESSED
load(file = "~/hmf/RData/reference/allFusions.RData")
load(file = "~/hmf/RData/reference/allFusionCodingRegions.RData")
load(file = "~/hmf/RData/reference/allFusionBioTypes.RData")

fusionConstraints = read.csv(file = "~/hmf/Resources/3PrimeFusionDomainConstraintsV1.csv", stringsAsFactors = F) %>%
  filter(!is.na(`Min.Exon`) | !is.na(`Max.Exon`)) %>%
  select(`3pTranscript` = Transcript, `3pMinExonConstraint` = `Min.Exon`, `3pMaxExonConstraint` = `Max.Exon`)

fusions = allFusions %>%
  # Exon constraints
  left_join(fusionConstraints, by = "3pTranscript") %>%
  filter(is.na(`3pMinExonConstraint`) | `3pExonRankDownstream` >= `3pMinExonConstraint`) %>%
  filter(is.na(`3pMaxExonConstraint`) | `3pExonRankDownstream` <= `3pMaxExonConstraint`) %>%
  select(-`3pMinExonConstraint`, -`3pMaxExonConstraint`) %>%
  
  # Coding Constraints
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
  
  # Clean up
  mutate(
    driver = "Fusion-Coding", 
    driver = ifelse(is.na(`5pPreCoding`) | `5pPreCoding`, "Fusion-UTR", driver), 
    driver = ifelse(intragenic, "Fusion-Intragenic", driver)) %>%
  select(-ends_with("IsStartEnd"),  -ends_with("Coding"), -filter, -intragenic)

load(file = "~/hmf/RData/reference/highestPurityCohort.RData")
hpcFusions = fusions %>% 
  filter(sampleId %in% highestPurityCohort$sampleId) %>%
  select(-starts_with("start"), -starts_with("end")) %>%
  left_join( allFusionBioTypes, by = '3pTranscript')
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

