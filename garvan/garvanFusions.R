library(RMySQL)
library(dplyr)

outputDir = "~/garvan/RData/"
resourceDir = "~/hmf/resources/"

outputDir = "~/Documents/LKCGP_projects/RData/"
resourceDir = "/Users/marwo2/Documents/LKCGP_projects/RData/Resources/"

referenceDir = paste0(outputDir, "reference/")
processedDir = paste0(outputDir, "processed/")
load(file = paste0(referenceDir, "cohortFusions.RData"))
load(file = paste0(referenceDir, "cohortFusionCodingRegions.RData"))
load(file = paste0(processedDir, "highestPurityCohortSummary.RData"))

fusionConstraints = read.csv(file = paste0(resourceDir, "3PrimeFusionDomainConstraintsV1.csv"), stringsAsFactors = F) %>%
    filter(!is.na(`Min.Exon`) | !is.na(`Max.Exon`)) %>%
    select(`3pTranscript` = Transcript, `3pMinExonConstraint` = `Min.Exon`, `3pMaxExonConstraint` = `Max.Exon`)

fusions = cohortFusions %>%
# Exon constraints
    left_join(fusionConstraints, by = "3pTranscript") %>%
    filter(is.na(`3pMinExonConstraint`) | `3pExonRankDownstream` >= `3pMinExonConstraint`) %>%
    filter(is.na(`3pMaxExonConstraint`) | `3pExonRankDownstream` <= `3pMaxExonConstraint`) %>%
    select(-`3pMinExonConstraint`, -`3pMaxExonConstraint`) %>%

# Coding Constraints
    mutate(`5pPosition` = ifelse(`5pIsStartEnd` == 1, startPosition, endPosition)) %>%
    left_join(cohortFusionCodingRegions %>% select(`5pTranscript` = transcriptId, `5pCodingStart` = coding_start, `5pCodingEnd` = coding_end), by = "5pTranscript")  %>%
    mutate(
    `5pPreCoding` = ifelse((`5pStrand` == 1 & `5pPosition` < `5pCodingStart`) | (`5pStrand` == -1 & `5pPosition` > `5pCodingEnd`), T, F),
    `5pPostCoding` = ifelse((`5pStrand` == 1 & `5pPosition` > `5pCodingEnd`) | (`5pStrand` == -1 & `5pPosition` < `5pCodingStart`), T, F)) %>%
    select(-`5pPosition`, -`5pCodingStart`, -`5pCodingEnd`) %>%
    mutate(`3pPosition` = ifelse(`3pIsStartEnd` == 1, startPosition, endPosition)) %>%
    left_join(cohortFusionCodingRegions %>% select(`3pTranscript` = transcriptId, `3pCodingStart` = coding_start, `3pCodingEnd` = coding_end), by = "3pTranscript") %>%
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

hpcFusions = fusions %>%
    filter(sampleId %in% highestPurityCohortSummary$sampleId) %>%
    select(-starts_with("start"), -starts_with("end"))

save(hpcFusions, file = paste0(processedDir, "hpcFusions.RData"))
