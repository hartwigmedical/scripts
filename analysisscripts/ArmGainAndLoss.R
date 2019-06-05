library(tidyr)
library(dplyr)

load(file = "~/hmf/analysis/genepanel/highestPurityCohort.RData")
load(file = "~/hmf/analysis/genepanel/hpcCopyNumbers.RData")

cohortCancerTypes = unique(highestPurityCohort$cancerType)

centromeres = hpcCopyNumbers %>% filter(segmentStartSupport == 'CENTROMERE') %>% select(chromosome, centromere = start) %>% distinct()

enrichedCopyNumbers = hpcCopyNumbers %>% 
  select(sampleId, chromosome, start, end, minorAllelePloidy, copyNumber) %>%
  left_join(centromeres, by = "chromosome") %>%
  left_join(highestPurityCohort %>% select(sampleId, ploidy, cancerType), by = "sampleId") %>% 
  filter(!is.na(cancerType)) %>%
  mutate(
    chromosome = factor(chromosome, levels = c(1:22, "X", "Y"), ordered = T),
    length = end - start + 1,
    arm = ifelse(start < centromere, "p", "q"),
    isLOH = minorAllelePloidy < 0.5, 
    isGain = copyNumber > 1.4 * ploidy, 
    isLoss = copyNumber < 0.6 * ploidy,
    isLossAndLOH = isLoss & isLOH) %>%
  group_by(cancerType, sampleId, chromosome, arm) %>%
  mutate(
    armLength = sum(length)) 

cancerTypeTotals = enrichedCopyNumbers %>% ungroup() %>% distinct(cancerType, sampleId) %>% group_by(cancerType) %>% summarise(observations = n())

summarisedData = enrichedCopyNumbers %>%
  group_by(cancerType, sampleId, chromosome, arm) %>%
  summarise(
    lohProportion = sum(isLOH * length) / max(armLength),
    gainProportion = sum(isGain * length) / max(armLength),
    lossProportion = sum(isLoss * length) / max(armLength),
    lossAndLOHProportion = sum(isLossAndLOH * length) / max(armLength)) %>%
  mutate(
    isLOH = lohProportion > 0.75,
    isGain = gainProportion > 0.75,
    isLoss = lossProportion > 0.75,
    isLossAndLOH = lossAndLOHProportion > 0.75
  ) %>% 
  group_by(cancerType, chromosome, arm) %>%
  summarise(
    loh = sum(isLOH), 
    gain = sum(isGain),
    loss = sum(isLoss),
    lossAndLOH = sum(isLossAndLOH)) %>%
  left_join(cancerTypeTotals, by = "cancerType")

#### VERIFICATION
nrow(centromeres) # 24
nrow(enrichedCopyNumbers %>% ungroup() %>% select(chromosome, arm, armLength) %>% distinct) # 48  
sum(cancerTypeTotals$observations) # 3524
length(unique(enrichedCopyNumbers$sampleId)) # 3524
sum(summarisedData %>% ungroup() %>% distinct(cancerType, observations) %>% select(observations)) # 3524
sum(summarisedData$loh > summarisedData$observations) #0
sum(summarisedData$gain > summarisedData$observations) #0
sum(summarisedData$gain > summarisedData$observations) #0
nrow(summarisedData %>% filter(lossAndLOH > loss)) #0



write.csv(summarisedData, file = "~/hmf/analysis/ArmGainAndLoss.csv", row.names = F)
