load(file = "~/hmf/analysis/svPaper/cohort.RData")
load(file = "~/hmf/analysis/svPaper/cohortCancerTypes.RData")
load(file = "~/hmf/analysis/svPaper/hpcCopyNumbers.RData")

cohortCancerTypes = cohortCancerTypes %>% ungroup()

centromeres = hpcCopyNumbers %>% filter(segmentStartSupport == 'CENTROMERE') %>% select(chromosome, centromere = start) %>% distinct()

enrichedCopyNumbers = hpcCopyNumbers %>% 
  select(sampleId, chromosome, start, end, minorAllelePloidy, copyNumber) %>%
  left_join(centromeres, by = "chromosome") %>%
  left_join(highestPurityCohort %>% select(sampleId, ploidy), by = "sampleId") %>% 
  left_join(cohortCancerTypes %>% select(sampleId, primaryTumorLocation), by = "sampleId") %>%
  filter(!is.na(primaryTumorLocation)) %>%
  mutate(
    chromosome = factor(chromosome, levels = c(1:22, "X", "Y"), ordered = T),
    length = end - start + 1,
    arm = ifelse(start < centromere, "p", "q"),
    isLOH = minorAllelePloidy < 0.5, 
    isGain = copyNumber > 1.4 * ploidy, 
    isLoss = copyNumber < 0.6 * ploidy) %>%
  group_by(primaryTumorLocation, sampleId, chromosome, arm) %>%
  mutate(
    armLength = sum(length)) 

primaryTumorLocationTotals = enrichedCopyNumbers %>% ungroup() %>% distinct(primaryTumorLocation, sampleId) %>% group_by(primaryTumorLocation) %>% summarise(observations = n())

summarisedData = enrichedCopyNumbers %>%
  group_by(primaryTumorLocation, sampleId, chromosome, arm) %>%
  summarise(
    lohProportion = sum(isLOH * length) / max(armLength),
    gainProportion = sum(isGain * length) / max(armLength),
    lossProportion = sum(isLoss * length) / max(armLength)) %>%
  mutate(
    isLOH = lohProportion > 0.75,
    isGain = gainProportion > 0.75,
    isLoss = lossProportion > 0.75
  ) %>% 
  group_by(primaryTumorLocation, chromosome, arm) %>%
  summarise(
    loh = sum(isLOH), 
    gain = sum(isGain),
    loss = sum(isLoss)) %>%
  left_join(primaryTumorLocationTotals, by = "primaryTumorLocation")

#### VERIFICATION
nrow(centromeres) # 24
nrow(enrichedCopyNumbers %>% ungroup() %>% select(chromosome, arm, armLength) %>% distinct) # 48  
sum(primaryTumorLocationTotals$observations) # 3462
length(unique(enrichedCopyNumbers$sampleId)) # 3462
sum(summarisedData %>% ungroup() %>% distinct(primaryTumorLocation, observations) %>% select(observations)) # 3462
sum(summarisedData$loh > summarisedData$observations) #0
sum(summarisedData$gain > summarisedData$observations) #0
sum(summarisedData$gain > summarisedData$observations) #0



write.csv(summarisedData, file = "~/hmf/analysis/ArmGainAndLoss.csv", row.names = F)
