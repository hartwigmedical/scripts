library(dplyr)
library(tidyr)
library(ggplot2)

load(file = "~/hmf/RData/input/geneCopyNumberDeletes.RData")
load(file = "~/hmf/RData/input/geneCopyNumberAmplifications.RData")
load(file = "~/hmf/RData/output/geneCopyNumberDeleteTargets.RData")
load(file = "~/hmf/RData/output/geneCopyNumberAmplificationTargets.RData")
load(file = "~/hmf/RData/output/driverGenes.RData")
load(file = "~/hmf/RData/output/genePanel.RData")
genePanel = genePanel %>% filter(martincorena | hmf | cosmicCurated)


allAmplifications = geneCopyNumberAmplifications %>%
  filter(gene %in% oncoGenes$gene_name | gene %in% geneCopyNumberAmplificationTargets$target) %>%
  select(sampleId = sampleId, gene) %>% mutate(amp = T)

targetAmplifications = geneCopyNumberAmplifications %>%
  filter(gene %in% geneCopyNumberAmplificationTargets$target) %>%
  select(sampleId = sampleId, gene) %>% mutate(amp = T)

allDeletions = geneCopyNumberDeletes %>%
  filter(gene %in% tsGenes$gene_name | gene %in% geneCopyNumberDeleteTargets$target) %>%
  select(sampleId = sampleId, gene) %>% mutate(del = T)

targetDeletions = geneCopyNumberDeletes %>%
  filter(gene %in% geneCopyNumberDeleteTargets$target) %>%
  select(sampleId = sampleId, gene) %>% mutate(del = T)

load(file = "~/hmf/RData/output/tsgDrivers.RData")
tsgDriverByGene = tsgDrivers %>%
  mutate(driver = ifelse(hotspot, 1, driver)) %>%
  group_by(sampleId, gene) %>%
  summarise(
    multihit = n() > 1,
    biallelic = sum(biallelic) > 0, 
    hotspot = sum(hotspot) > 0, 
    nearHotspot = F,
    impact = paste(impact, collapse = ","),
    driver = max(driver)) %>%
  full_join(allDeletions, by = c("sampleId" ,"gene")) %>% 
  left_join(targetAmplifications, by = c("sampleId" ,"gene")) %>% 
  mutate(driver = ifelse(!is.na(del), 1, driver)) %>%
  mutate(type = 'TSG')

load(file = "~/hmf/RData/output/oncoDrivers.RData")
oncoDriverByGene = oncoDrivers %>%
  group_by(sampleId, gene) %>%
  summarise(
    multihit = F,
    biallelic = F,
    hotspot = sum(hotspot) > 0,
    nearHotspot = any(nearHotspot),
    impact = paste(impact, collapse = ","),
    driver = max(driver))  %>%
  full_join(allAmplifications, by = c("sampleId" ,"gene")) %>% 
  left_join(targetDeletions, by = c("sampleId" ,"gene")) %>% 
  mutate(driver = ifelse(!is.na(amp), 1, driver)) %>%
  mutate(type = 'ONCO')

load(file = "~/hmf/RData/input/tertPromoters.Rdata")
tertPromoters = tertPromoters %>% group_by(sampleId, gene) %>% summarise(promoter = T)

driversByGene = bind_rows(oncoDriverByGene, tsgDriverByGene) %>%
  full_join(tertPromoters, by = c("sampleId", "gene")) %>%
  mutate(driver = ifelse(!is.na(promoter), 1, driver)) %>%
  mutate(type = ifelse(is.na(type) & gene == 'TERT', "TSG", type))

load(file = "~/hmf/RData/input/highestPurityCohort.RData")
cancerTypes = highestPurityCohort %>% select(sampleId = sampleId, primaryTumorLocation)
driversByGene = left_join(driversByGene, cancerTypes, by = "sampleId")
save(driversByGene, file = "~/hmf/RData/output/driversByGene.RData")

driversByGene %>% group_by(sampleId, gene) %>% filter(any(del) & any(amp)) 
