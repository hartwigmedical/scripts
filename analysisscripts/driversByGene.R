library(dplyr)
library(tidyr)
library(ggplot2)

load(file = "~/hmf/RData/Reference/canonicalTranscripts.RData")
load(file = "~/hmf/RData/Reference/hpcTertPromoters.Rdata")
load(file = "~/hmf/RData/Reference/hpcGeneCopyNumberDeletes.RData")
load(file = "~/hmf/RData/Reference/hpcGeneCopyNumberAmplifications.RData")
load(file = "~/hmf/RData/Processed/geneCopyNumberDeleteTargets.RData")
load(file = "~/hmf/RData/Processed/geneCopyNumberAmplificationTargets.RData")
load(file = "~/hmf/RData/Processed/driverGenes.RData")
load(file = "~/hmf/RData/processed/tsgDrivers.RData")
load(file = "~/hmf/RData/processed/oncoDrivers.RData")
load(file = "~/hmf/RData/reference/highestPurityCohort.RData")
load(file = "~/hmf/RData/processed/genePanel.RData")
load(file = "~/hmf/RData/Processed/fragileGenes.RData")
load(file = "~/hmf/RData/Processed/hpcFusions.RData")

fusions = hpcFusions %>% select(sampleId, gene = `3pGene`, driver) %>% mutate(type = "FUSION")
fusions$type = ifelse(fusions$gene %in% tsGenes$gene_name, "TSG", fusions$type)
fusions$type = ifelse(fusions$gene %in% oncoGenes$gene_name, "ONCO", fusions$type)

amplifications = hpcGeneCopyNumberAmplifications %>%
  filter(gene %in% oncoGenes$gene_name | gene %in% geneCopyNumberAmplificationTargets$target) %>%
  group_by(sampleId = sampleId, gene) %>% summarise(driver = "Amp") %>%
  mutate(driverLikelihood = 1, type = ifelse(gene %in% tsGenes$gene, "TSG", "ONCO"))

deletions = hpcGeneCopyNumberDeletes %>%
  filter(gene %in% tsGenes$gene_name | gene %in% geneCopyNumberDeleteTargets$target) %>%
  filter(germlineHetRegions == 0, germlineHomRegions == 0)
  group_by(sampleId = sampleId, gene) %>% summarise(driver = "Del", partial = somaticRegions > 1) %>% 
  left_join(fragileGenes %>% select(gene = gene_name, fragile), by = "gene") %>%
  mutate(fragile = ifelse(is.na(fragile), F, T)) %>%
  mutate(driverLikelihood = 1, type = ifelse(gene %in% oncoGenes$gene, "ONCO", "TSG"))%>%
  mutate(driver = ifelse(fragile, "FragileDel", driver))

tertPromoters = hpcTertPromoters %>% 
  group_by(sampleId, gene) %>% 
  summarise(driver = "Promoter", driverLikelihood = 1, type = "ONCO") 

tsgDriverByGene = tsgDrivers %>%
  group_by(sampleId, gene) %>%
  summarise(
    multihit = n() > 1,
    biallelic = sum(biallelic) > 0, 
    hotspot = sum(hotspot) > 0, 
    impact = paste(impact, collapse = ","),
    driver = max(driver),
    driverLikelihood = max(driver)) %>%
  mutate(driver = ifelse(multihit, "Multihit", impact), type = 'TSG') %>%
  select(sampleId, gene, impact, driver, driverLikelihood, type, biallelic)

oncoDriverByGene = oncoDrivers %>%
  group_by(sampleId, gene) %>%
  summarise(
    hotspot = sum(hotspot) > 0,
    nearHotspot = any(nearHotspot),
    impact = paste(impact, collapse = ","),
    driver = max(driver),
    driverLikelihood = max(driver))  %>%
  mutate(driver = impact, type = 'ONCO') %>%
  select(sampleId, gene, impact, driver, driverLikelihood, type)

driverFactors = c("IntragenicFusion","3PrimeFusion","5PrimeFusion","Del","PartialDel","FragileDel","Multihit","Promoter","Frameshift","Nonsense","Splice","Missense","Inframe","Amp")
driversByGene = bind_rows(oncoDriverByGene, tsgDriverByGene) %>% 
  bind_rows(amplifications) %>% 
  bind_rows(deletions) %>% 
  bind_rows(tertPromoters) %>% 
  bind_rows(fusions) %>%
  mutate(driver = factor(driver, driverFactors)) %>%
  ungroup() %>% group_by(sampleId, gene) %>% 
  top_n(1, driver)

driversByGene%>%  ungroup() %>% group_by(sampleId, gene) %>% summarise(n = n()) %>% filter(n > 1)

cancerTypes = highestPurityCohort %>% select(sampleId = sampleId, cancerType)
driversByGene = left_join(driversByGene, cancerTypes, by = "sampleId")

driversByGene = left_join(driversByGene, canonicalTranscripts %>% select(gene, chromosome, start = geneStart, end = geneEnd), by  = "gene")
save(driversByGene, file = "~/hmf/RData/processed/driversByGene.RData")

