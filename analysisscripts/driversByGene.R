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
load(file = "~/hmf/RData/Reference/highestPurityCohort.RData")
load(file = "~/hmf/RData/processed/genePanel.RData")
load(file = "~/hmf/RData/Processed/fragileGenes.RData")
load(file = "~/hmf/RData/Processed/hpcFusions.RData")
load(file = "~/hmf/RData/Processed/hpcOncoDrivers.RData")
load(file = "~/hmf/RData/Processed/hpcTsgDrivers.RData")

fusions = hpcFusions %>% select(sampleId, gene = `3pGene`, driver) %>% 
  mutate(type = "FUSION") %>% 
  distinct() %>%
  mutate(
    type = ifelse(gene %in% tsGenes$gene_name, "TSG", type),
    type = ifelse(gene %in% oncoGenes$gene_name, "ONCO", type),
    driverLikelihood = 1)

amplifications = hpcGeneCopyNumberAmplifications %>%
  filter(gene %in% oncoGenes$gene_name | gene %in% geneCopyNumberAmplificationTargets$target) %>%
  group_by(sampleId = sampleId, gene) %>% summarise(driver = "Amp") %>%
  mutate(
    driverLikelihood = 1, 
    type = ifelse(gene %in% tsGenes$gene, "TSG", "ONCO"),
    biallelic = T)

deletionArms = geneCopyNumberDeleteTargets %>% mutate(arm = coalesce(telomere, centromere)) %>% filter(!is.na(arm) ) %>% select(gene, arm)
deletions = hpcGeneCopyNumberDeletes %>%
  filter(gene %in% tsGenes$gene_name | gene %in% geneCopyNumberDeleteTargets$target) %>%
  filter(germlineHetRegions == 0, germlineHomRegions == 0) %>%
  group_by(sampleId = sampleId, gene) %>% summarise(driver = "Del", partial = somaticRegions > 1) %>% 
  left_join(fragileGenes %>% select(gene = gene_name, fragile), by = "gene") %>%
  mutate(fragile = ifelse(is.na(fragile), F, T)) %>%
  mutate(
    driverLikelihood = 1, 
    type = ifelse(gene %in% oncoGenes$gene, "ONCO", "TSG"),
    driver = ifelse(fragile, "FragileDel", driver),
    biallelic = T) %>%
  left_join(deletionArms, by = "gene") %>%
  mutate(gene = coalesce(arm, gene)) %>%
  select(-arm)

tertPromoters = hpcTertPromoters %>% 
  group_by(sampleId, gene) %>% 
  summarise(
    driver = "Promoter", 
    driverLikelihood = 1, 
    type = "ONCO") 

tsgDriverByGene = hpcTsgDrivers %>% select(sampleId, gene, impact, driver, driverLikelihood = driverLikelihoodAdjusted, type, biallelic, hotspot, clonality)
oncoDriverByGene = hpcOncoDrivers %>% select(sampleId, gene, impact, driver, driverLikelihood = driverLikelihoodAdjusted, type, hotspot, clonality)

driverFactors = c("Fusion-Intragenic","Fusion-Coding","Fusion-UTR","Del","FragileDel","Multihit","Promoter","Frameshift","Nonsense","Splice","Missense","Inframe","Indel","Amp")
hpcDriversByGene = bind_rows(oncoDriverByGene, tsgDriverByGene) %>% 
  bind_rows(amplifications) %>% 
  bind_rows(deletions) %>% 
  bind_rows(tertPromoters) %>% 
  bind_rows(fusions) %>%
  mutate(driver = factor(driver, driverFactors)) %>%
  ungroup() %>% group_by(sampleId, gene) %>% 
  top_n(1, driver)

hpcDriversByGene%>%  ungroup() %>% group_by(sampleId, gene) %>% summarise(n = n()) %>% filter(n > 1)

cancerTypes = highestPurityCohort %>% select(sampleId = sampleId, cancerType)
hpcDriversByGene = left_join(hpcDriversByGene, cancerTypes, by = "sampleId")

hpcDriversByGene = left_join(hpcDriversByGene, canonicalTranscripts %>% select(gene, chromosome, start = geneStart, end = geneEnd), by  = "gene")
save(hpcDriversByGene, file = "~/hmf/RData/processed/hpcDriversByGene.RData")

