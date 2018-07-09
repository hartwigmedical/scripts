detach("package:purple", unload=TRUE)
library(purple)
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
load(file = "~/hmf/RData/Processed/hpcDndsTsgDrivers.RData")
load(file = "~/hmf/RData/Processed/hpcDndsOncoDrivers.RData")


fusions = purple::driver_fusions(hpcFusions, tsGenes, oncoGenes) %>% mutate(category = "Fusion")
amplifications = purple::driver_amplifications(hpcGeneCopyNumberAmplifications, tsGenes, oncoGenes, geneCopyNumberAmplificationTargets) %>% mutate(category = "CNA")
deletions = purple::driver_deletions(hpcGeneCopyNumberDeletes, tsGenes, oncoGenes, geneCopyNumberDeleteTargets, fragileGenes) %>% mutate(category = "CNA")
tertPromoters = purple::driver_promoters(hpcTertPromoters) %>% mutate(category = "Promoter")
tsgDriverByGene = hpcDndsTsgDrivers %>% 
  mutate(hotspot = ifelse(hotspot, "Hotspot", "None")) %>%
  select(sampleId, coordinate, gene, category = variant, driver, impact, driverLikelihood = driverLikelihoodAdjusted, type, biallelic, hotspot, clonality, shared, pHGVS)
oncoDriverByGene = hpcDndsOncoDrivers %>% 
  mutate(hotspot = ifelse(hotspot, "Hotspot", "None"),
         hotspot = ifelse(nearHotspot, "NearHotspot", hotspot)) %>%
  select(sampleId, coordinate, gene, category = variant, driver, impact, driverLikelihood = driverLikelihoodAdjusted, type, hotspot, clonality, shared, pHGVS)

hotspotFactors = c("Hotspot","NearHotspot","None")
driverFactors = c("Fusion-Intragenic","Fusion-Coding","Fusion-UTR","Del","FragileDel","Multihit","Promoter","Frameshift","Nonsense","Splice","Missense","Inframe","Indel","Amp")
hpcDriversByGene = bind_rows(oncoDriverByGene, tsgDriverByGene) %>% 
  bind_rows(amplifications) %>% 
  bind_rows(deletions) %>% 
  bind_rows(tertPromoters) %>% 
  bind_rows(fusions) %>%
  mutate(
    hotspot = factor(hotspot, rev(hotspotFactors)),
    driver = factor(driver, rev(driverFactors))) %>%
  ungroup() %>% group_by(sampleId, gene) %>% 
  top_n(1, driverLikelihood) %>% 
  top_n(1, driver) %>%
  ungroup()

hpcDriversByGene %>%ungroup() %>% group_by(sampleId, gene) %>% summarise(n = n()) %>% filter(n > 1)

cancerTypes = highestPurityCohort %>% select(sampleId = sampleId, cancerType)
hpcDriversByGene = left_join(hpcDriversByGene, cancerTypes, by = "sampleId")

driverCatalogue = left_join(hpcDriversByGene, canonicalTranscripts %>% select(gene, chromosome, geneStart = geneStart, geneEnd = geneEnd), by  = "gene") %>%
  select(sampleId, cancerType, gene, type, chromosome, geneStart, geneEnd, coordinate, category, driver, impact, driverLikelihood,  clonality, pHGVS, hotspot, biallelic)

write.csv(driverCatalogue, file = "~/hmf/RData/DriverCatalogue.csv", row.names = F) 
save(hpcDriversByGene, file = "~/hmf/RData/processed/hpcDriversByGene.RData")


