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

fusions = purple::driver_fusions(hpcFusions, tsGenes, oncoGenes)
amplifications = purple::driver_amplifications(hpcGeneCopyNumberAmplifications, tsGenes, oncoGenes, geneCopyNumberAmplificationTargets)
deletions = purple::driver_deletions(hpcGeneCopyNumberDeletes, tsGenes, oncoGenes, geneCopyNumberDeleteTargets, fragileGenes)
tertPromoters = purple::driver_promoters(hpcTertPromoters)
tsgDriverByGene = hpcDndsTsgDrivers %>% select(sampleId, gene, impact, driver, driverLikelihood = driverLikelihoodAdjusted, type, biallelic, hotspot, clonality, shared)
oncoDriverByGene = hpcDndsOncoDrivers %>% 
  mutate(hotspot = hotspot | nearHotspot) %>%
  select(sampleId, gene, impact, driver, driverLikelihood = driverLikelihoodAdjusted, type, hotspot, clonality, shared)


driverFactors = c("Fusion-Intragenic","Fusion-Coding","Fusion-UTR","Del","FragileDel","Multihit","Promoter","Frameshift","Nonsense","Splice","Missense","Inframe","Indel","Amp")
hpcDriversByGene = bind_rows(oncoDriverByGene, tsgDriverByGene) %>% 
  bind_rows(amplifications) %>% 
  bind_rows(deletions) %>% 
  bind_rows(tertPromoters) %>% 
  bind_rows(fusions) %>%
  mutate(driver = factor(driver, rev(driverFactors))) %>%
  ungroup() %>% group_by(sampleId, gene) %>% 
  top_n(1, driverLikelihood) %>% 
  top_n(1, driver)

hpcDriversByGene%>%  ungroup() %>% group_by(sampleId, gene) %>% summarise(n = n()) %>% filter(n > 1)

cancerTypes = highestPurityCohort %>% select(sampleId = sampleId, cancerType)
hpcDriversByGene = left_join(hpcDriversByGene, cancerTypes, by = "sampleId")

hpcDriversByGene = left_join(hpcDriversByGene, canonicalTranscripts %>% select(gene, chromosome, start = geneStart, end = geneEnd), by  = "gene")
save(hpcDriversByGene, file = "~/hmf/RData/processed/hpcDriversByGene.RData")

