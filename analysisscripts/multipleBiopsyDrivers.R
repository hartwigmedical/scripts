library(dplyr)
detach("package:purple", unload=TRUE);
library(purple)

load(file = "~/hmf/RData/Processed/fragileGenes.RData")
load(file = "~/hmf/RData/Processed/driverGenes.RData")
load(file = "~/hmf/RData/processed/genePanel.RData")
load(file = "~/hmf/RData/Processed/geneCopyNumberDeleteTargets.RData")
load(file = "~/hmf/RData/Processed/geneCopyNumberAmplificationTargets.RData")
load(file = "~/hmf/RData/reference/canonicalTranscripts.RData")

load(file = "~/hmf/RData/reference/multipleBiopsyCohort.RData")
load(file = "~/hmf/RData/reference/multipleBiopsyScope.RData")
load(file = "~/hmf/RData/Reference/mbTertPromoters.Rdata")
load(file = "~/hmf/RData/Processed/mbFusions.RData")
load(file = "~/hmf/RData/Reference/mbGeneCopyNumberAmplifications.RData")
load(file = "~/hmf/RData/Reference/mbGeneCopyNumberDeletes.RData")
load(file = "~/hmf/RData/Processed/mbDndsTsgDrivers.RData")
load(file = "~/hmf/RData/Processed/mbDndsOncoDrivers.RData")

fusions = purple::driver_fusions(mbFusions %>% mutate(shared = scope == 'Shared'), tsGenes, oncoGenes)
amplifications = purple::driver_amplifications(mbGeneCopyNumberAmplifications %>% mutate(shared = scope == 'Shared'), tsGenes, oncoGenes, geneCopyNumberAmplificationTargets)
deletions = purple::driver_deletions(mbGeneCopyNumberDeletes %>% mutate(shared = scope == 'Shared'), tsGenes, oncoGenes, geneCopyNumberDeleteTargets, fragileGenes)
tertPromoters = purple::driver_promoters(mbTertPromoters %>% mutate(shared = scope == 'Shared'))
tsgDriverByGene = mbDndsTsgDrivers %>% select(sampleId, gene, impact, driver, driverLikelihood = driverLikelihoodAdjusted, type, biallelic, hotspot, clonality, shared)
oncoDriverByGene = mbDndsOncoDrivers %>% 
  mutate(hotspot = hotspot | nearHotspot) %>%
  select(sampleId, gene, impact, driver, driverLikelihood = driverLikelihoodAdjusted, type, hotspot, clonality, shared)

driverFactors = c("Fusion-Intragenic","Fusion-Coding","Fusion-UTR","Del","FragileDel","Multihit","Promoter","Frameshift","Nonsense","Splice","Missense","Inframe","Indel","Amp")
mbDriversByGene = bind_rows(oncoDriverByGene, tsgDriverByGene) %>% 
  bind_rows(amplifications) %>% 
  bind_rows(deletions) %>% 
  bind_rows(tertPromoters) %>% 
  bind_rows(fusions) %>%
  mutate(driver = factor(driver, driverFactors)) %>%
  ungroup() %>% group_by(sampleId, gene) %>% 
  top_n(1, driver)

cancerTypes = multipleBiopsyCohort %>% select(sampleId = sampleId, cancerType, patientId)
mbDriversByGene = left_join(mbDriversByGene, cancerTypes, by = "sampleId")

mbDriversByGene = left_join(mbDriversByGene, canonicalTranscripts %>% select(gene, chromosome, start = geneStart, end = geneEnd), by  = "gene")
save(mbDriversByGene, file = "~/hmf/RData/processed/mbDriversByGene.RData")

