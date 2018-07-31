detach("package:purple", unload=TRUE)
library(purple)
library(dplyr)
library(tidyr)
library(ggplot2)

outputDir = "~/garvan/RData/"
outputDir = "~/Documents/LKCGP_projects/RData/"
referenceDir = paste0(outputDir, "reference/")
processedDir = paste0(outputDir, "processed/")

load(file = paste0(referenceDir, "canonicalTranscripts.RData"))
load(file = paste0(referenceDir, "cohortTertPromoters.Rdata"))
load(file = paste0(referenceDir, "cohortGeneCopyNumberDeletes.RData"))
load(file = paste0(referenceDir, "cohortGeneCopyNumberAmplifications.RData"))

load(file = paste0(processedDir, "highestPurityCohortSummary.RData"))
load(file = paste0(processedDir, "geneCopyNumberDeleteTargets.RData"))
load(file = paste0(processedDir, "geneCopyNumberAmplificationTargets.RData"))
load(file = paste0(processedDir, "driverGenes.RData"))
load(file = paste0(processedDir, "fragileGenes.RData"))
load(file = paste0(processedDir, "hpcFusions.RData"))
load(file = paste0(processedDir, "hpcDndsTsgDrivers.RData"))
load(file = paste0(processedDir, "hpcDndsOncoDrivers.RData"))

#Clean up data
geneCopyNumberDeleteTargets$centromere <- as.character(geneCopyNumberDeleteTargets$centromere)

# Run
fusions = purple::driver_fusions(hpcFusions, tsGenes, oncoGenes) %>% mutate(category = "Fusion")
amplifications = purple::driver_amplifications(cohortGeneCopyNumberAmplifications, tsGenes, oncoGenes, geneCopyNumberAmplificationTargets) %>% mutate(category = "CNA")
deletions = purple::driver_deletions(cohortGeneCopyNumberDeletes, tsGenes, oncoGenes, geneCopyNumberDeleteTargets, fragileGenes) %>% mutate(category = "CNA")
tertPromoters = purple::driver_promoters(cohortTertPromoters) %>% mutate(category = "Promoter")
tsgDriverByGene = hpcDndsTsgDrivers %>%
    mutate(hotspot = ifelse(hotspot, "Hotspot", "NonHotspot")) %>%
    select(sampleId, coordinate, gene, category = variant, driver, impact, driverLikelihood = driverLikelihoodAdjusted, type, biallelic, hotspot, clonality, shared, pHGVS)
oncoDriverByGene = hpcDndsOncoDrivers %>%
    mutate(hotspot = ifelse(hotspot, "Hotspot", "NonHotspot"),
    hotspot = ifelse(nearHotspot, "NearHotspot", hotspot)) %>%
    select(sampleId, coordinate, gene, category = variant, driver, impact, driverLikelihood = driverLikelihoodAdjusted, type, hotspot, clonality, shared, pHGVS)

hotspotFactors = c("Hotspot","NearHotspot","NonHotspot")
driverFactors = c("Fusion-Intragenic","Fusion-Coding","Fusion-UTR","Del","FragileDel","Multihit","Promoter","Frameshift","Nonsense","Splice","Missense","Inframe","Indel","Amp")
hpcDriversByGene = bind_rows(oncoDriverByGene, tsgDriverByGene) %>%
    bind_rows(amplifications) %>%
    bind_rows(deletions) %>%
    bind_rows(fusions) %>%
    mutate(
    hotspot = factor(hotspot, rev(hotspotFactors)),
    driver = factor(driver, rev(driverFactors))) %>%
    ungroup() %>% group_by(sampleId, gene) %>%
    top_n(1, driverLikelihood) %>%
    top_n(1, driver) %>%
    ungroup() %>%
    filter(sampleId %in% highestPurityCohortSummary$sampleId)

hpcDriversByGene %>%ungroup() %>% group_by(sampleId, gene, type) %>% summarise(n = n()) %>% filter(n > 1)

cancerTypes = highestPurityCohortSummary %>% select(sampleId = sampleId, cancerType)
hpcDriversByGene = left_join(hpcDriversByGene, cancerTypes, by = "sampleId") %>%
left_join(canonicalTranscripts %>% select(gene, chromosome, geneStart = geneStart, geneEnd = geneEnd), by  = "gene")
save(hpcDriversByGene, file = paste0(processedDir, "hpcDriverCatalog.RData"))

