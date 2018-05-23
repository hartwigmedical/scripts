library(IRanges)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
detach("package:purple", unload=TRUE);
library(purple)
library(RMySQL)

####### Variant Annotation
load(file = "~/hmf/RData/reference/hpcExonicSomatics.RData")
somatics = hpcExonicSomatics %>% 
   select(sampleId, chromosome, position, ref, alt, type, worstCodingEffect, canonicalCodingEffect, hotspot, biallelic, clonality)

load(file = "~/hmf/RData/processed/dndsFilteredAnnotatedMutations.RData")
filteredMutations = dnds_annotate_somatics(dndsFilteredAnnotatedMutations, somatics)

# Gene Panel
load(file = "~/hmf/RData/processed/genePanel.RData")
genePanel = genePanel %>% filter(martincorena | hmf | cosmicCurated)
filteredMutations = filteredMutations %>% filter(gene %in% genePanel$gene_name, impact != "")

filteredTsgMutations = tsg_mutations(filteredMutations)
filteredOncoMutations = onco_mutations(filteredMutations)

summarise_annotation <- function(annotatedMutations) {
  annotatedMutations$impact <- ifelse(annotatedMutations$impact %in% c("Frameshift","Inframe"), "INDEL", annotatedMutations$impact)
  annotatedMutations$geneStatus <- ifelse(annotatedMutations$redundant, "Redundant", annotatedMutations$geneStatus)
  annotatedMutations$impact <- as.character(annotatedMutations$impact)
  annotatedMutations = annotatedMutations %>% group_by(gene, impact, geneStatus) %>% summarise(n = n()) %>% spread(geneStatus, n)
  annotatedMutations[is.na(annotatedMutations)] <- 0
  return (annotatedMutations)
}

filteredTsgMutationsByGene = summarise_annotation(filteredTsgMutations)
filteredOncoMutationsByGene = summarise_annotation(filteredOncoMutations)

load(file = "~/hmf/Rdata/Processed/HmfRefCDSCv.RData")
HmfRefCDSCv = purple::dnds_excess(HmfRefCDSCv)

excessVariants = HmfRefCDSCv %>% 
  filter(cancerType == 'All') %>%
  select(gene = gene_name, Missense = excess_mis, Nonsense = excess_non, Splice = excess_spl, INDEL = excess_ind) %>%
  gather(impact, Excess, Missense, Nonsense, Splice, INDEL) %>%
  mutate(Excess = round(Excess))


excessTsgRates = left_join(filteredTsgMutationsByGene, excessVariants, by = c("gene", "impact"))
excessTsgRates$BiallelicDriverRate = pmax(0,pmin(1.0, (excessTsgRates$Excess) / excessTsgRates$Biallelic))
excessTsgRates$MultiHitDriverRate = pmax(0,pmin(1.0, (excessTsgRates$Excess - excessTsgRates$Biallelic) / excessTsgRates$MultiHit))
excessTsgRates$SingleHitDriverRate = pmax(0,pmin(1.0, (excessTsgRates$Excess - excessTsgRates$Biallelic - excessTsgRates$MultiHit) / (excessTsgRates$SingleHit)))
excessTsgRates[is.na(excessTsgRates)] <- 0
excessTsgInframeRates = excessTsgRates %>% ungroup() %>% filter(impact == "INDEL") %>% mutate(impact = "Inframe")
excessTsgFrameshiftRates = excessTsgRates %>% ungroup() %>% filter(impact == "INDEL") %>% mutate(impact = "Frameshift")
excessTsgRates = excessTsgRates %>% filter(impact != "INDEL") %>% bind_rows(excessTsgInframeRates) %>% bind_rows(excessTsgFrameshiftRates) %>% arrange(gene, impact) %>%
  select(gene, impact, BiallelicDriverRate, MultiHitDriverRate, SingleHitDriverRate)


excessOncoRates = left_join(filteredOncoMutationsByGene, excessVariants, by = c("gene", "impact"))
excessOncoRates$HotspotDriverRate = pmax(0,pmin(1.0, (excessOncoRates$Excess) / (excessOncoRates$Hotspot + excessOncoRates$NearHotspot)))
excessOncoRates$HitDriverRate = pmax(0,pmin(1.0, (excessOncoRates$Excess - excessOncoRates$Hotspot - excessOncoRates$NearHotspot) / (excessOncoRates$Hit)))
excessOncoRates[is.na(excessOncoRates)] <- 0
excessOncoInframeRates = excessOncoRates %>% ungroup() %>% filter(impact == "INDEL") %>% mutate(impact = "Inframe")
excessOncoFrameshiftRates = excessOncoRates %>% ungroup() %>% filter(impact == "INDEL") %>% mutate(impact = "Frameshift")
excessOncoRates = excessOncoRates %>% filter(impact != "INDEL") %>% bind_rows(excessOncoInframeRates) %>% bind_rows(excessOncoFrameshiftRates) %>% arrange(gene, impact) %>%
  select(gene, impact, HotspotDriverRate, HitDriverRate)

save(excessTsgRates, excessOncoRates, file = "~/hmf/RData/processed/excessRates.RData")


####### APPLY EXCESS RATES TO MUTATIONS
load("~/hmf/RData/processed/excessRates.RData")

load(file = "~/hmf/RData/reference/hpcExonicSomatics.RData")
somatics = hpcExonicSomatics %>% 
  select(sampleId, chromosome, position, ref, alt, type, worstCodingEffect, canonicalCodingEffect, hotspot, biallelic, clonality)

load(file = "~/hmf/RData/processed/dndsUnfilteredAnnotatedMutations.RData")
mutations = dnds_annotate_somatics(dndsUnfilteredAnnotatedMutations, somatics)

# Gene Panel
load(file = "~/hmf/RData/processed/genePanel.RData")
genePanel = genePanel %>% filter(martincorena | hmf | cosmicCurated)
mutations = mutations %>% filter(gene %in% genePanel$gene_name, impact != "")

tsgAnnotatedMutations = tsg_mutations(mutations)
oncoAnnotatedMutations = onco_mutations(mutations)

load(file = "~/hmf/RData/Processed/driverGenes.RData")
tsgAnnotatedMutations$geneStatus <- ifelse(tsgAnnotatedMutations$redundant, "Redundant", tsgAnnotatedMutations$geneStatus)
tsgAnnotatedMutations = left_join(tsgAnnotatedMutations, excessTsgRates %>% select(gene, impact, MultiHitDriverRate, SingleHitDriverRate), by = c("gene", "impact"))
tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$geneStatus == "SingleHit", tsgAnnotatedMutations$SingleHitDriverRate, NA )
tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$geneStatus == "MultiHit", tsgAnnotatedMutations$MultiHitDriverRate, tsgAnnotatedMutations$driver )
tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$geneStatus == "Biallelic", 1, tsgAnnotatedMutations$driver)
tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$hotspot > 0, 1, tsgAnnotatedMutations$driver)
tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$redundant, 0, tsgAnnotatedMutations$driver)

tsgDrivers = tsgAnnotatedMutations %>% filter(gene %in% tsGenes$gene_name, redundant == F)
save(tsgDrivers, file = "~/hmf/RData/Processed/tsgDrivers.RData")

oncoAnnotatedMutations$geneStatus <- as.character(oncoAnnotatedMutations$geneStatus)
oncoAnnotatedMutations$geneStatus <- ifelse(oncoAnnotatedMutations$redundant, "Redundant", oncoAnnotatedMutations$geneStatus)
oncoAnnotatedMutations$impact <- as.character(oncoAnnotatedMutations$impact)
oncoAnnotatedMutations = left_join(oncoAnnotatedMutations, excessOncoRates %>% select(gene, impact, HitDriverRate), by = c("gene", "impact"))
oncoAnnotatedMutations$driver = ifelse(oncoAnnotatedMutations$geneStatus == "Hit", oncoAnnotatedMutations$HitDriverRate, NA )
oncoAnnotatedMutations$driver = ifelse(oncoAnnotatedMutations$geneStatus == "Hotspot", 1, oncoAnnotatedMutations$driver )
oncoAnnotatedMutations$driver = ifelse(oncoAnnotatedMutations$geneStatus == "NearHotspot", 1, oncoAnnotatedMutations$driver )
oncoAnnotatedMutations$driver = ifelse(oncoAnnotatedMutations$redundant, 0, oncoAnnotatedMutations$driver)

oncoDrivers = oncoAnnotatedMutations %>% filter(gene %in% oncoGenes$gene_name, redundant == F)
save(oncoDrivers, file = "~/hmf/RData/Processed/oncoDrivers.RData")

