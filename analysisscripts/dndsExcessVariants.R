library(IRanges)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
detach("package:purple", unload=TRUE);
library(purple)
library(RMySQL)

####### Variant Annotation
nearHotspot <-function(mutations, distance = 10) { 
  hotspots = mutations %>% filter(hotspot > 0) %>% select(chromosome, position) %>% distinct
  hrange <- GRanges(hotspots$chromosome, IRanges(hotspots$position, hotspots$position + distance))
  mrange <- GRanges(mutations$chromosome, IRanges(mutations$position, mutations$position + distance))
  
  ol = as.matrix(findOverlaps(hrange, mrange, type="any", select="all"))
  mutations$nearHotspot <- FALSE
  mutations[ol[,2], c("nearHotspot")] <- TRUE
  return (mutations$nearHotspot & !mutations$hotspot)
}

tsGeneStatus <-function(hotspot, biallelic, n) {
  levels = c("Biallelic", "MultiHit", "SingleHit")
    
  result = ifelse(n > 1, factor("MultiHit", levels) , factor("SingleHit", levels))
  result = ifelse(biallelic, factor("Biallelic", levels), result)
  return (result)
}

tsGeneStatus <-function(hotspot, biallelic, n) {
  result = ifelse(n > 1, "MultiHit" , "SingleHit")
  result = ifelse(biallelic, "Biallelic", result)
  return (result)
}

tsGeneStatusPrimaryPositions <- function(pos, geneStatus, hotspot, impact) {
  df = data.frame(pos = pos, geneStatus = geneStatus, hotspot = hotspot, impact = impact) %>% 
    arrange(geneStatus, -hotspot, impact)
  
  if (df[1, "geneStatus"] == "MultiHit") {
    return (paste(df[1:2, c("pos")], collapse =","))
  }
  
  return (as.character(df[1, c("pos")]))
}

tsGeneStatusRedundant <- function(position, driverPositions) {
  result = vector(mode = "character", length = length(position))
  for (i in 1:length(position)) {
    drivers = as.integer(strsplit(driverPositions[i], ",")[[1]])
    if (position[i] %in% drivers) {
      result[i] <- FALSE
    } else {
      result[i] <- TRUE
    }
  }
  return (result)
}

oncoGeneStatus <-function(hotspot, nearHotspot, n) {
  result = ifelse(nearHotspot, "NearHotspot", "Hit")
  result = ifelse(hotspot, "Hotspot", result)
  return (result)
}

oncoGeneStatusPrimaryPosition <- function(pos, hotspot, nearHotspot, biallelic, impact) {
  df = data.frame(pos = pos, hotspot = hotspot, nearHotspot = nearHotspot, biallelic = biallelic, impact = impact)
  ordered = df %>% arrange(-hotspot, -nearHotspot, -biallelic, impact)
  if (nrow(ordered) > 0) {
    return (ordered[1, c("pos")])
  }
  
  return (NA)
}

tsg_mutations <- function(mutations) {
  result = mutations %>%
    #filter(gene %in% c("MET", "NRAS","CDKN2A", "KRAS")) %>%
    filter(impact != "Synonymous") %>% 
    group_by(sampleId, gene) %>% 
    mutate(
      n = n(), 
      impact = factor(impact, levels = c("MNV", "Frameshift", "Nonsense", "Splice", "Missense", "Inframe")),
      geneStatus = factor(tsGeneStatus(hotspot, biallelic, n),levels = c("Biallelic", "MultiHit", "SingleHit")),
      geneStatusPrimaryPositions = tsGeneStatusPrimaryPositions(position, geneStatus, hotspot, impact),
      redundant = tsGeneStatusRedundant(position, geneStatusPrimaryPositions)) %>%
    select(-geneStatusPrimaryPositions) %>% 
    ungroup() %>%
    mutate(
      impact = as.character(impact), 
      geneStatus = as.character(geneStatus))
  
  return (result)
}

onco_mutations <- function(mutations) {
  result = mutations %>% 
    #oncoAnnotatedMutations = mutations[mutations$gene %in% c("APC", "NRAS","CDKN2A", "KRAS"), ] %>%
    filter(impact %in% c("MNV", "Frameshift", "Missense", "Inframe")) %>% 
    mutate(impact = factor(impact, levels = c("Inframe", "MNV", "Missense", "Frameshift"))) %>%
    group_by(sampleId, gene) %>%
    mutate(n = n(), geneStatus =oncoGeneStatus(hotspot, nearHotspot, n), redundant = position != oncoGeneStatusPrimaryPosition(position, hotspot, nearHotspot, biallelic, impact)) %>%
    ungroup() %>%
    mutate(
      impact = as.character(impact), 
      geneStatus = as.character(geneStatus))
  
  return (result)
}

annotate_somatics <- function(annotmuts, somatics) {
  result = left_join(annotmuts, somatics, by = c("sampleID","chr","pos","ref", "mut")) %>% 
    select(sampleId = sampleID, chromosome = chr, position = pos, ref = ref, alt = mut, gene, type, impact, worstCodingEffect, canonicalCodingEffect, hotspot, biallelic, clonality)
  
  result = result %>%
    mutate(impact = ifelse(impact == "Essential_Splice", "Splice", impact)) %>%
    mutate(impact = ifelse(impact == "no-SNV", paste0(substr(canonicalCodingEffect, 1, 1), tolower(substring(canonicalCodingEffect, 2))), impact)) %>%
    mutate(impact = ifelse(impact == "None", paste0(substr(worstCodingEffect, 1, 1), tolower(substring(worstCodingEffect, 2))), impact)) %>%
    mutate(impact = ifelse(impact %in% c("Stop_loss","Nonsense_or_frameshift"), "Nonsense", impact)) %>%
    mutate(impact = ifelse(type == "INDEL" & impact == "Missense", "Inframe", impact)) %>%
    mutate(impact = ifelse(type == "INDEL" & impact == "Nonsense", "Frameshift", impact)) %>%
    filter(!is.na(impact), impact != "None")
  
  return (result) 
}

load(file = "~/hmf/RData/reference/hpcExonicSomatics.RData")
somatics = hpcExonicSomatics %>% 
  select(sampleID = sampleId, chr = chromosome, pos = position, ref = ref, mut = alt, type, worstCodingEffect, canonicalCodingEffect, hotspot, biallelic, clonality)

load(file = "~/hmf/RData/processed/dndsFilteredAnnotatedMutations.RData")
filteredMutations = annotate_somatics(dndsFilteredAnnotatedMutations, somatics)
filteredMutations$nearHotspot <- nearHotspot(filteredMutations)

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
  select(sampleID = sampleId, chr = chromosome, pos = position, ref = ref, mut = alt, type, worstCodingEffect, canonicalCodingEffect, hotspot, biallelic, clonality)

load(file = "~/hmf/RData/processed/dndsUnfilteredAnnotatedMutations.RData")
mutations = annotate_somatics(dndsUnfilteredAnnotatedMutations, somatics)
mutations$nearHotspot <- nearHotspot(mutations)

# Gene Panel
load(file = "~/hmf/RData/processed/genePanel.RData")
genePanel = genePanel %>% filter(martincorena | hmf | cosmicCurated)
mutations = mutations %>% filter(gene %in% genePanel$gene_name, impact != "")

tsgAnnotatedMutations = tsg_mutations(mutations)
oncoAnnotatedMutations = onco_mutations(mutations)

load(file = "~/hmf/RData/processed/driverGenes.RData")
tsgAnnotatedMutations$geneStatus <- ifelse(tsgAnnotatedMutations$redundant, "Redundant", tsgAnnotatedMutations$geneStatus)
tsgAnnotatedMutations = left_join(tsgAnnotatedMutations, excessTsgRates %>% select(gene, impact, MultiHitDriverRate, SingleHitDriverRate), by = c("gene", "impact"))
tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$geneStatus == "SingleHit", tsgAnnotatedMutations$SingleHitDriverRate, NA )
tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$geneStatus == "MultiHit", tsgAnnotatedMutations$MultiHitDriverRate, tsgAnnotatedMutations$driver )
tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$geneStatus == "Biallelic", 1, tsgAnnotatedMutations$driver)
tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$hotspot > 0, 1, tsgAnnotatedMutations$driver)
tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$redundant, 0, tsgAnnotatedMutations$driver)

tsgDrivers = tsgAnnotatedMutations %>% filter(gene %in% tsGenes$gene_name, driver > 0)
save(tsgDrivers, file = "~/hmf/RData/processed/tsgDrivers.RData")

oncoAnnotatedMutations$geneStatus <- as.character(oncoAnnotatedMutations$geneStatus)
oncoAnnotatedMutations$geneStatus <- ifelse(oncoAnnotatedMutations$redundant, "Redundant", oncoAnnotatedMutations$geneStatus)
oncoAnnotatedMutations$impact <- as.character(oncoAnnotatedMutations$impact)
oncoAnnotatedMutations = left_join(oncoAnnotatedMutations, excessOncoRates %>% select(gene, impact, HitDriverRate), by = c("gene", "impact"))
oncoAnnotatedMutations$driver = ifelse(oncoAnnotatedMutations$geneStatus == "Hit", oncoAnnotatedMutations$HitDriverRate, NA )
oncoAnnotatedMutations$driver = ifelse(oncoAnnotatedMutations$geneStatus == "Hotspot", 1, oncoAnnotatedMutations$driver )
oncoAnnotatedMutations$driver = ifelse(oncoAnnotatedMutations$geneStatus == "NearHotspot", 1, oncoAnnotatedMutations$driver )
oncoAnnotatedMutations$driver = ifelse(oncoAnnotatedMutations$redundant, 0, oncoAnnotatedMutations$driver)

oncoDrivers = oncoAnnotatedMutations %>% filter(gene %in% oncoGenes$gene_name, driver > 0)
save(oncoDrivers, file = "~/hmf/RData/processed/oncoDrivers.RData")