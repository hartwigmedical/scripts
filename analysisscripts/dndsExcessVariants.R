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

load(file = "~/hmf/RData/highestPurityExonicSomaticsFilteredImpact.RData")
filteredMutations = highestPurityExonicSomaticsFilteredImpact
filteredMutations$nearHotspot <- nearHotspot(filteredMutations)

# Gene Panel
load(file = "~/hmf/RData/genePanel.RData")
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
excessOncoInframeRates = excessTsgRates %>% ungroup() %>% filter(impact == "INDEL") %>% mutate(impact = "Inframe")
excessOncoFrameshiftRates = excessTsgRates %>% ungroup() %>% filter(impact == "INDEL") %>% mutate(impact = "Frameshift")
excessOncoRates = excessOncoRates %>% filter(impact != "INDEL") %>% bind_rows(excessOncoInframeRates) %>% bind_rows(excessOncoFrameshiftRates) %>% arrange(gene, impact) %>%
  select(gene, impact, HotspotDriverRate, HitDriverRate)

save(excessTsgRates, excessOncoRates, file = "~/hmf/RData/excessRates.RData")


####### APPLY EXCESS RATES TO MUTATIONS
load("~/hmf/RData/excessRates.RData")
load(file = "~/hmf/RData/highestPurityExonicSomaticsImpact.RData")
mutations = highestPurityExonicSomaticsImpact
mutations$nearHotspot <- nearHotspot(mutations)

# Gene Panel
load(file = "~/hmf/RData/genePanel.RData")
genePanel = genePanel %>% filter(martincorena | hmf | cosmicCurated)
mutations = mutations %>% filter(gene %in% genePanel$gene_name, impact != "")

tsgAnnotatedMutations = tsg_mutations(mutations)
oncoAnnotatedMutations = onco_mutations(mutations)
load(file = "~/hmf/RData/driverGenes.RData")

tsgAnnotatedMutations$geneStatus <- ifelse(tsgAnnotatedMutations$redundant, "Redundant", tsgAnnotatedMutations$geneStatus)

tsgAnnotatedMutations = left_join(tsgAnnotatedMutations, excessTsgRates %>% select(gene, impact, MultiHitDriverRate, SingleHitDriverRate), by = c("gene", "impact"))
tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$geneStatus == "SingleHit", tsgAnnotatedMutations$SingleHitDriverRate, NA )
tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$geneStatus == "MultiHit", tsgAnnotatedMutations$MultiHitDriverRate, tsgAnnotatedMutations$driver )
tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$geneStatus == "Biallelic", 1, tsgAnnotatedMutations$driver)
tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$hotspot > 0, 1, tsgAnnotatedMutations$driver)
tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$redundant, 0, tsgAnnotatedMutations$driver)

tsgDrivers = tsgAnnotatedMutations %>% filter(gene %in% tsGenes$gene_name, driver > 0)
save(tsgDrivers, file = "~/hmf/RData/tsgDrivers.RData")

oncoAnnotatedMutations$geneStatus <- as.character(oncoAnnotatedMutations$geneStatus)
oncoAnnotatedMutations$geneStatus <- ifelse(oncoAnnotatedMutations$redundant, "Redundant", oncoAnnotatedMutations$geneStatus)
oncoAnnotatedMutations$impact <- as.character(oncoAnnotatedMutations$impact)
#oncoAnnotatedMutations$impact <- ifelse(oncoAnnotatedMutations$type == "INDEL", "INDEL", oncoAnnotatedMutations$impact)
#oncoAnnotatedMutations$impact <- ifelse(oncoAnnotatedMutations$type == "MNV", "INDEL", oncoAnnotatedMutations$impact)

oncoAnnotatedMutations = left_join(oncoAnnotatedMutations, excessOncoRates %>% select(gene, impact, HitDriverRate), by = c("gene", "impact"))

oncoAnnotatedMutations$driver = ifelse(oncoAnnotatedMutations$geneStatus == "Hit", oncoAnnotatedMutations$HitDriverRate, NA )
oncoAnnotatedMutations$driver = ifelse(oncoAnnotatedMutations$geneStatus == "Hotspot", 1, oncoAnnotatedMutations$driver )
oncoAnnotatedMutations$driver = ifelse(oncoAnnotatedMutations$geneStatus == "NearHotspot", 1, oncoAnnotatedMutations$driver )
oncoAnnotatedMutations$driver = ifelse(oncoAnnotatedMutations$redundant, 0, oncoAnnotatedMutations$driver)

oncoDrivers = oncoAnnotatedMutations %>% filter(gene %in% oncoGenes$gene_name, driver > 0)
save(oncoDrivers, file = "~/hmf/RData/oncoDrivers.RData")






load(file = "~/hmf/RData/geneCopyNumberDeletes.RData")
load(file = "~/hmf/RData/geneCopyNumberAmplifications.RData")
load(file = "~/hmf/RData/geneCopyNumberDeleteTargets.RData")
load(file = "~/hmf/RData/geneCopyNumberAmplificationTargets.RData")
load(file = "~/hmf/RData/driverGenes.RData")

#load(file = "~/hmf/RData/genePanel.RData")
#genePanel = genePanel %>% filter(hmf | martincorena | cosmicCurated)

#geneCopyNumberAmplificationTargets = geneCopyNumberAmplificationTargets %>% filter(N >= 20)

allAmplifications = geneCopyNumberAmplifications %>%
  filter(gene %in% genePanel$gene_name | gene %in% geneCopyNumberAmplificationTargets$target) %>%
  select(sampleId = sampleId, gene) %>% mutate(amp = T)

targetAmplifications = geneCopyNumberAmplifications %>%
  filter(gene %in% geneCopyNumberAmplificationTargets$target) %>%
  select(sampleId = sampleId, gene) %>% mutate(amp = T)

#geneCopyNumberDeleteTargets = geneCopyNumberDeleteTargets %>% filter(N > 5)
allDeletions = geneCopyNumberDeletes %>%
  filter(gene %in% genePanel$gene_name | gene %in% geneCopyNumberDeleteTargets$target) %>%
  select(sampleId = sampleId, gene) %>% mutate(del = T)

targetDeletions = geneCopyNumberDeletes %>%
  filter(gene %in% geneCopyNumberDeleteTargets$target) %>%
  select(sampleId = sampleId, gene) %>% mutate(del = T)

#significantCopyNumber = merge(significantAmplifications, significantDeletions, by = c("sampleId", "gene"))
load(file = "~/hmf/RData/tsgDrivers.RData")
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

load(file = "~/hmf/RData/oncoDrivers.RData")
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

driversByGene = rbind(oncoDriverByGene, tsgDriverByGene)

load(file = "~/hmf/RData/highestPurityCohort.RData")
cancerTypes = highestPurityCohort %>% select(sampleId = sampleId, primaryTumorLocation)
driversByGene = left_join(driversByGene, cancerTypes, by = "sampleId")
save(driversByGene, file = "~/hmf/RData/driversByGene.RData")

driversByGene %>% group_by(sampleId, gene) %>% filter(any(del) & any(amp)) 



load("~/hmf/RData/driversByGene.RData")
load("~/hmf/RData/cohortByPrimaryTumorLocation.RData")
load("~/hmf/RData/PrimaryTumorLocationColours.RData")

driversByGene$impact = ifelse(driversByGene$multihit, "MultiHit", driversByGene$impact)
driversByGene$impact = ifelse(!is.na(driversByGene$amp), "Amp", driversByGene$impact)
driversByGene$impact = ifelse(!is.na(driversByGene$del), "Del", driversByGene$impact)

driversByGene = driversByGene %>% left_join(cohortByPrimaryTumorLocation %>% select(primaryTumorLocation, samples = N), by = "primaryTumorLocation")
driversByGene = driversByGene %>%  mutate(driverRate = driver / samples)

tidyDriversByCancerType = driversByGene %>% group_by(primaryTumorLocation, impact, type) %>% summarise(driver = sum(driverRate))
tidyDriversByCancerTypeLevels = driversByGene %>% group_by(primaryTumorLocation) %>% summarise(driver = sum(driverRate)) %>% arrange(-driver)
tidyDriversByCancerType$primaryTumorLocation = factor(tidyDriversByCancerType$primaryTumorLocation, levels= tidyDriversByCancerTypeLevels$primaryTumorLocation)

ggplot(data=tidyDriversByCancerType, aes(x = primaryTumorLocation, y = driver)) +
  geom_bar(aes(fill = impact), stat = "identity") + 
  ggtitle("Drivers rate per cancer type") + xlab("Primary Tumor Location") + ylab("Driver Rate") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4), legend.position="bottom") + facet_wrap(~type)

tidyDriversByGene = driversByGene %>% group_by(gene, primaryTumorLocation, type, impact) %>% summarise(driver = sum(driver))
tidyDriversByGeneLevels = driversByGene %>% group_by(gene, type) %>% summarise(driver = sum(driver)) %>% arrange(-driver)
tidyDriversByGeneLevels = tidyDriversByGeneLevels[1:50, ]
tidyDriversByGene = tidyDriversByGene %>% filter(gene %in% tidyDriversByGeneLevels$gene) %>% ungroup() %>%
  mutate(gene = factor(gene, levels= tidyDriversByGeneLevels$gene)) 
  
typeColour <- ifelse(tidyDriversByGeneLevels$type == "ONCO", "red", "blue")
ggplot(data=tidyDriversByGene %>% filter(!is.na(primaryTumorLocation)), aes(x = gene, y = driver)) +
  geom_bar(aes(fill = primaryTumorLocation), stat = "identity") + 
  ggtitle("Drivers per gene") + xlab("Gene") + ylab("Drivers") +
  scale_fill_manual( values= primaryTumorLocationColours) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = typeColour), legend.position="bottom")

ggplot(data=tidyDriversByGene, aes(x = gene, y = driver)) +
  geom_bar(aes(fill = impact), stat = "identity") + 
  ggtitle("Drivers per gene") + xlab("Gene") + ylab("Drivers") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = typeColour), legend.position="bottom")




View(driversByGene %>% group_by(sampleId) %>% summarise(n = sum(driver)))


data = (driversByGene %>% mutate(CNV = !is.na(amp) | !is.na(del)) %>% group_by(sampleId, cancerType, CNV) %>% summarise(n = sum(driver)))
data = (driversByGene %>% mutate(CNV = !is.na(amp) | !is.na(del)) %>% group_by(sampleId, cancerType) %>% summarise(n = sum(driver)))

View(driversByGene %>% filter(sampleId == "CPCT02220030T"))

sum(driversByGene$driver)
sum(tsgDriverByGene$driver)

ggplot(data = data %>% filter(CNV), aes(n)) + stat_ecdf(geom="step", pad = F) 
+ facet_wrap(~cancerType)

View(data %>% filter(cancerType == "Breast"))
