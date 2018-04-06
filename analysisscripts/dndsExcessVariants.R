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
  hotspots = mutations %>% filter(hotspot > 0) %>% select(chr, pos) %>% distinct
  hrange <- GRanges(hotspots$chr, IRanges(hotspots$pos, hotspots$pos + distance))
  mrange <- GRanges(mutations$chr, IRanges(mutations$pos, mutations$pos + distance))
  
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

#load(file = "~/hmf/RData/allHighestPuritySomaticsProd.RData")
#head(minorAllelePloidy)
#minorAllelePloidy = highestPuritySomaticsProd %>% select(sampleId, chromosome, position, ref, alt, minorAllelePloidy)
#mutations = left_join(mutations, minorAllelePloidy, by = c("sampleID" = "sampleId", "chr" = "chromosome", "pos" = "position", "ref", "mut" = "alt"))
#save(mutations, file="~/hmf/RData/mutations.RData")


load(file = "~/hmf/RData/allHotspots.RData")
allHotspots$hotspot <- 1
colnames(allHotspots) <- c("chr", "pos", "ref", "mut", "hotspot")

load(file="~/hmf/RData/mutations.RData")
mutations$pid <- NULL
mutations$hotspot <- NULL
mutations = dplyr::left_join(mutations, allHotspots, by = c("chr", "pos", "ref", "mut"))
mutations$hotspot <- ifelse(is.na(mutations$hotspot), 0, 1)

mutations = mutations[!is.na(mutations$impact), ]
mutations[mutations$type == "MNP", c("impact")] <- "MNV"
mutations[mutations$type == "INDEL" & abs(nchar(mutations$ref) - nchar(mutations$mut)) %% 3 == 0 , c("impact")] <- "Inframe"
mutations[mutations$type == "INDEL" & abs(nchar(mutations$ref) - nchar(mutations$mut)) %% 3 != 0 , c("impact")] <- "Frameshift"
mutations$impact <- ifelse(mutations$impact == "Stop_loss", "Nonsense", mutations$impact)
mutations$nearHotspot <- nearHotspot(mutations)

#View(tsgAnnotatedMutations %>% filter(sampleID == "CPCT02220030T", gene == "MET"))

# Gene Panel
load(file = "~/hmf/RData/GenePanel.RData")
genePanel = genePanel %>% filter(martincorena | hmf | cosmicCurated)
mutations = mutations %>% filter(gene %in% genePanel$gene_name)

#prodDB = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
#geneCopyNumbers = purple::query_gene_copy_number_by_gene(prodDB, genePanel)
#save(geneCopyNumbers, file = '~/hmf/RData/geneCopyNumbersByGene.RData')
#dbDisconnect(prodDB)
#rm(prodDB)

#load('~/hmf/RData/geneCopyNumbersByGene.RData')
#geneCopyNumbers[, c(1:3)] <- NULL
#minGeneCopyNumbers = geneCopyNumbers %>% 
#  filter(germlineHomRegions == 0, germlineHetRegions == 0) %>%
#  mutate(relativeMinCopyNumber = minCopyNumber/ploidy) %>%
#  select(sampleID = sampleId, gene, minCopyNumber, relativeMinCopyNumber)

#minGeneCopyNumbers %>% filter(sampleID == 'DRUP01020012T')
#mutations = dplyr::left_join(mutations, minGeneCopyNumbers, by = c("sampleID", "gene"))
#mutations$del <- ifelse(!is.na(mutations$minCopyNumber) & mutations$minCopyNumber < 0.5, T, F)
#mutations$amp <-  ifelse(!is.na(mutations$relativeMinCopyNumber) & mutations$relativeMinCopyNumber >= 3, T, F)
#mutations$potentiallyBiallelic <-ifelse(!is.na(mutations$minCopyNumber) & mutations$minorAllelePloidy >  mutations$minCopyNumber, T, F)

tsgAnnotatedMutations = mutations %>%
  #filter(gene %in% c("MET", "NRAS","CDKN2A", "KRAS")) %>%
  filter(impact != "Synonymous") %>% 
  group_by(sampleID, gene) %>% 
  mutate(
    n = n(), 
    impact = factor(impact, levels = c("MNV", "Frameshift", "Nonsense", "Essential_Splice", "Missense", "Inframe")),
    geneStatus = factor(tsGeneStatus(hotspot, biallelic, n),levels = c("Biallelic", "MultiHit", "SingleHit")),
    geneStatusPrimaryPositions = tsGeneStatusPrimaryPositions(pos, geneStatus, hotspot, impact), 
    redundant = tsGeneStatusRedundant(pos, geneStatusPrimaryPositions)) %>% 
  select(-geneStatusPrimaryPositions) %>% 
  ungroup() %>%
  mutate(
    impact = as.character(impact), 
    geneStatus = as.character(geneStatus))

oncoAnnotatedMutations = mutations %>% 
#oncoAnnotatedMutations = mutations[mutations$gene %in% c("APC", "NRAS","CDKN2A", "KRAS"), ] %>%
  filter(impact %in% c("MNV", "Frameshift", "Missense", "Inframe")) %>% 
  mutate(impact = factor(impact, levels = c("Inframe", "MNV", "Missense", "Frameshift"))) %>%
  group_by(sampleID, gene) %>% 
  mutate(n = n(), geneStatus =oncoGeneStatus(hotspot, nearHotspot, n), redundant = pos != oncoGeneStatusPrimaryPosition(pos, hotspot, nearHotspot, biallelic, impact)) %>% 
  ungroup() %>%
  mutate(
    impact = as.character(impact), 
    geneStatus = as.character(geneStatus))

save(tsgAnnotatedMutations, file = "~/hmf/RData/tsgAnnotatedMutations.RData")
save(oncoAnnotatedMutations, file = "~/hmf/RData/oncoAnnotatedMutations.RData")

summarise_annotation <- function(annotatedMutations) {
  annotatedMutations$geneStatus <- ifelse(annotatedMutations$redundant, "Redundant", annotatedMutations$geneStatus)
  annotatedMutations$impact <- as.character(annotatedMutations$impact)
  #annotatedMutations$impact <- ifelse(annotatedMutations$type == "INDEL", "INDEL", annotatedMutations$impact)
  #annotatedMutations$impact <- ifelse(annotatedMutations$type == "MNV", "INDEL", annotatedMutations$impact)
  annotatedMutations = annotatedMutations %>% group_by(gene, impact, geneStatus) %>% summarise(n = n()) %>% spread(geneStatus, n)
  annotatedMutations[is.na(annotatedMutations)] <- 0
  return (annotatedMutations)
}

tsgAnnotatedMutationsByGene = summarise_annotation(tsgAnnotatedMutations)
oncoAnnotatedMutationsByGene = summarise_annotation(oncoAnnotatedMutations)

save(tsgAnnotatedMutationsByGene, file = "~/hmf/RData/tsgAnnotatedMutationsByGene.RData")
save(oncoAnnotatedMutationsByGene, file = "~/hmf/RData/oncoAnnotatedMutationsByGene.RData")

excessVariants = HmfRefCDSCv %>% 
  filter(cancerType == 'All') %>%
  select(gene = gene_name, Missense = excess_mis, Nonsense = excess_non, Essential_Splice = excess_spl, INDEL = excess_ind) %>%
  gather(impact, Excess, Missense, Nonsense, Essential_Splice, INDEL) %>%
  mutate(Excess = round(Excess))

excessTsgRates = left_join(tsgAnnotatedMutationsByGene, excessVariants, by = c("gene", "impact"))
excessTsgRates$BiallelicDriverRate = pmax(0,pmin(1.0, (excessTsgRates$Excess) / excessTsgRates$Biallelic))
excessTsgRates$MultiHitDriverRate = pmax(0,pmin(1.0, (excessTsgRates$Excess - excessTsgRates$Biallelic) / excessTsgRates$MultiHit))
excessTsgRates$SingleHitDriverRate = pmax(0,pmin(1.0, (excessTsgRates$Excess - excessTsgRates$Biallelic - excessTsgRates$MultiHit) / (excessTsgRates$SingleHit)))
excessTsgRates[is.na(excessTsgRates)] <- 0

excessOncoRates = left_join(oncoAnnotatedMutationsByGene, excessVariants, by = c("gene", "impact"))
excessOncoRates$HotspotDriverRate = pmax(0,pmin(1.0, (excessOncoRates$Excess) / (excessOncoRates$Hotspot + excessOncoRates$NearHotspot)))
excessOncoRates$HitDriverRate = pmax(0,pmin(1.0, (excessOncoRates$Excess - excessOncoRates$Hotspot - excessOncoRates$NearHotspot) / (excessOncoRates$Hit)))
excessOncoRates[is.na(excessOncoRates)] <- 0

save(excessTsgRates, excessOncoRates, file = "~/hmf/RData/ExcessRates.RData")


load("~/hmf/RData/ExcessRates.RData")
load(file = "~/hmf/RData/tsgAnnotatedMutations.RData")
load(file = "~/hmf/RData/oncoAnnotatedMutations.RData")
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
oncoAnnotatedMutations$impact <- ifelse(oncoAnnotatedMutations$type == "INDEL", "INDEL", oncoAnnotatedMutations$impact)
oncoAnnotatedMutations$impact <- ifelse(oncoAnnotatedMutations$type == "MNV", "INDEL", oncoAnnotatedMutations$impact)

oncoAnnotatedMutations = left_join(oncoAnnotatedMutations, excessOncoRates %>% select(gene, impact, HitDriverRate), by = c("gene", "impact"))

oncoAnnotatedMutations$driver = ifelse(oncoAnnotatedMutations$geneStatus == "Hit", oncoAnnotatedMutations$HitDriverRate, NA )
oncoAnnotatedMutations$driver = ifelse(oncoAnnotatedMutations$geneStatus == "Hotspot", 1, oncoAnnotatedMutations$driver )
oncoAnnotatedMutations$driver = ifelse(oncoAnnotatedMutations$geneStatus == "NearHotspot", 1, oncoAnnotatedMutations$driver )
oncoAnnotatedMutations$driver = ifelse(oncoAnnotatedMutations$redundant, 0, oncoAnnotatedMutations$driver)

oncoDrivers = oncoAnnotatedMutations %>% filter(gene %in% oncoGenes$gene_name, driver > 0)
save(oncoDrivers, file = "~/hmf/RData/oncoDrivers.RData")



load(file = "~/hmf/RData/geneCopyNumberDeletesData.RData")
load(file = "~/hmf/RData/geneCopyNumberDeleteTargets.RData")
load(file = "~/hmf/RData/geneCopyNumberAmplificationsData.RData")
load(file = "~/hmf/RData/geneCopyNumberAmplificationTargets.RData")
load(file = "~/hmf/RData/genePanel.RData")
genePanel = genePanel %>% filter(hmf | martincorena | cosmicCurated)

geneCopyNumberAmplificationTargets = geneCopyNumberAmplificationTargets %>% filter(N >= 20)
significantAmplifications = geneCopyNumberAmplifications %>%
  filter(gene %in% genePanel$gene_name | gene %in% geneCopyNumberAmplificationTargets$target) %>%
  select(sampleID = sampleId, gene) %>% mutate(amp = T)

geneCopyNumberDeleteTargets = geneCopyNumberDeleteTargets %>% filter(N > 5)
significantDeletions = geneCopyNumberDeletes %>%
  filter(gene %in% genePanel$gene_name | gene %in% geneCopyNumberDeleteTargets$target) %>%
  select(sampleID = sampleId, gene) %>% mutate(del = T)

#significantCopyNumber = merge(significantAmplifications, significantDeletions, by = c("sampleId", "gene"))


load(file = "~/hmf/RData/tsgDrivers.RData")
tsgDriverByGene = tsgDrivers %>%
  mutate(driver = ifelse(hotspot, 1, driver)) %>%
  group_by(sampleID, gene) %>%
  summarise(
    multihit = n() > 1,
    biallelic = sum(biallelic) > 0, 
    hotspot = sum(hotspot) > 0, 
    nearHotspot = F,
    impact = paste(impact, collapse = ","),
    driver = max(driver)) %>%
  full_join(significantDeletions, by = c("sampleID" ,"gene")) %>% 
  left_join(significantAmplifications, by = c("sampleID" ,"gene")) %>% 
  mutate(driver = ifelse(!is.na(del), 1, driver)) %>%
  mutate(type = 'TSG')

load(file = "~/hmf/RData/oncoDrivers.RData")
oncoDriverByGene = oncoDrivers %>%
  group_by(sampleID, gene) %>%
  summarise(
    multihit = F,
    biallelic = F,
    hotspot = sum(hotspot) > 0,
    nearHotspot = any(nearHotspot),
    impact = paste(impact, collapse = ","),
    driver = max(driver))  %>%
  full_join(significantAmplifications, by = c("sampleID" ,"gene")) %>% 
  left_join(significantDeletions, by = c("sampleID" ,"gene")) %>% 
  mutate(driver = ifelse(!is.na(amp), 1, driver)) %>%
  mutate(type = 'ONCO')

driversByGene = rbind(oncoDriverByGene, tsgDriverByGene)

cancerTypes = cohort %>% select(sampleID = sampleId, cancerType)
driversByGene = left_join(driversByGene, cancerTypes, by = "sampleID")
save(driversByGene, file = "~/hmf/RData/driversByGene.RData")


driversByGene %>% group_by(sampleID, gene) %>% filter(any(del) & any(amp)) 


load("~/hmf/RData/driversByGene.RData")
load("~/hmf/RData/cohortByCancerType.RData")
load("~/hmf/RData/cancerTypeColours.RData")

driversByGene$impact = ifelse(driversByGene$multihit, "MultiHit", driversByGene$impact)
driversByGene$impact = ifelse(!is.na(driversByGene$amp), "Amp", driversByGene$impact)
driversByGene$impact = ifelse(!is.na(driversByGene$del), "Del", driversByGene$impact)

driversByGene = driversByGene %>% left_join(cohortByCancerType %>% select(cancerType, samples = N), by = "cancerType")
driversByGene = driversByGene %>%  mutate(driverRate = driver / samples)

tidyDriversByCancerType = driversByGene %>% group_by(cancerType, impact, type) %>% summarise(driver = sum(driverRate))
tidyDriversByCancerTypeLevels = driversByGene %>% group_by(cancerType) %>% summarise(driver = sum(driverRate)) %>% arrange(-driver)
tidyDriversByCancerType$cancerType = factor(tidyDriversByCancerType$cancerType, levels= tidyDriversByCancerTypeLevels$cancerType)

ggplot(data=tidyDriversByCancerType, aes(x = cancerType, y = driver)) +
  geom_bar(aes(fill = impact), stat = "identity") + 
  ggtitle("Drivers rate per cancer type") + xlab("Cancer Type") + ylab("Driver Rate") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4), legend.position="bottom") + facet_wrap(~type)

tidyDriversByGene = driversByGene %>% group_by(gene, cancerType, type, impact) %>% summarise(driver = sum(driver))
tidyDriversByGeneLevels = driversByGene %>% group_by(gene, type) %>% summarise(driver = sum(driver)) %>% arrange(-driver)
tidyDriversByGeneLevels = tidyDriversByGeneLevels[1:50, ]
tidyDriversByGene = tidyDriversByGene %>% filter(gene %in% tidyDriversByGeneLevels$gene) %>% ungroup() %>%
  mutate(gene = factor(gene, levels= tidyDriversByGeneLevels$gene)) 
  
typeColour <- ifelse(tidyDriversByGeneLevels$type == "ONCO", "red", "blue")
ggplot(data=tidyDriversByGene %>% filter(!is.na(cancerType)), aes(x = gene, y = driver)) +
  geom_bar(aes(fill = cancerType), stat = "identity") + 
  ggtitle("Drivers per gene") + xlab("Gene") + ylab("Drivers") +
  scale_fill_manual( values= cancerTypeColours) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = typeColour), legend.position="bottom")

ggplot(data=tidyDriversByGene, aes(x = gene, y = driver)) +
  geom_bar(aes(fill = impact), stat = "identity") + 
  ggtitle("Drivers per gene") + xlab("Gene") + ylab("Drivers") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = typeColour), legend.position="bottom")



View(driversByGene %>% group_by(sampleID) %>% summarise(n = sum(driver)))


data = (driversByGene %>% mutate(CNV = !is.na(amp) | !is.na(del)) %>% group_by(sampleID, cancerType, CNV) %>% summarise(n = sum(driver)))
data = (driversByGene %>% mutate(CNV = !is.na(amp) | !is.na(del)) %>% group_by(sampleID, cancerType) %>% summarise(n = sum(driver)))

View(driversByGene %>% filter(sampleID == "CPCT02220030T"))

sum(driversByGene$driver)
sum(tsgDriverByGene$driver)

ggplot(data = data %>% filter(CNV), aes(n)) + stat_ecdf(geom="step", pad = F) 
+ facet_wrap(~cancerType)

View(data %>% filter(cancerType == "Breast"))
