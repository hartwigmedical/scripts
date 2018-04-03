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

tsGeneStatus <-function(del, biallelic, n) {
  result = ifelse(n > 1, "MultiHit", "SingleHit")
  result = ifelse(biallelic, "Biallelic", result)
  result = ifelse(del, "Del", result)
  return (result)
}

tsGeneStatusPrimaryPositions <- function(pos, geneStatus, impact) {
  df = data.frame(pos = pos, geneStatus = geneStatus, impact = impact)
  del = df %>% filter(geneStatus == "Del") %>% arrange(impact)
  if (nrow(del)) {
    return (as.character(del[1, c("pos")]))
  }
  
  biallelic = df %>% filter(geneStatus == "Biallelic") %>% arrange(impact)
  if (nrow(biallelic)) {
    return (as.character(biallelic[1, c("pos")]))
  }
  
  multihit = df %>% filter(geneStatus == "MultiHit") %>% arrange(impact)
  if (nrow(multihit) > 1) {
    return (paste(multihit[1:2, c("pos")], collapse =","))
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

oncoGeneStatus <-function(hotspot, nearHotspot, amp, n) {
  result = ifelse(hotspot, "Hotspot", "Hit")
  result = ifelse(nearHotspot, "NearHotspot", result)
  result = ifelse(amp, "Amp", result)
  return (result)
}

oncoGeneStatusPrimaryPosition <- function(pos, hotspot, nearHotspot, amp, impact) {
  df = data.frame(pos = pos, hotspot = hotspot, nearHotspot = nearHotspot, impact = impact)
  ordered = df %>% arrange(-hotspot, -nearHotspot, -amp, impact)
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

# Gene Panel
load(file = "~/hmf/RData/GenePanel.RData")
genePanel = genePanel %>% filter(martincorena | hmf | cosmicCurated)
mutations = mutations %>% filter(gene %in% genePanel$gene_name)

#prodDB = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
#geneCopyNumbers = purple::query_gene_copy_number_by_gene(prodDB, genePanel)
#save(geneCopyNumbers, file = '~/hmf/RData/geneCopyNumbersByGene.RData')
#dbDisconnect(prodDB)
#rm(prodDB)

load('~/hmf/RData/geneCopyNumbersByGene.RData')
geneCopyNumbers[, c(1:3)] <- NULL
minGeneCopyNumbers = geneCopyNumbers %>% 
  filter(germlineHomRegions == 0, germlineHetRegions == 0) %>%
  mutate(relativeMinCopyNumber = minCopyNumber/ploidy) %>%
  select(sampleID = sampleId, gene, minCopyNumber, relativeMinCopyNumber)


#minGeneCopyNumbers %>% filter(sampleID == 'DRUP01020012T')
mutations = dplyr::left_join(mutations, minGeneCopyNumbers, by = c("sampleID", "gene"))
mutations$del <- ifelse(!is.na(mutations$minCopyNumber) & mutations$minCopyNumber < 0.5, T, F)
mutations$amp <-  ifelse(!is.na(mutations$relativeMinCopyNumber) & mutations$relativeMinCopyNumber >= 3, T, F)
mutations$potentiallyBiallelic <-ifelse(!is.na(mutations$minCopyNumber) & mutations$minorAllelePloidy >  mutations$minCopyNumber, T, F)


tsgAnnotatedMutations = mutations %>% 
#tsgAnnotatedMutations = mutations[mutations$gene %in% c("APC", "NRAS","CDKN2A", "KRAS"), ] %>%
  filter(impact != "Synonymous") %>% 
  mutate(impact = factor(impact, levels = c("MNV", "Frameshift", "Nonsense", "Essential_Splice", "Missense", "Inframe"))) %>%
  group_by(sampleID, gene) %>% 
  mutate(n = n(), geneStatus = tsGeneStatus(del, biallelic, n), geneStatusPrimaryPositions = tsGeneStatusPrimaryPositions(pos, geneStatus, impact), redundant = tsGeneStatusRedundant(pos, geneStatusPrimaryPositions)) %>% 
  select(-geneStatusPrimaryPositions) %>% ungroup()

oncoAnnotatedMutations = mutations %>% 
#oncoAnnotatedMutations = mutations[mutations$gene %in% c("APC", "NRAS","CDKN2A", "KRAS"), ] %>%
  filter(impact %in% c("MNV", "Frameshift", "Missense", "Inframe")) %>% 
  mutate(impact = factor(impact, levels = c("MNV", "Frameshift", "Missense", "Inframe"))) %>%
  group_by(sampleID, gene) %>% 
  mutate(n = n(), geneStatus =oncoGeneStatus(hotspot, nearHotspot, amp, n), redundant = pos != oncoGeneStatusPrimaryPosition(pos, hotspot, nearHotspot, amp, impact)) %>% 
  ungroup()

save(tsgAnnotatedMutations, file = "~/hmf/RData/tsgAnnotatedMutations.RData")
save(oncoAnnotatedMutations, file = "~/hmf/RData/oncoAnnotatedMutations.RData")

summarise_annotation <- function(annotatedMutations) {
  annotatedMutations$geneStatus <- ifelse(annotatedMutations$redundant, "Redundant", annotatedMutations$geneStatus)
  annotatedMutations$impact <- as.character(annotatedMutations$impact)
  annotatedMutations$impact <- ifelse(annotatedMutations$type == "INDEL", "INDEL", annotatedMutations$impact)
  annotatedMutations$impact <- ifelse(annotatedMutations$type == "MNV", "INDEL", annotatedMutations$impact)
  annotatedMutations = annotatedMutations %>% group_by(gene, impact, geneStatus) %>% summarise(n = n()) %>% spread(geneStatus, n)
  annotatedMutations[is.na(annotatedMutations)] <- 0
  return (annotatedMutations)
}

tsgAnnotatedMutationsByGene = summarise_annotation(tsgAnnotatedMutations)
oncoAnnotatedMutationsByGene = summarise_annotation(oncoAnnotatedMutations)

save(tsgAnnotatedMutationsByGene, file = "~/hmf/RData/tsgAnnotatedMutationsByGene.RData")
save(oncoAnnotatedMutationsByGene, file = "~/hmf/RData/oncoAnnotatedMutationsByGene.RData")

