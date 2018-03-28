library(IRanges)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)

nearHotspot <-function(mutations, distance = 10) { 
  hotspots = mutations %>% filter(hotspot > 0) %>% select(chr, pos) %>% distinct
  hrange <- GRanges(hotspots$chr, IRanges(hotspots$pos, hotspots$pos + distance))
  mrange <- GRanges(mutations$chr, IRanges(mutations$pos, mutations$pos + distance))
  
  ol = as.matrix(findOverlaps(hrange, mrange, type="any", select="all"))
  mutations$nearHotspot <- FALSE
  mutations[ol[,2], c("nearHotspot")] <- TRUE
  return (mutations$nearHotspot & !mutations$hotspot)
}

tsGeneStatus <-function(biallelic, n) {
  result = ifelse(n > 1, "MultiHit", "SingleHit")
  result = ifelse(biallelic, "Biallelic", result)
  return (result)
}

tsGeneStatusPrimaryPositions <- function(pos, geneStatus, impact) {
  df = data.frame(pos = pos, geneStatus = geneStatus, impact = impact)
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

oncoGeneStatus <-function(hotspot, nearHotspot, n) {
  result = ifelse(n > 1, "MultiHit", "SingleHit")
  result = ifelse(hotspot, "Hotspot", result)
  result = ifelse(nearHotspot, "NearHotspot", result)
  return (result)
}

oncoGeneStatusPrimaryPosition <- function(pos, hotspot, nearHotspot, impact) {
  df = data.frame(pos = pos, hotspot = hotspot, nearHotspot = nearHotspot, impact = impact)
  ordered = df %>% arrange(-hotspot, -nearHotspot, impact)
  if (nrow(ordered) > 0) {
    return (ordered[1, c("pos")])
  }
  
  return (NA)
}

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


tsgSomatics = mutations %>% 
#tsgByVariant = mutations[mutations$gene %in% c("APC", "NRAS","CDKN2A", "KRAS"), ] %>%
  filter(impact != "Synonymous") %>% 
  mutate(impact = factor(impact, levels = c("MNV", "Frameshift", "Nonsense", "Essential_Splice", "Missense", "Inframe"))) %>%
  group_by(sampleID, gene) %>% 
  mutate(n = n(), geneStatus = tsGeneStatus(biallelic, n), geneStatusPrimaryPositions = tsGeneStatusPrimaryPositions(pos, geneStatus, impact), redundant = tsGeneStatusRedundant(pos, geneStatusPrimaryPositions)) %>% 
  select(-geneStatusPrimaryPositions) %>% ungroup()


oncoSomatics = mutations %>% 
#oncoByVariant = mutations[mutations$gene %in% c("APC", "NRAS","CDKN2A", "KRAS"), ] %>%
  filter(impact %in% c("MNV", "Frameshift", "Missense", "Inframe")) %>% 
  mutate(impact = factor(impact, levels = c("MNV", "Frameshift", "Missense", "Inframe"))) %>%
  group_by(sampleID, gene) %>% 
  mutate(n = n(), geneStatus =oncoGeneStatus(hotspot, nearHotspot, n), redundant = pos != oncoGeneStatusPrimaryPosition(pos, hotspot, nearHotspot, impact)) %>% 
  ungroup()

save(tsgSomatics, file = "~/hmf/RData/tsgSomatics.RData")
save(oncoSomatics, file = "~/hmf/RData/oncoSomatics.RData")

