library(IRanges)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)


#Excess
load("~/hmf/RData/HmfRefCDSCv.RData")
HmfRefCDSCvNullPcawg$gene_name <- as.character(HmfRefCDSCvNullPcawg$gene_name)
HmfRefCDSCv$prob_mis = ifelse(HmfRefCDSCv$n_mis>0,pmax(0,(HmfRefCDSCv$wmis_cv-1)/HmfRefCDSCv$wmis_cv),0)
HmfRefCDSCv$prob_non = ifelse(HmfRefCDSCv$n_non,pmax(0,(HmfRefCDSCv$wnon_cv-1)/HmfRefCDSCv$wnon_cv),0)
HmfRefCDSCv$prob_spl = ifelse(HmfRefCDSCv$n_spl>0,pmax(0,(HmfRefCDSCv$wspl_cv-1)/HmfRefCDSCv$wspl_cv),0)
HmfRefCDSCv$prob_ind = ifelse(HmfRefCDSCv$n_ind>0,pmax(0,(HmfRefCDSCv$wind_cv-1)/HmfRefCDSCv$wind_cv),0)
HmfRefCDSCv$excess_mis = HmfRefCDSCv$prob_mis*HmfRefCDSCv$n_mis
HmfRefCDSCv$excess_non = HmfRefCDSCv$prob_non*HmfRefCDSCv$n_non
HmfRefCDSCv$excess_spl = HmfRefCDSCv$prob_spl*HmfRefCDSCv$n_spl
HmfRefCDSCv$excess_ind = HmfRefCDSCv$prob_ind*HmfRefCDSCv$n_ind






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

# Copy Number
load(file = "~/hmf/RData/geneCopyNumberAmplificationsData.RData")

load(file = "~/hmf/RData/geneCopyNumberDeletesData.RData")
geneCopyNumberDeletes = geneCopyNumberDeletes %>% select(sampleID = sampleId, gene) %>% mutate(del = T)
#mutations = dplyr::left_join(mutations, geneCopyNumberDeletes, by = c("sampleID", "gene"))

tsgByVariant = mutations %>% 
#tsgByVariant = mutations[mutations$gene %in% c("APC", "NRAS","CDKN2A", "KRAS"), ] %>%
  filter(impact != "Synonymous") %>% 
  mutate(impact = factor(impact, levels = c("MNV", "Frameshift", "Nonsense", "Essential_Splice", "Missense", "Inframe"))) %>%
  group_by(sampleID, gene) %>% 
  mutate(n = n(), geneStatus = tsGeneStatus(biallelic, n), geneStatusPrimaryPositions = tsGeneStatusPrimaryPositions(pos, geneStatus, impact), redundant = tsGeneStatusRedundant(pos, geneStatusPrimaryPositions)) %>% 
  select(-geneStatusPrimaryPositions) %>% ungroup()


oncoByVariant = mutations %>% 
#oncoByVariant = mutations[mutations$gene %in% c("APC", "NRAS","CDKN2A", "KRAS"), ] %>%
  filter(impact %in% c("MNV", "Frameshift", "Missense", "Inframe")) %>% 
  mutate(impact = factor(impact, levels = c("MNV", "Frameshift", "Missense", "Inframe"))) %>%
  group_by(sampleID, gene) %>% 
  mutate(n = n(), geneStatus =oncoGeneStatus(hotspot, nearHotspot, n), redundant = pos != oncoGeneStatusPrimaryPosition(pos, hotspot, nearHotspot, impact)) %>% 
  ungroup()

save(tsgByVariant, file = "~/hmf/RData/tsgByVariantNewHotspots.RData")
save(oncoByVariant, file = "~/hmf/RData/oncoByVariantNewHotspots.RData")

View(tsgByVariant %>% unite(combined,c(driver,redundant, impact)) %>% group_by(gene,combined) %>% 
       summarise(count=n()) %>% 
       spread(combined,count) )




# Signficant Genes
library(RMySQL)
pilotDB = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
genePanel = (dbGetQuery(pilotDB, "SELECT distinct gene as gene_name FROM genePanel"))
dbDisconnect(pilotDB)
rm(pilotDB)

load("~/hmf/RData/HmfRefCDSCv.RData")
load("~/hmf/RData/PcawgRefCDSCv.RData")
HMFSigGenes = HmfRefCDSCv %>% filter(qglobal_cv < 0.05) %>% group_by(gene_name) %>% summarise() %>% as.data.frame
PCAWGSigGenes = PcawgRefCDSCv %>% filter(qglobal < 0.05) %>% group_by(gene_name) %>% summarise() %>% as.data.frame
genes = rbind(HMFSigGenes, PCAWGSigGenes)
genes = rbind(genes, genePanel)
genes = unique(genes)
load(file = "~/hmf/RData/knownCancerGenes.RData")
colnames(knownCancerGenes) <- c("gene_name","oncogene","tsg")
GenePanel = merge(genes, knownCancerGenes, all.x = T)
save(GenePanel, file = "~/hmf/RData/GenePanel.RData")

load("~/hmf/RData/knownCancerGenes.RData")
colnames(knownCancerGenes) <- c("gene_name", "oncogene","tsg")
knownCancerGenes$cosmic_type <- ifelse (knownCancerGenes$oncogene, "oncogene", "neither")
knownCancerGenes$cosmic_type <- ifelse (knownCancerGenes$tsg, "tsg", knownCancerGenes$cosmic_type)
knownCancerGenes$cosmic_type <- ifelse (knownCancerGenes$tsg & knownCancerGenes$oncogene, "both", knownCancerGenes$cosmic_type)
save(knownCancerGenes, file = "~/hmf/RData/knownCancerGenes.RData")




