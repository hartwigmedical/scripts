library(IRanges)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)

load("~/hmf/RData/HmfRefCDSCv.RData")
HmfRefCDSCv$probmis = ifelse(HmfRefCDSCv$n_mis>0,pmax(0,(HmfRefCDSCv$wmis_cv-1)/HmfRefCDSCv$wmis_cv),0)
HmfRefCDSCv$probnon = ifelse(HmfRefCDSCv$n_non,pmax(0,(HmfRefCDSCv$wnon_cv-1)/HmfRefCDSCv$wnon_cv),0)
HmfRefCDSCv$probspl = ifelse(HmfRefCDSCv$n_spl>0,pmax(0,(HmfRefCDSCv$wspl_cv-1)/HmfRefCDSCv$wspl_cv),0)
HmfRefCDSCv$probind = ifelse(HmfRefCDSCv$n_ind>0,pmax(0,(HmfRefCDSCv$wind_cv-1)/HmfRefCDSCv$wind_cv),0)
HmfRefCDSCv$excess_mis = HmfRefCDSCv$probmis*HmfRefCDSCv$n_mis
HmfRefCDSCv$excess_non = HmfRefCDSCv$probnon*HmfRefCDSCv$n_non
HmfRefCDSCv$excess_spl = HmfRefCDSCv$probspl*HmfRefCDSCv$n_spl
HmfRefCDSCv$excess_ind = HmfRefCDSCv$probind*HmfRefCDSCv$n_ind


load("~/hmf/RData/geneClassification.RData")
load("~/hmf/RData/oncoByVariantNewHotspots.RData")
load("~/hmf/RData/tsgByVariantNewHotspots.RData")
load("~/hmf/RData/canonicalTranscripts.RData")
canonicalTranscripts$codons = round(canonicalTranscripts$codingBases /3)


oncoByVariant$geneStatus <- ifelse(oncoByVariant$redundant, "Redundant", oncoByVariant$geneStatus)
oncoByVariant$impact <- as.character(oncoByVariant$impact)
oncoByVariant$impact <- ifelse(oncoByVariant$type == "INDEL", "INDEL", oncoByVariant$impact)
oncoByVariant$impact <- ifelse(oncoByVariant$type == "MNV", "INDEL", oncoByVariant$impact)
oncoByVariantSummary = oncoByVariant %>% filter(gene %in% oncoGenes$gene_name) %>% group_by(gene, impact, geneStatus) %>% summarise(n = n()) %>% spread(geneStatus, n)
oncoByVariantSummary[is.na(oncoByVariantSummary)] <- 0

oncoExcess = HmfRefCDSCv %>% filter(gene_name %in% oncoGenes$gene_name, cancerType == "All")  

oncoMissense = oncoExcess  %>% select(gene = gene_name, excess_mis) %>% filter(excess_mis > 0)
oncoMissense = left_join(oncoMissense, oncoByVariantSummary, by = "gene" ) %>% filter(impact == "Missense") %>% arrange(-excess_mis)
oncoMissense$Hit <- oncoMissense$MultiHit + oncoMissense$SingleHit
oncoMissense$impact <- NULL
oncoMissense$excess_mis <- round(oncoMissense$excess_mis)
oncoMissense$Unexplained <-  oncoMissense$excess_mis - oncoMissense$Hotspot - oncoMissense$NearHotspot
oncoMissense$RequiredHits <- pmax(0, oncoMissense$Unexplained)
oncoMissense$InsufficientHits <- pmin(0, oncoMissense$Hit - oncoMissense$RequiredHits)
oncoMissense$ExcessHits <- pmax(0, oncoMissense$Hit - oncoMissense$RequiredHits)




###### TUMOR SUPRESSOR

tsgByVariant$geneStatus <- ifelse(tsgByVariant$redundant, "Redundant", tsgByVariant$geneStatus)
tsgByVariant$impact <- as.character(tsgByVariant$impact)
tsgByVariant$impact <- ifelse(tsgByVariant$type == "INDEL", "INDEL", tsgByVariant$impact)
tsgByVariant$impact <- ifelse(tsgByVariant$type == "MNV", "INDEL", tsgByVariant$impact)

#tsgByVariant$impact <- ifelse(tsgByVariant$impact == "Frameshift", "INDEL", tsgByVariant$impact)
#tsgByVariant$impact <- ifelse(tsgByVariant$impact == "Inframe", "INDEL", tsgByVariant$impact)

tsgSummary = tsgByVariant %>% filter(gene %in% tsGenes$gene_name) %>% group_by(gene, impact, geneStatus) %>% summarise(n = n()) %>% spread(geneStatus, n)
tsgSummary[is.na(tsgSummary)] <- 0
tsgExcess = HmfRefCDSCv %>% filter(gene_name %in% tsGenes$gene_name, cancerType == "All")  

#missense = tsgExcess %>% select(gene = gene_name, excess_mis)
#missense = left_join(missense, tsgSummary, by = "gene" ) %>% filter(impact == "Missense") %>% arrange(-excess_mis)
##missense$impact <- NULL
#missense$excess_mis <- round(missense$excess_mis)
#missense$Unexplained <-  missense$excess_mis - missense$Biallelic - missense$MultiHit
#missense$RequiredSingles <- pmax(0, missense$Unexplained)
#missense$InsufficientSingles <- pmin(0, missense$SingleHit - missense$RequiredSingles)
#missense$ExcessSingles <- pmax(0, missense$SingleHit - missense$RequiredSingles)

#identical(missense, jon)

combine_tsg_with_excess <- function(tsgCounts, excessSummary) {
  tsgCounts$gene <- as.character(tsgCounts$gene)
  excessSummary$gene <- as.character(excessSummary$gene)

  combined = left_join(excessSummary, tsgCounts, by = "gene" ) %>% arrange(-Excess)
  combined$Excess <- round(combined$Excess)
  combined = combined %>% filter(Excess > 0)
  combined$Unexplained <-  combined$Excess - combined$Biallelic - combined$MultiHit
  combined$RequiredSingles <- pmax(0, combined$Unexplained)
  combined$InsufficientSingles <- pmin(0, combined$SingleHit - combined$RequiredSingles)
  combined$ExcessSingles <- pmax(0, combined$SingleHit - combined$RequiredSingles)
  
  combined = left_join(combined, canonicalTranscripts[, c("gene","codons")], by = "gene")

  return (combined)
}

combined = missense
combine_lengths <- function(combined) {
  combined$N = combined$Biallelic + combined$MultiHit + combined$Redundant + combined$SingleHit
  combined$BiallelicOrMultihit = (combined$Biallelic + combined$MultiHit)
  
  result = combined %>% select(gene, codons, ExcessSingles, MultiHit, BiallelicOrMultihit) %>% arrange(-codons)
  #result$gene = factor(result$gene, levels = result$gene)
  
  
  #result = result %>% gather(type, N, ExcessSingles, MultiHit)
  return(result)
}

missense = combine_tsg_with_excess(tsgSummary %>% filter(impact == "Missense"), tsgExcess %>% select(gene = gene_name, Excess = excess_mis))
nonsense = combine_tsg_with_excess(tsgSummary %>% filter(impact == "Nonsense"), tsgExcess %>% select(gene = gene_name, Excess = excess_non))
splice = combine_tsg_with_excess(tsgSummary %>% filter(impact == "Essential_Splice"), tsgExcess %>% select(gene = gene_name, Excess = excess_spl))


tsgMissenseLengths = tsgMissenseSummary %>% select(gene, codons, ExcessSingles, MultiHit) %>% arrange(-codons)  %>% gather(type, N, ExcessSingles, MultiHit)

missenseLengths = combine_lengths(missense)
nonsenseLengths = combine_lengths(nonsense)
spliceLengths = combine_lengths(splice)

ggplot(data=missenseLengths %>% filter (gene != "TP53" & gene != "FAT4"), aes(MultiHit, ExcessSingles, label = gene)) +
  geom_point(stat = "identity") + 
  ggtitle(paste0("TSG - " ,"Missense")) + xlab("Multihit") + ylab("ExcessSingles") +geom_text(size=2,hjust = 0, nudge_x = 0.15)

ggplot(data=missenseLengths %>% filter (gene != "TP53" & gene != "FAT4"), aes(codons, N)) +
  geom_point(aes(color = type), stat = "identity") + 
  ggtitle(paste0("TSG - " ,"Missense")) + xlab("Codons") + ylab("")

ggplot(data=missenseLengths, aes(codons, N)) +
  geom_point(aes(color = type), stat = "identity") + 
  ggtitle(paste0("TSG - " ,"Missense")) + xlab("Codons") + ylab("")

ggplot(data=tsgMissenseLengths, aes(codons, N)) +
  geom_point(aes(color = type), stat = "identity") + 
  ggtitle(paste0("TSG - " ,"Missense")) + xlab("Codons") + ylab("")

ggplot(data=nonsenseLengths, aes(codons, N)) +
  geom_point(aes(color = type), stat = "identity") + 
  ggtitle(paste0("TSG - " ,"Nonsense")) + xlab("Codons") + ylab("")



ggplot(data=nonsenseLengths, aes(gene, N)) +
  geom_bar(stat = "identity") + 
  ggtitle(paste0("TSG - " ,"Nonsense")) + xlab("Genes") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(type ~ ., scales = "free_y")

ggplot(data=spliceLengths, aes(gene, N)) +
  geom_bar(stat = "identity") + 
  ggtitle(paste0("TSG - " ,"Splice")) + xlab("Genes") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(type ~ ., scales = "free_y")




indel = tsgExcess %>% select(gene = gene_name, excess_ind)
indel = left_join(indel, tsgSummary, by = "gene" ) %>% filter(impact == "INDEL")
indel$impact <- NULL
indel$excess_ind <- round(indel$excess_ind)
indel$Unexplained <-  indel$excess_ind - indel$Biallelic - indel$MultiHit


tsgCounts = tsgSummary %>% filter(impact == "Essential_Splice")
excessSummary = tsgExcess %>% select(gene = gene_name, Excess = excess_spl)

#jon = missense %>% group_by(gene) %>% select(-excess_mis, -Redundant, -SingleHit) %>% gather("type", "N", 2:4)
#jon$N = jon$N / jon$excess_mis




do_plot_tsg(missense[1:3,], "Missense")
do_plot_tsg(missense[4:123,], "Missense")

do_plot_tsg(nonsense[1:2,], "Nonsense")
do_plot_tsg(nonsense[3:131,], "Nonsense")

do_plot_tsg(splice, "Splice")
#data= missense

jon = plot_tsg_excess(tsgMissenseSummary, "Missense")


plot_tsg_excess <- function(data, title) {
  data = data %>% group_by(gene) %>% gather("type", "N", Biallelic, MultiHit, Unexplained)
  data = data %>% unite(gene, gene, codons)
  
  dataGeneLevels = data %>% group_by(gene) %>% summarise(N = sum(N)) %>% arrange(-N)
  data$gene <- factor(data$gene, dataGeneLevels$gene)
  
  p1 = ggplot(data=data, aes(gene, N)) +
    geom_bar(aes(fill = type), stat = "identity") + 
    ggtitle(paste0("TSG - ", title)) + xlab("Genes") + ylab("Variants") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  
  data2 = data %>% select(gene, InsufficientSingles, ExcessSingles, RequiredSingles) %>% gather("type", "N", InsufficientSingles, ExcessSingles, RequiredSingles)
  data2 = data2 %>% unite(gene, gene, codons)
  
  data2$gene <- factor(data2$gene, jonGeneLevels$gene)
  
  p2 = ggplot(data=data2, aes(gene, N)) +
    geom_bar(aes(fill = type), stat = "identity") + 
    ggtitle(paste0("TSG - Singles - ", title)) + xlab("Genes") + ylab("Variants") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  
  return(list(p1, p2))
}

data = oncoMissense
title = "Missense"
onco1 = do_plot_onco(oncoMissense[1:3, ], "Missense")
onco2 = do_plot_onco(oncoMissense[4:nrow(oncoMissense), ], "Missense")

multiplot(onco1[[1]], onco1[[2]], onco2[[1]], onco2[[2]], cols = 2)
multiplot(onco1[[1]], onco1[[2]], cols = 1)
multiplot(onco2[[1]], onco2[[2]], cols = 1)

do_plot_onco <- function(data, title) {
  data1 = data %>% group_by(gene) %>% gather("type", "N", Hotspot, NearHotspot, Unexplained)
  data1 = data1 %>% unite(gene, gene, codons)
  
  data1GeneLevels = data %>% group_by(gene) %>% summarise(N = sum(N)) %>% arrange(-N)
  data1$gene <- factor(data1$gene, data1GeneLevels$gene)
  
  
  p1 = ggplot(data1=jon, aes(gene, N)) +
    geom_bar(aes(fill = type), stat = "identity") + 
    ggtitle(paste0("Onco - ", title)) + xlab("Genes") + ylab("Variants") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

  data2 = data %>% select(gene, InsufficientHits, ExcessHits, RequiredHits) %>% gather("type", "N", InsufficientHits, ExcessHits, RequiredHits)
  data2 = left_join(data2, canonicalTranscripts[, c("gene", "codons")], by = "gene")
  data2 = data2 %>% unite(gene, gene, codons)
  
  data2$gene <- factor(data2$gene, jonGeneLevels$gene)
  
  p2 = ggplot(data=data2, aes(gene, N)) +
    geom_bar(aes(fill = type), stat = "identity") + 
    ggtitle(paste0("Onco - Hits - ", title)) + xlab("Genes") + ylab("Variants") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  
  return (multiplot(p1, p2))
}

