library(dplyr)
library(stringr)
library(tidyr)

oncoFile <- read.csv("~/data/r/oncoKb.tsv", sep="\t", stringsAsFactors=FALSE)[,c("Gene", "Alteration", "Oncogenicity")]
civicFile <- read.csv("~/data/r/civic_variants.tsv", sep="\t", stringsAsFactors=FALSE)[,c("gene", "variant", "variant_types")]
cgiFile <- read.csv("~/data/r/cgi_biomarkers_per_variant.tsv", sep="\t", stringsAsFactors=FALSE)[,c("Gene", "Alteration", "Alteration.type")]
cosmicFile <- read.csv("~/data/r/cosmic_gene_fusions.csv", stringsAsFactors=FALSE)[,c(1, 2)]
colnames(cosmicFile) <- c("H_gene", "T_gene")

# --- functions ---

genePattern = "([A-Za-z0-9-]+)"
separatorPattern = "(?:\ ?[-\\?]\ ?)"

extractGene <- function(string){
  if(grepl("_", string))
    return(str_match(string, paste(genePattern, "_", sep=""))[1,2])
  else return(string)
}

isHGene <- function(gene, fusion, separator){
  return(!grepl(paste(separator, geneStartLetters(gene), sep=""), fusion))
}

geneStartLetters <- function(gene){
  return(substr(gene, 0, 3))
}

hGene <- function(gene, fusion, separator){
  if(isHGene(gene, fusion, separator)){
    return(gene)
  } else{
    pattern <- paste(genePattern, separator, geneStartLetters(gene), sep = "" )
    tGene <- str_match(fusion, pattern)
    return(tGene[1,2])
  }
}

tGene <- function(gene, fusion, separator){
  if(!isHGene(gene, fusion, separator)){
    return(gene)
  } else{
    pattern <- paste(gene, separator, genePattern, sep = "" )
    tGene <- str_match(fusion, pattern)
    return(tGene[1,2])
  }
}

extractColumn <- function(index, columns) {
  sortedColumns <- columns %>% sort(decreasing = TRUE)
  columnValue <- sortedColumns[index]
  return(columnValue)
}

extractColumnName <- function(index, columns) {
  sortedColumns <- columns %>% sort(decreasing = TRUE)
  return(names(sortedColumns)[index])
}

flipGenePair <- function(data, first, second){
  rows <- data %>% filter(H_gene == first & T_gene == second) %>%  mutate(H_gene = second, T_gene = first)
  return(data %>% filter(H_gene != first | T_gene != second) %>% rbind(rows))
}

removePair <- function(data, first, second){
  return(data %>% filter(H_gene != first | T_gene != second))
}

correctCivicVariant <- function(data, geneString, variantString, newVariantString){
  rows <- data %>% filter(gene == geneString & variant == variantString) %>%  mutate(variant = newVariantString)
  return(data %>% filter(gene != geneString | variant != variantString) %>% rbind(rows))
}

# -----------------
civicFile <- civicFile %>% correctCivicVariant("KMT2A", "MLL-MLLT3", "KMT2A-MLLT3")

oncoPairs <- oncoFile %>% filter(grepl("Fusion$", Alteration)) %>% rowwise() %>%
  mutate(H_gene = hGene(Gene, Alteration, separatorPattern), T_gene = tGene(Gene, Alteration, separatorPattern), Source ="oncoKb") %>%
  select(H_gene, T_gene, Source) %>% flipGenePair("ROS1", "CD74") %>% flipGenePair("EP300", "MLL") %>% flipGenePair("EP300", "MOZ") %>% flipGenePair("RET", "CCDC6") %>% distinct

oncoPromiscuous <- oncoFile %>% filter(grepl("Fusions", Alteration)) %>% mutate(H_gene = Gene, T_gene = NA, Source ="oncoKb") %>%
  select(H_gene, T_gene, Source) %>% distinct

civicFusions <- civicFile %>% filter(grepl("fusion", variant_types)) %>% rowwise() %>%
  mutate(H_gene = hGene(gene, variant, separatorPattern), T_gene = tGene(gene, variant,separatorPattern), Source = "civic") %>%
  select(H_gene, T_gene, Source) %>% distinct

civicPairs <- civicFusions %>% filter(!is.na(H_gene) && !is.na(T_gene)) %>% removePair("BRAF", "CUL1")
civicPromiscuous <- civicFusions %>% filter(is.na(H_gene) || is.na(T_gene))

cgiFusions <- cgiFile %>% filter(Alteration.type=="FUS") %>% rowwise() %>%
  mutate(H_gene = hGene(Gene, Alteration, "__"), T_gene = tGene(Gene, Alteration, "__"), Source = "cgi") %>%
  select(H_gene, T_gene, Source) %>% distinct

cgiPairs <- cgiFusions %>% filter(!is.na(H_gene) && !is.na(T_gene)) %>% 
  flipGenePair("ABL1", "BCR") %>% flipGenePair("PDGFRA", "FIP1L1") %>% flipGenePair("PDGFB", "COL1A1") %>% removePair("RET", "TPCN1") %>% distinct

cgiPromiscuous <- cgiFusions %>% filter(is.na(H_gene) || is.na(T_gene))

cosmicPairs <- cosmicFile %>% rowwise() %>% mutate(H_gene = extractGene(H_gene), T_gene = extractGene(T_gene), Source = "cosmic") %>% distinct

allFusionPairs <- rbind(cosmicPairs, oncoPairs, civicPairs, cgiPairs)
allExternalPromiscuous <- rbind(cgiPromiscuous, oncoPromiscuous, civicPromiscuous)
allDistinctFusionPairs <- allFusionPairs %>% select(H_gene, T_gene) %>% distinct
allDistinctExternalPromiscuous <- allExternalPromiscuous[,c("H_gene", "T_gene")] %>% distinct
allDistinctFusionPairsWide <- allFusionPairs %>% mutate(value = TRUE) %>% spread(Source, value)
# all fusions in knowledgebases
allFusionsWide <- rbind(allFusionPairs, allExternalPromiscuous) %>% mutate(value = TRUE) %>% spread(Source, value)

promiscuousHFusions <- allDistinctFusionPairs %>% ungroup() %>% group_by(H_gene) %>% count %>% filter(n > 2) %>%
  mutate(T_gene = NA_character_, Source = "promiscuous_H") %>% select(H_gene, T_gene, Source)
promiscuousTFusions <- allDistinctFusionPairs %>% ungroup() %>% group_by(T_gene) %>% count %>% filter(n > 2) %>%
  mutate(H_gene = NA_character_, Source = "promiscuous_T") %>% select(H_gene, T_gene, Source)
verifiedExternalPromiscuousH <- allDistinctExternalPromiscuous %>% merge(allDistinctFusionPairs, by="H_gene") %>% group_by(H_gene) %>% select(H_gene) %>%
  mutate(T_gene = NA_character_, Source = "external") %>% distinct
verifiedExternalPromiscuousT <- allDistinctExternalPromiscuous %>% mutate(T_gene = H_gene, H_gene = NA_character_) %>% merge(allDistinctFusionPairs, by="T_gene") %>%
  group_by(T_gene) %>% select(T_gene) %>% mutate(H_gene = NA_character_, Source = "external") %>% distinct
allPromiscuousFusions <- rbind(verifiedExternalPromiscuousH, verifiedExternalPromiscuousT, promiscuousHFusions, promiscuousTFusions) %>% ungroup()

# all fusion pairs + verified promiscuous + inferred promiscuous
allKnownFusionsWide <- rbind(allFusionPairs, allPromiscuousFusions) %>% mutate(value = TRUE) %>% spread(Source, value) %>% arrange(H_gene, T_gene)
allKnownFusions <- allKnownFusionsWide %>% select(H_gene, T_gene)
knownPromiscuousH <- allKnownFusions %>% filter(is.na(T_gene)) %>% select(H_gene)
knownPromiscuousT <- allKnownFusions %>% filter(is.na(H_gene)) %>% select(T_gene)
knownFusionPairs <- allKnownFusions %>% filter(!is.na(H_gene) & !is.na(T_gene))
write.csv(allKnownFusions, "~/data/r/allKnownFusions.csv", row.names = FALSE)
write.csv(allKnownFusionsWide, "~/data/r/allKnownFusionsWide.csv", row.names = FALSE)
write.csv(knownPromiscuousH, "~/data/r/knownPromiscuousH.csv", row.names = FALSE)
write.csv(knownPromiscuousT, "~/data/r/knownPromiscuousT.csv", row.names = FALSE)
write.csv(knownFusionPairs, "~/data/r/knownFusionPairs.csv", row.names = FALSE)
