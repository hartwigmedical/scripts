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

allFusionPairs <- rbind(cosmicPairs, oncoPairs, civicPairs, cgiPairs) %>% mutate(value = TRUE) %>% spread(Source, value)
allExternalPromiscuous <- rbind(cgiPromiscuous, oncoPromiscuous, civicPromiscuous) %>% mutate(value = TRUE) %>% spread(Source, value) %>% mutate(cosmic = NA)
# all fusions in knowledgebases
allFusions <- rbind(allFusionPairs, allExternalPromiscuous)

promiscuousHFusions <- allFusionPairs %>% group_by(H_gene) %>%
  summarise(oncoKb = any(oncoKb), civic = any(civic), cgi = any(cgi), cosmic = any(cosmic), count = n()) %>%
  filter(count > 2) %>% mutate(T_gene = NA_character_) %>% select(H_gene, T_gene, oncoKb, civic, cgi, cosmic)

verifiedExternalPromiscuousH <- allExternalPromiscuous %>% merge(allFusionPairs, by="H_gene") %>% group_by(H_gene) %>%
  summarise(oncoKb = any(oncoKb.x | oncoKb.y), civic = any(civic.x | civic.y), cgi = any(cgi.x | cgi.y), cosmic = any(cosmic.x | cosmic.y)) %>%
  mutate(T_gene = NA_character_)

knownPromiscuousH <- rbind(promiscuousHFusions, verifiedExternalPromiscuousH) %>% group_by(H_gene) %>% 
  summarise(cgi = any(cgi), civic = any(civic), cosmic = any(cosmic), oncoKb = any(oncoKb)) %>% mutate(gene = H_gene) %>%
  select(gene, cgi, civic, cosmic, oncoKb)

promiscuousTFusions <- allFusionPairs %>% group_by(T_gene) %>%
  summarise(oncoKb = any(oncoKb), civic = any(civic), cgi = any(cgi), cosmic = any(cosmic), count = n()) %>%
  filter(count > 2) %>% mutate(H_gene = NA_character_) %>% select(H_gene, T_gene, oncoKb, civic, cgi, cosmic)

verifiedExternalPromiscuousT <- allExternalPromiscuous %>% mutate(T_gene = H_gene, H_gene = NA_character_) %>% 
  merge(allFusionPairs, by="T_gene") %>% group_by(T_gene) %>% 
  summarise(oncoKb = any(oncoKb.x | oncoKb.y), civic = any(civic.x | civic.y), cgi = any(cgi.x | cgi.y), cosmic = any(cosmic.x | cosmic.y)) %>% 
  mutate(H_gene = NA_character_)

knownPromiscuousT <- rbind(promiscuousTFusions, verifiedExternalPromiscuousT) %>% group_by(T_gene) %>%
  summarise(cgi = any(cgi), civic = any(civic), cosmic = any(cosmic), oncoKb = any(oncoKb)) %>% mutate(gene = T_gene) %>%
  select(gene, cgi, civic, cosmic, oncoKb)

allPromiscuousFusions <- rbind(knownPromiscuousH %>% mutate(T_gene = NA_character_, H_gene = gene) %>% select(-gene),
                               knownPromiscuousT %>% mutate(H_gene = NA_character_, T_gene = gene) %>% select(-gene))  %>% ungroup()

# all fusion pairs + verified promiscuous + inferred promiscuous
allKnownFusions <- rbind(allFusionPairs, allPromiscuousFusions) %>% arrange(H_gene, T_gene)
write.csv(allKnownFusions, "~/data/r/allKnownFusions.csv", row.names = FALSE)

knownPromiscuousFive <- knownPromiscuousH %>% mutate(OncoKB = oncoKb, COSMIC= cosmic, CGI = cgi, CIViC = civic) %>% select(gene, OncoKB, COSMIC, CGI, CIViC)
write.csv(knownPromiscuousFive, "~/data/r/knownPromiscuousFive.csv", row.names = FALSE)
knownPromiscuousThree <- knownPromiscuousT %>% mutate(OncoKB = oncoKb, COSMIC= cosmic, CGI = cgi, CIViC = civic) %>% select(gene, OncoKB, COSMIC, CGI, CIViC)
write.csv(knownPromiscuousThree, "~/data/r/knownPromiscuousThree.csv", row.names = FALSE)
knownFusionPairs <- allKnownFusions %>% filter(!is.na(H_gene) & !is.na(T_gene)) %>% mutate(OncoKB = oncoKb, COSMIC= cosmic, CGI = cgi, CIViC = civic) %>%
  select(H_gene, T_gene, OncoKB, COSMIC, CGI, CIViC)
write.csv(knownFusionPairs, "~/data/r/knownFusionPairs.csv", row.names = FALSE)
