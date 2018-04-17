library(dplyr)
library(tidyr)

oncoFile <- read.csv("~/data/r/oncoKb.tsv", sep="\t", stringsAsFactors=FALSE)[,c("Gene", "Alteration", "Oncogenicity")]
civicFile <- read.csv("~/data/r/civic_variants.tsv", sep="\t", stringsAsFactors=FALSE)[,c("gene", "variant", "variant_types")]
cgiFile <- read.csv("~/data/r/cgi_biomarkers_per_variant.tsv", sep="\t", stringsAsFactors=FALSE)#[,c("Gene", "Alteration", "Alteration.type")]

extractType <- function(string){
  if(grepl("amp", string) | grepl("AMPLIFICATION", string))
    return("Amplification")
  else return("Deletion")
}

oncoAmpsDels <- oncoFile %>% filter(Alteration == "Amplification" | Alteration == "Deletion") %>% mutate(gene = Gene, type = Alteration) %>% select(gene, type) %>% 
  mutate(source = "oncoKb") %>% distinct()
civicAmpsDels <- civicFile %>% filter(variant == "AMPLIFICATION" | variant == "DELETION" | variant == "LOH") %>% rowwise() %>% mutate(type = extractType(variant)) %>%
  select(gene, type) %>% mutate(source = "civic") %>% distinct()
cgiAmpsDels <- cgiFile %>% filter(Alteration.type == "CNA") %>% rowwise() %>% mutate (type = extractType(Alteration)) %>% mutate(gene = Gene) %>% select(gene, type) %>%
  mutate(source = "cgi") %>% distinct()

ampsDels = rbind(oncoAmpsDels, civicAmpsDels, cgiAmpsDels) %>% mutate(value = TRUE) %>% spread(source, value)

write.csv(ampsDels, "~/data/r/ampsDels.csv", row.names = FALSE)